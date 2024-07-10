# Created on Fri Sep 7 19:15:32 2018

# Author: XiaoTao Wang

import os, subprocess, io, pairtools
from runHiC.utilities import cleanFile, sleep, extract_chrom_sizes
from runHiC.filtering import create_frag, stats_pairs, stats_samfrag
from runHiC.quality import outStatsCache

if pairtools.__version__.startswith('0'):
    from pairtools import _fileio, _headerops
else:
    from pairtools.lib import fileio as _fileio
    from pairtools.lib import headerops as _headerops

##### functions handling with sequencing reads
def commandExists(command):
    "Check if the bash command exists"
    command = command.split()[0]
    if subprocess.call(['which', command]) != 0:
        return False
    return True

def gzipWriter(filename):
    """
    Create a writing process with gzip or parallel gzip (pigz) attached to
    it.
    
    """
    filename = os.path.abspath(filename)
    
    with open(filename, 'wb') as outFile:
        if commandExists("pigz"):
            writer = ["pigz", "-c", "-4"]
        else:
            writer = ["gzip", "-c", "-1"]

        pwrite = subprocess.Popen(writer, stdin = subprocess.PIPE,
                                  stdout = outFile, shell = False,
                                  bufsize = -1)
    return pwrite

def uncompressSRA(filename, folder):
    
    inFile = os.path.abspath(filename)
    
    pread = subprocess.Popen(['fastq-dump', inFile, "-Z", "--split-files"],
                              stdout = subprocess.PIPE, bufsize = -1)
    
    inStream = pread.stdout
    
    outFile = os.path.split(filename)[1].replace('.sra', '') + '_{0}.fastq.gz'
    outProc1 = gzipWriter(os.path.join(folder, outFile).format(1))
    outProc2 = gzipWriter(os.path.join(folder, outFile).format(2))
    outStream1 = outProc1.stdin
    outStream2 = outProc2.stdin
    
    count = 0
    while True:
        line = inStream.readline()
        textline = line.decode()
        
        try:
            assert textline[0] == '@'
        except AssertionError:
            raise IOError('{0} is an invalid fastq file'.format(filename))
        except IndexError:
            break
        
        fastq_entry_1 = (line, inStream.readline(),
                         inStream.readline(), inStream.readline())
        outStream1.writelines(fastq_entry_1)
        
        fastq_entry_2 = (inStream.readline(), inStream.readline(),
                         inStream.readline(), inStream.readline())
        outStream2.writelines(fastq_entry_2)
        
        count += 1
    
    outProc1.communicate()
    outProc2.communicate()

def splitSRA(filename, folder, splitBy = 4000000):

    inFile = os.path.abspath(filename)
    outFile = os.path.split(inFile)[1].replace('.sra', '') + '_chunk{0}_{1}.fastq.gz'
    pread = subprocess.Popen(['fastq-dump', inFile, "-Z", "--split-files"],
                              stdout = subprocess.PIPE, bufsize = -1)
    inStream = pread.stdout

    halted = False
    counters = []
    for counter in range(1000000):

        outProc1 = gzipWriter(os.path.join(folder, outFile).format(counter, 1))
        outProc2 = gzipWriter(os.path.join(folder, outFile).format(counter, 2))
        outStream1 = outProc1.stdin
        outStream2 = outProc2.stdin

        for j in range(splitBy):

            line = inStream.readline()
            textline = line.decode()

            try:
                assert textline[0] == '@'
            except AssertionError:
                raise IOError('{0} is an invalid fastq file'.format(filename))
            except IndexError:
                halted = True
                counters.append(j)
                break

            fastq_entry = (line, inStream.readline(),
                           inStream.readline(), inStream.readline())

            outStream1.writelines(fastq_entry)
            outStream2.writelines((inStream.readline(), inStream.readline(),
                       inStream.readline(), inStream.readline()))

        outProc1.communicate()
        outProc2.communicate()
        if halted:
            return counters
        counters.append(splitBy)
    return counters


def splitSingleFastq(filename, pre, pair_index, folder, splitBy=4000000):

    inFile = os.path.abspath(filename)
    outFile = pre + '_chunk{0}_{1}.fastq.gz'
    
    if inFile.endswith('.gz'):
        pread = subprocess.Popen(['gunzip', inFile, '-c'],
                                  stdout = subprocess.PIPE, bufsize = -1)
    else:
        pread = subprocess.Popen(['cat', inFile],
                                  stdout = subprocess.PIPE, bufsize = -1)
                                  
    inStream = pread.stdout

    halted = False
    counters = []
    for counter in range(1000000):

        outProc1 = gzipWriter(os.path.join(folder, outFile).format(counter, pair_index))
        outStream1 = outProc1.stdin

        for j in range(splitBy):

            line = inStream.readline()
            textline = line.decode()

            try:
                assert textline[0] == '@'
            except AssertionError:
                raise IOError('{0} is not a fastq file'.format(filename))
            except IndexError:
                halted = True
                counters.append(j)
                break

            fastq_entry = (line, inStream.readline(), inStream.readline(),
                           inStream.readline())

            outStream1.writelines(fastq_entry)

        outProc1.communicate()
        counters.append(splitBy)
        
        if halted:
            return counters
        counters.append(splitBy)
        
    return counters

##### functions that invoke different read aligners
def buildMapIndex(aligner, genomeFolder, genomeName):
    """
    Build bwa-mem/bwa-mem2/chromap index files.
    """
    lockFile = os.path.join(genomeFolder, '.'.join([genomeName, aligner, 'lock']))
    # Aquire lock
    lock = open(lockFile, 'wb')
    lock.close()
    
    wholeGenome = os.path.join(genomeFolder, '.'.join([genomeName, 'fa']))

    if aligner=='chromap':
        indexOut = os.path.join(genomeFolder, '.'.join([genomeName, 'chromap-runhic', 'mmi']))
        build_command = ['chromap', '-i', '-r', wholeGenome, '-o', indexOut]
    elif aligner=='bwa-mem':
        build_command = ['bwa', 'index', wholeGenome]
    else:
        build_command = ['bwa-mem2', 'index', wholeGenome]
    
    subprocess.check_call(' '.join(build_command), shell=True)

    cleanFile(lockFile)
    
def map_core(fastq_1, fastq_2, ref_fa, ref_index, outdir, tmpdir, aligner='chromap', nthread=1):

    if aligner=='chromap':
        outformat = '.pairs'
    else:
        outformat = '.bam'

    # output file name
    if fastq_1.endswith('_1.fastq.gz'):
        outpath = os.path.join(outdir,
                        os.path.split(fastq_1)[1].replace('_1.fastq.gz', outformat))
    else:
        outpath = os.path.join(outdir,
                        os.path.split(fastq_1)[1].replace('_1.fastq', outformat))
    
    if aligner=='chromap':
        map_command = ['chromap', '--preset', 'hic', '--low-mem', '-x', ref_index,
                       '-r', ref_fa, '-t', str(nthread), '-1', fastq_1, '-2', fastq_2,
                       '-o', outpath]
        bam_command = []
    elif aligner=='bwa-mem':
        map_command = ['bwa', 'mem', '-SP', '-t', str(nthread), ref_index, fastq_1, fastq_2]    
        bam_command = ['samtools', 'view', '-bS', '-']
    else:
        map_command = ['bwa-mem2', 'mem', '-SP', '-t', str(nthread), ref_index, fastq_1, fastq_2]    
        bam_command = ['samtools', 'view', '-bS', '-']
    
    pipeline = []
    try:
        # Mapping
        if aligner=='chromap':
            pipeline.append(subprocess.Popen(map_command,
                        stderr=subprocess.PIPE,
                        bufsize=-1))
        else:
            pipeline.append(
                    subprocess.Popen(map_command,
                        stdout=subprocess.PIPE,
                        bufsize=-1))

            pipeline.append(
                    subprocess.Popen(bam_command,
                        stdin=pipeline[-1].stdout,
                        stdout=open(outpath, 'wb'),
                        bufsize=-1))

        pipeline[-1].wait()

    finally:
        sleep()
        for process in pipeline:
            if process.poll() is None:
                process.terminate()
    
    #### collect mapping statistics from stderr of chromap
    if aligner=='chromap':
        chromap_stderr = pipeline[-1].stderr
        outlog = os.path.join(tmpdir, os.path.split(outpath)[1].replace(outformat, '.chromap.log'))
        stats = _collect_chromap_stats(chromap_stderr, outlog)
    else:
        stats = {}
    
    return outpath, stats

def _collect_chromap_stats(chromap_stderr, outlog):

    from collections import defaultdict

    stats = defaultdict(int)
    with open(outlog, 'w') as out:
        for line in chromap_stderr:
            tmp = line.decode().rstrip()
            if tmp.startswith('Number of reads:'):
                stats['000_SequencedReads'] = int(tmp.split(':')[1].strip().rstrip('.')) // 2
            if tmp.startswith('Number of output mappings (passed filters):'):
                stats['010_DoubleSideMappedReads'] = int(tmp.split(':')[1].strip().rstrip('.'))
                stats['100_NormalPairs'] = stats['010_DoubleSideMappedReads']
            out.write(tmp+'\n')
    
    return stats


##### functions used for parsing/sorting/filtering original alignments
def has_correct_order(loci1, loci2, chrom_index):

    check = (chrom_index[loci1[0]], loci1[1]) <= (chrom_index[loci2[0]], loci2[1])

    return check

def _pairs_write(outstream, line, chrom_index):

    parse = line.rstrip().split()
    if not len(parse):
        return
    
    readID, c1, p1, c2, p2, strand1, strand2, pair_type = parse[:8]
    p1, p2 = int(p1), int(p2)
    loci1 = (c1, p1)
    loci2 = (c2, p2)
    if not has_correct_order(loci1, loci2, chrom_index):
        loci1, loci2 = loci2, loci1
        strand1, strand2 = strand2, strand1
    
    cols = [readID, loci1[0], str(loci1[1]), loci2[0], str(loci2[1]), strand1, strand2, pair_type]
    outstream.write('\t'.join(cols) + '\n')

def parse_chromap(align_path, out_total, chrom_path, assembly, nproc_in, nproc_out, memory, tmpdir):
    """
    Customized function to flip/sort original pairs outputed by chromap.
    """
    instream = _fileio.auto_open(align_path, mode='r', nproc=nproc_in, command=None)
    outstream = _fileio.auto_open(out_total, mode='w', nproc=nproc_out, command=None)
    chromsizes = extract_chrom_sizes(chrom_path)
    chrom_index = dict(_headerops.get_chrom_order(chrom_path))

    # writer header, chromosome sizes must be in correct order
    header = _headerops.make_standard_pairsheader(
        assembly=assembly, chromsizes=chromsizes,
        columns=['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type'],
        shape='upper triangle'
    )
    outstream.writelines((l+'\n' for l in header))
    outstream.flush()

    _, body_stream = _headerops.get_header(instream)
    # sort command
    command = r'''/bin/bash -c 'export LC_COLLATE=C; export LANG=C; sort -k 2,2 -k 4,4 -k 3,3n -k 5,5n --stable {0} {1} -S {2} {3}'''.format(
        '--parallel={0}'.format(nproc_out), '--temporary-directory={0}'.format(tmpdir), memory,'--compress-program=lz4c'
    )
    command += "'"
    with subprocess.Popen(command, stdin=subprocess.PIPE, bufsize=-1, shell=True, stdout=outstream) as process:
        stdin_wrapper = io.TextIOWrapper(process.stdin, 'utf-8')
        for line in body_stream:
            _pairs_write(stdin_wrapper, line, chrom_index)
        
        stdin_wrapper.flush()
        process.communicate()
    
    instream.close()
    outstream.close()

def parse_align(align_path, align_stats, outfile, genomepath, phased_snp, chromsizes, assembly, min_mapq,
                max_molecule_size, max_inter_align_gap, walks_policy, include_readid, include_sam, drop_seq,
                tmpdir, enzyme, memory, add_frag):
    
    out_total = outfile.replace('.pairsam.gz', '.total.pairsam.gz')
    
    nproc_in = 3
    nproc_out = 8

    #### step 1
    if align_path.endswith('.pairs'):
        parse_chromap(align_path, out_total, chromsizes, assembly, nproc_in, nproc_out, memory, tmpdir)
    else:
        if phased_snp is None:
            basic_command = ['pairtools', 'parse', '-c', chromsizes, '--assembly', assembly,
                            '--min-mapq', str(min_mapq), '--max-molecule-size', str(max_molecule_size),
                            '--max-inter-align-gap', str(max_inter_align_gap), '--add-columns', 'mapq',
                            '--walks-policy', walks_policy, '--nproc-in', str(nproc_in), '--nproc-out',
                            str(nproc_out)]
            
            if not include_readid:
                basic_command.append('--drop-readid')
            
            if not include_sam:
                basic_command.append('--drop-sam')
            
            if drop_seq:
                basic_command.append('--drop-seq')
            basic_command.append(align_path)
        
            pipeline = []
            try:
                pipeline.append(
                    subprocess.Popen(basic_command,
                        stdout=subprocess.PIPE,
                        bufsize=-1)
                )

                sort_command = ['pairtools', 'sort', '-o', out_total, '--nproc', str(nproc_out), '--memory', memory, '--tmpdir', tmpdir,
                                '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out)]
                pipeline.append(
                    subprocess.Popen(sort_command,
                        stdin=pipeline[-1].stdout,
                        stdout=None,
                        bufsize=-1)
                )

                pipeline[-1].wait()

            finally:
                sleep()
                for process in pipeline:
                    if process.poll() is None:
                        process.terminate()
        else:
            basic_command = ['pairtools', 'parse', '-c', chromsizes, '--assembly', assembly,
                            '--min-mapq', str(min_mapq), '--max-molecule-size', str(max_molecule_size),
                            '--max-inter-align-gap', str(max_inter_align_gap), '--add-columns', 'mapq',
                            '--walks-policy', walks_policy, '--nproc-in', str(nproc_in), '--nproc-out',
                            str(nproc_out), align_path]
        
            pipeline = []
            try:
                pipeline.append(
                    subprocess.Popen(basic_command,
                        stdout=subprocess.PIPE,
                        bufsize=-1)
                )

                phase_command = [
                    'runHiC-phase', '--phased-SNPs', phased_snp, '--nproc-in', str(nproc_in),
                    '--nproc-out', str(nproc_out), '--max-molecule-size', str(max_molecule_size),
                ]

                if include_readid:
                    phase_command.append('--include-readid')
                
                if include_sam:
                    phase_command.append('--include-sam')

                pipeline.append(
                    subprocess.Popen(phase_command,
                        stdin=pipeline[-1].stdout,
                        stdout=subprocess.PIPE,
                        bufsize=-1)
                )

                sort_command = ['pairtools', 'sort', '-o', out_total, '--nproc', str(nproc_out),
                                '--memory', memory, '--tmpdir', tmpdir, '--nproc-in', str(nproc_in),
                                '--nproc-out', str(nproc_out)]
                
                pipeline.append(
                    subprocess.Popen(sort_command,
                        stdin=pipeline[-1].stdout,
                        stdout=None,
                        bufsize=-1)
                )

                pipeline[-1].wait()

            finally:
                sleep()
                for process in pipeline:
                    if process.poll() is None:
                        process.terminate()
    
    # mapping stats
    if align_path.endswith('.pairs'):
        stats = align_stats
    else:
        refkey = {'total':'000_SequencedReads',
                'total_mapped':'010_DoubleSideMappedReads',
                'total_single_sided_mapped':'020_SingleSideMappedReads',
                'total_unmapped':'030_UnmappedReads'
                }
        stats = stats_pairs(out_total, refkey, nproc_in=nproc_in, nproc_out=nproc_out)
        stats['100_NormalPairs'] = stats['010_DoubleSideMappedReads']

    #### step 2
    outpath_1 = outfile.replace('.pairsam.gz', '.select.pairsam.gz')
    if align_path.endswith('.pairs'):
        mv_command = ['mv', out_total, outpath_1]
        subprocess.check_call(' '.join(mv_command), shell=True)
    else:
        pipeline = []
        try:
            select_command = ['pairtools', 'select', '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out), '-o', outpath_1,
                            '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu")',
                            out_total]
            pipeline.append(
                subprocess.Popen(select_command,
                    stdout=None,
                    bufsize=-1
                )
            )

            pipeline[-1].wait()
        
        finally:
            sleep()
            for process in pipeline:
                if process.poll() is None:
                    process.terminate()
        
        os.remove(out_total)

    #### step 3
    if add_frag:
        frag_path = create_frag(genomepath, chromsizes, enzyme, tmpdir)
        outpath_2 = outfile.replace('.pairsam.gz', '.select.samefrag.pairsam.gz')
        pipeline = []
        try:
            # assign fragment
            restrict_command = ['pairtools', 'restrict', '-f', frag_path, '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out), outpath_1]
            pipeline.append(
                subprocess.Popen(restrict_command,
                    stdout=subprocess.PIPE,
                    bufsize=-1)
            )

            ####### COLS[-6]==COLS[-3], the index may change to follow pairtools
            select_command = ['pairtools', 'select', '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out),
                            '--output-rest', outfile, '-o', outpath_2,
                            '(COLS[-6]==COLS[-3]) and (chrom1==chrom2)']
            pipeline.append(
                subprocess.Popen(select_command,
                    stdin=pipeline[-1].stdout,
                    stdout=None,
                    bufsize=-1)
            )

            pipeline[-1].wait()
        finally:
            sleep()
            for process in pipeline:
                if process.poll() is None:
                    process.terminate()
        
        os.remove(outpath_1)

        substats, libsize = stats_samfrag(outpath_2)
        stats['110_AfterFilteringReads'] = stats['100_NormalPairs'] - substats['120_SameFragmentReads']
        stats['400_TotalContacts'] = stats['110_AfterFilteringReads']
        stats.update(substats)
        stats['libsize'] = libsize
    else:
        mv_command = ['mv', outpath_1, outfile]
        subprocess.check_call(' '.join(mv_command), shell=True)
        stats['110_AfterFilteringReads'] = stats['100_NormalPairs']
        stats['400_TotalContacts'] = stats['110_AfterFilteringReads']

    stats_pool = {'pseudo': stats}
    stats_pre = outfile.replace('.pairsam.gz', '.pstats') # pickled stats
    outStatsCache(stats_pool, stats_pre)