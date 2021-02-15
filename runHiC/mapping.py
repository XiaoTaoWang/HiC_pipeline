# Created on Fri Sep 7 19:15:32 2018

# Author: XiaoTao Wang

import os, subprocess, atexit
from runHiC.utilities import cleanFile, sleep
from runHiC.filtering import create_frag, stats_pairs, stats_samfrag
from runHiC.quality import outStatsCache

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

def buildMapIndex(aligner, genomeFolder, genomeName):
    """
    Build bwa/chromap/minimap2 index files.
    """
    lockFile = os.path.join(genomeFolder, '.'.join([genomeName, aligner, 'lock']))
    # Aquire lock
    lock = open(lockFile, 'wb')
    lock.close()

    atexit.register(cleanFile, lockFile)

    wholeGenome = os.path.join(genomeFolder, '.'.join([genomeName, 'fa']))

    if aligner=='minimap2':
        indexOut = os.path.join(genomeFolder, '.'.join([genomeName, 'minimap2', 'mmi']))
        build_command = ['minimap2', '-x', 'sr', '-d', indexOut, wholeGenome]
    elif aligner=='chromap':
        indexOut = os.path.join(genomeFolder, '.'.join([genomeName, 'chromap', 'mmi']))
        build_command = ['chromap', '-i', '-r', wholeGenome, '-o', indexOut]
    else:
        build_command = ['bwa', 'index', wholeGenome]
    
    subprocess.check_call(' '.join(build_command), shell=True)
    

def map_core(fastq_1, fastq_2, ref_fa, ref_index, outdir, aligner='minimap2', min_mapq=1, nthread=1):

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
        map_command = ['chromap', '-m', '-r', ref_fa, '-x', ref_index, '-t', str(nthread), '-1', fastq_1, '-2', fastq_2,
                       '-o', outpath, '--pairs', '--split-alignment', '-e', str(8), '-f', '500,1000', '-q', str(min_mapq)]
        bam_command = []
    else:
        if aligner=='minimap2':
            map_command = ['minimap2', '-ax', 'sr', '-t', str(nthread), ref_index, fastq_1, fastq_2]
        else:
            map_command = ['bwa', 'mem', '-SP5M', '-t', str(nthread), ref_index, fastq_1, fastq_2]
            
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
        stats = _collect_chromap_stats(chromap_stderr)
    else:
        stats = {}
    
    return outpath, stats

def _collect_chromap_stats(chromap_stderr):

    from collections import defaultdict

    stats = defaultdict(int)
    for line in chromap_stderr:
        tmp = line.decode().rstrip()
        if tmp.startswith('Number of reads:'):
            stats['000_SequencedReads'] = int(tmp.split(':')[1].strip().rstrip('.')) // 2
        if tmp.startswith('Number of output mappings (passed filters):'):
            stats['010_DoubleSideMappedReads'] = int(tmp.split(':')[1].strip().rstrip('.'))
            stats['100_NormalPairs'] = stats['010_DoubleSideMappedReads']
    
    return stats

def parse_align(align_path, align_stats, outfile, genomepath, chromsizes, assembly, min_mapq, max_molecule_size,
              max_inter_align_gap, walks_policy, include_readid, include_sam, drop_seq, tmpdir, enzyme, nproc_in,
              nproc_out, memory):
    
    frag_path = create_frag(genomepath, chromsizes, enzyme, tmpdir)
    out_total = outfile.replace('.pairsam.gz', '.total.pairsam.gz')
    
    #### step 1
    if align_path.endswith('.pairs'):
        basic_command = ['pairtools', 'flip', '-c', chromsizes, '--nproc-in', str(nproc_in),
                         '--nproc-out', str(nproc_out), align_path]
    else:
        basic_command = ['pairtools', 'parse', '-c', chromsizes, '--assembly', assembly,
                        '--min-mapq', str(min_mapq), '--max-molecule-size', str(max_molecule_size),
                        '--max-inter-align-gap', str(max_inter_align_gap), '--walks-policy', walks_policy,
                        '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out)]
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

    stats_pool = {'pseudo': stats}
    stats_pre = outfile.replace('.pairsam.gz', '.pstats') # pickled stats
    outStatsCache(stats_pool, stats_pre)