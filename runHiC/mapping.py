# Created on Fri Sep 7 19:15:32 2018

# Author: XiaoTao Wang

import os, subprocess, io, pairtools
from runHiC.utilities import cleanFile, extract_chrom_sizes
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

#### bam check
def validate_bam(bam_path, full_scan=True):
    """
    Validate BAM file:
    1. File exists and is not empty.
    2. samtools quickcheck checks BAM header and EOF.
    3. samtools view -c reads the entire BAM to detect
       internal corruption or truncated records.
    """

    if not os.path.isfile(bam_path):
        raise RuntimeError(
            "BAM validation failed: file does not exist: {}".format(bam_path)
        )

    if os.path.getsize(bam_path) == 0:
        raise RuntimeError(
            "BAM validation failed: empty BAM file: {}".format(bam_path)
        )

    checks = [
        (
            "samtools quickcheck",
            ["samtools", "quickcheck", "-v", bam_path]
        )
    ]

    if full_scan:
        checks.append(
            (
                "samtools full scan",
                ["samtools", "view", "-c", bam_path]
            )
        )

    for label, command in checks:

        proc = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        stdout, stderr = proc.communicate()

        if proc.returncode != 0:

            error_message = stderr.decode(
                "utf-8", "replace"
            ).strip()

            if not error_message:
                error_message = stdout.decode(
                    "utf-8", "replace"
                ).strip()

            raise RuntimeError(
                "{} failed for BAM: {}\n{}".format(
                    label, bam_path, error_message
                )
            )

        # Make sure samtools view -c returned a valid count.
        if label == "samtools full scan":
            try:
                count = int(stdout.decode("utf-8").strip())
            except (ValueError, AttributeError):
                raise RuntimeError(
                    "BAM validation failed: "
                    "invalid samtools view output: {}".format(bam_path)
                )

            if count <= 0:
                raise RuntimeError(f"BAM validation failed: no alignment records: {bam_path}")

    return True

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
    
    if aligner == 'chromap':
        outlog = os.path.join(tmpdir, os.path.basename(outpath).replace(outformat, '.chromap.log'))
        if os.path.exists(outpath):
            os.remove(outpath)
        try:
            # Write stderr to disk while Chromap runs; a PIPE can fill before wait().
            with open(outlog, 'wb') as stderr:
                result = subprocess.run(map_command, stderr=stderr)
            if result.returncode:
                raise RuntimeError(f"Chromap failed (exit {result.returncode}); see {outlog}")
            if not os.path.isfile(outpath) or os.path.getsize(outpath) == 0:
                raise RuntimeError(f"Chromap produced no pairs file; see {outlog}")
            stats = _collect_chromap_stats(outlog)
        except BaseException:
            if os.path.exists(outpath):
                os.remove(outpath)
            raise
    else:
        # Remove an old BAM to prevent accidental reuse.
        if os.path.exists(outpath):
            os.remove(outpath)

        align_proc = None
        bam_proc = None

        try:

            with open(outpath, 'wb') as bam_out:

                # Step 1: Run BWA/BWA-MEM2.
                align_proc = subprocess.Popen(
                    map_command,
                    stdout=subprocess.PIPE,
                    bufsize=-1
                )

                # Step 2: Convert SAM stream to BAM.
                try:
                    bam_proc = subprocess.Popen(
                        bam_command,
                        stdin=align_proc.stdout,
                        stdout=bam_out,
                        bufsize=-1
                    )

                finally:
                    # Important: close parent's copy of the pipe.
                    align_proc.stdout.close()

                # Check both processes.
                bam_returncode = bam_proc.wait()
                align_returncode = align_proc.wait()

            # Both processes must finish successfully.
            if align_returncode != 0 or bam_returncode != 0:

                raise RuntimeError(
                    "Mapping-to-BAM pipeline failed!\n"
                    "Aligner exit code: {}\n"
                    "Samtools exit code: {}\n"
                    "BAM: {}".format(
                        align_returncode,
                        bam_returncode,
                        outpath
                    )
                )

            # Step 3: Validate the complete BAM.
            print("[BAM QC] Validating BAM: {}".format(outpath))

            validate_bam(outpath, full_scan=True)

            print("[BAM QC] BAM validation PASSED: {}".format(outpath))

        except BaseException:

            # Stop any remaining subprocesses.
            for process in (bam_proc, align_proc):

                if process is not None and process.poll() is None:
                    try:
                        process.terminate()
                    except OSError:
                        pass

            for process in (bam_proc, align_proc):

                if process is not None:
                    try:
                        process.wait()
                    except Exception:
                        pass

            # Delete incomplete or invalid BAM.
            if os.path.exists(outpath):
                os.remove(outpath)

            print(
                "[BAM QC] FAILED! Invalid BAM removed: {}".format(
                    outpath
                )
            )

            raise

        stats = {}
    
    return outpath, stats

def _collect_chromap_stats(outlog):
    from collections import defaultdict

    stats = defaultdict(int)
    with open(outlog, encoding='utf-8', errors='replace') as stream:
        for line in stream:
            if line.startswith('Number of reads:'):
                stats['000_SequencedReads'] = int(line.split(':')[1].strip().rstrip('.')) // 2
            if line.startswith('Number of output mappings (passed filters):'):
                stats['010_DoubleSideMappedReads'] = int(line.split(':')[1].strip().rstrip('.'))
                stats['100_NormalPairs'] = stats['010_DoubleSideMappedReads']
    return stats


def _run_checked_pipeline(*commands):
    """Abort if any stage fails, even when its downstream consumer exits zero."""
    processes = []
    try:
        for i, command in enumerate(commands):
            upstream = processes[-1].stdout if processes else None
            try:
                processes.append(subprocess.Popen(
                    command, stdin=upstream,
                    stdout=subprocess.PIPE if i < len(commands) - 1 else None,
                ))
            finally:
                if upstream is not None:
                    upstream.close()
        for process, command in reversed(list(zip(processes, commands))):
            code = process.wait()
            if code:
                raise subprocess.CalledProcessError(code, command)
    finally:
        for process in processes:
            if process.stdout is not None:
                process.stdout.close()
            if process.poll() is None:
                process.terminate()
        for process in processes:
            process.wait()


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
        if process.returncode:
            raise subprocess.CalledProcessError(process.returncode, command)
    
    instream.close()
    outstream.close()

def parse_align(align_path, align_stats, outfile, genomepath, phased_snp, chromsizes, assembly, min_mapq,
                max_molecule_size, max_inter_align_gap, walks_policy, include_readid, include_sam, drop_seq,
                tmpdir, enzyme, memory, add_frag):

    # Validate BAM before starting pairtools.
    if not align_path.endswith('.pairs'):

        print("[BAM QC] Checking BAM before pair parsing...")

        validate_bam(align_path, full_scan=False)

        print("[BAM QC] BAM passed pre-parsing check.")

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
        
            _run_checked_pipeline(
                basic_command,
                ['pairtools', 'sort', '-o', out_total, '--nproc', str(nproc_out), '--memory', memory,
                 '--tmpdir', tmpdir, '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out)],
            )
        else:
            basic_command = ['pairtools', 'parse', '-c', chromsizes, '--assembly', assembly,
                             '--min-mapq', str(min_mapq), '--max-molecule-size', str(max_molecule_size),
                             '--max-inter-align-gap', str(max_inter_align_gap), '--add-columns', 'mapq',
                             '--walks-policy', walks_policy, '--nproc-in', str(nproc_in),
                             '--nproc-out', str(nproc_out), align_path]
            phase_command = [
                'runHiC-phase', '--phased-SNPs', phased_snp, '--nproc-in', str(nproc_in),
                '--nproc-out', str(nproc_out), '--max-molecule-size', str(max_molecule_size),
            ]
            if include_readid:
                phase_command.append('--include-readid')
            if include_sam:
                phase_command.append('--include-sam')
            _run_checked_pipeline(
                basic_command, phase_command,
                ['pairtools', 'sort', '-o', out_total, '--nproc', str(nproc_out), '--memory', memory,
                 '--tmpdir', tmpdir, '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out)],
            )

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
        subprocess.run(
            ['pairtools', 'select', '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out),
             '-o', outpath_1, '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu")',
             out_total],
            check=True,
        )
        os.remove(out_total)

    #### step 3
    if add_frag:
        frag_path = create_frag(genomepath, chromsizes, enzyme, tmpdir)
        outpath_2 = outfile.replace('.pairsam.gz', '.select.samefrag.pairsam.gz')
        _run_checked_pipeline(
            ['pairtools', 'restrict', '-f', frag_path, '--nproc-in', str(nproc_in),
             '--nproc-out', str(nproc_out), outpath_1],
            ['pairtools', 'select', '--nproc-in', str(nproc_in), '--nproc-out', str(nproc_out),
             '--output-rest', outfile, '-o', outpath_2,
             '(COLS[-6]==COLS[-3]) and (chrom1==chrom2)'],
        )
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

    if not os.path.isfile(outfile) or os.path.getsize(outfile) == 0:
        raise RuntimeError(f"Pair parsing did not produce a nonempty output file: {outfile}")

    stats_pool = {'pseudo': stats}
    stats_pre = outfile.replace('.pairsam.gz', '.pstats') # pickled stats
    outStatsCache(stats_pool, stats_pre)