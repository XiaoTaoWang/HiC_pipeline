# Created on Fri Sep 7 19:15:32 2018

# Author: XiaoTao Wang

import os, subprocess, atexit
from runHiC.utilities import cleanFile, sleep

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


def splitSingleFastq(filename, folder, splitBy = 4000000):

    inFile = os.path.abspath(filename)
    
    parse = os.path.split(inFile)[1].split('.')[0].split('_')
    outFile = '_'.join(parse[:-1]) + '_chunk{0}_{1}.fastq.gz'
    
    if inFile.endswith('.fastq.gz'):
        pread = subprocess.Popen(['gunzip', inFile, '-c'],
                                  stdout = subprocess.PIPE, bufsize = -1)
    else:
        pread = subprocess.Popen(['cat', inFile],
                                  stdout = subprocess.PIPE, bufsize = -1)
                                  
    inStream = pread.stdout

    halted = False
    counters = []
    for counter in range(1000000):

        outProc1 = gzipWriter(os.path.join(folder, outFile).format(counter, parse[-1]))
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
    Build bwa/minimap2 index files.
    """
    lockFile = os.path.join(genomeFolder, '.'.join([genomeName, aligner, 'lock']))
    # Aquire lock
    lock = open(lockFile, 'wb')
    lock.close()

    atexit.register(cleanFile, lockFile)

    wholeGenome = os.path.join(genomeFolder, '.'.join([genomeName, 'fa']))

    if aligner=='minimap2':
        indexOut = os.path.join(genomeFolder, '.'.join([genomeName, 'mmi']))
        build_command = ['minimap2', '-x', 'sr', '-d', indexOut, wholeGenome]
    else:
        build_command = ['bwa', 'index', wholeGenome]
    
    subprocess.check_call(' '.join(build_command), shell=True)
    

def map_core(fastq_1, fastq_2, ref, outdir, aligner='minimap2', outformat='SAM', 
             nthread=1):

    outformat = '.' + outformat.lower()
    # output file name
    if fastq_1.endswith('_1.fastq.gz'):
        outpath = os.path.join(outdir,
                        os.path.split(fastq_1)[1].replace('_1.fastq.gz',outformat))
    else:
        outpath = os.path.join(outdir,
                        os.path.split(fastq_1)[1].replace('_1.fastq',outformat))
    
    # ref: reference genome index
    if aligner=='minimap2':
        map_command = ['minimap2', '-ax', 'sr', '-t', str(nthread), ref, fastq_1, fastq_2]
    else:
        map_command = ['bwa', 'mem', '-SP', '-t', str(nthread), ref, fastq_1, fastq_2]
    
    if outformat=='.sam':
        bam_command = []
    else:
        bam_command = ['samtools', 'view', '-bS', '-']
    
    pipeline = []
    try:
        # Mapping
        pipeline.append(
                subprocess.Popen(map_command,
                    stdout=subprocess.PIPE if bam_command else open(outpath, 'wb'),
                    bufsize=-1))
        
        if bam_command:
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
    
    return outpath

def parse_bam(bam, outfile, chromsizes, assembly, min_mapq, max_molecule_size, max_inter_align_gap,
              walks_policy, include_readid, include_sam, drop_seq):
    
    basic_command = ['pairtools', 'parse', '-c', chromsizes, '--assembly', assembly,
                     '--min-mapq', str(min_mapq), '--max-molecule-size', str(max_molecule_size),
                     '--max-inter-align-gap', str(max_inter_align_gap), '--walks-policy', walks_policy]
    if not include_readid:
        basic_command.append('--drop-readid')
    
    if not include_sam:
        basic_command.append('--drop-sam')
    
    if drop_seq:
        basic_command.append('--drop-seq')
    basic_command.append(bam)
    
    pipeline = []
    try:
        pipeline.append(
            subprocess.Popen(basic_command,
                stdout=subprocess.PIPE,
                bufsize=-1)
        )

        sort_command = ['pairtools', 'sort', '-o', outfile, '--nproc', '8', '--memory', '2G']
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

        


