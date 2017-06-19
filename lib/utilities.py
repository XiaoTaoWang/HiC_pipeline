# Created on Tue Dec 23 21:15:19 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

import logging, os, subprocess, time, gc
import numpy as np
from mirnylib.h5dict import h5dict

log = logging.getLogger(__name__)

def sleep():
    
    for _ in range(3):
        time.sleep(0.1)
    gc.collect()
    for _ in range(3):
        time.sleep(0.1)
        
def juncSeqCountSRA(fastqPath, bash_reader, enzyme):
    """
    Scan each read to match so-called junction sequence.
    
    A paucity of ligation junctions in a Hi-C library suggests that the
    ligation failed.
    
    """
    import Bio.Restriction
    
    enzyme_site = eval('Bio.Restriction.%s.site' % enzyme)
    cutsite = eval('Bio.Restriction.%s.charac' % enzyme)[:2]
    
    Dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    reverse = ''.join([Dict[i] for i in enzyme_site])
    
    if cutsite[-1]:
        jplus = enzyme_site[:cutsite[-1]] + enzyme_site[cutsite[0]:]
        jminus = reverse[:cutsite[-1]] + reverse[cutsite[0]:]
        jminus = jminus[::-1]
    else:
        jplus = enzyme_site[:None] + enzyme_site[cutsite[0]:]
        jminus = reverse[:None] + reverse[cutsite[0]:]
        jminus = jminus[::-1]
    
    
    reading_command = bash_reader.split() + [fastqPath,]
    reading_process = subprocess.Popen(reading_command,
                                       stdout = subprocess.PIPE,
                                       bufsize = -1)
    
    Len = 0
    Count = 0
    lineNum = 0
    
    if jplus == jminus:
        for line in reading_process.stdout:
            if (lineNum % 8 == 1): 
                Len += 1
                check_1 = (jplus in line)
            if (lineNum % 8 == 5):
                check_2 = (jplus in line)
                check = check_1 or check_2
                if check:
                    Count += 1
            lineNum += 1
    else:
        for line in reading_process.stdout:
            if (lineNum % 8 == 1): 
                Len += 1
                check_1 = (jplus in line) or (jminus in line)
            if (lineNum % 8 == 5):
                check_2 = (jplus in line) or (jminus in line)
                check = check_1 or check_2
                if check:
                    Count += 1
            lineNum += 1
    
    sleep()
    
    return Len, Count

def juncSeqCountFASTQ(Fastq_1, Fastq_2, enzyme):
    """
    Scan each read to match so-called junction sequence.
    
    A paucity of ligation junctions in a Hi-C library suggests that the
    ligation failed.
    
    """
    import Bio.Restriction
    
    enzyme_site = eval('Bio.Restriction.%s.site' % enzyme)
    cutsite = eval('Bio.Restriction.%s.charac' % enzyme)[:2]
    
    Dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    reverse = ''.join([Dict[i] for i in enzyme_site])
    
    if cutsite[-1]:
        jplus = enzyme_site[:cutsite[-1]] + enzyme_site[cutsite[0]:]
        jminus = reverse[:cutsite[-1]] + reverse[cutsite[0]:]
        jminus = jminus[::-1]
    else:
        jplus = enzyme_site[:None] + enzyme_site[cutsite[0]:]
        jminus = reverse[:None] + reverse[cutsite[0]:]
        jminus = jminus[::-1]
    
    if Fastq_1.endswith('.fastq.gz'):
        pread_1 = subprocess.Popen(['gunzip', Fastq_1, '-c'],
                                    stdout = subprocess.PIPE, bufsize = -1)
    else:
        pread_1 = subprocess.Popen(['cat', Fastq_1],
                                    stdout = subprocess.PIPE, bufsize = -1)
    
    if Fastq_2.endswith('.fastq.gz'):
        pread_2 = subprocess.Popen(['gunzip', Fastq_2, '-c'],
                                    stdout = subprocess.PIPE, bufsize = -1)
    else:
        pread_2 = subprocess.Popen(['cat', Fastq_2],
                                    stdout = subprocess.PIPE, bufsize = -1)
                                  
    inStream_1 = pread_1.stdout
    inStream_2 = pread_2.stdout
    
    Len = 0
    Count = 0
    
    if jplus == jminus:
        while True:
            line = inStream_1.readline()
            
            try:
                assert line[0] == '@'
            except IndexError:
                break
            
            Block_1 = [line, inStream_1.readline(), inStream_1.readline(),
                       inStream_1.readline()]
            
            Block_2 = [inStream_2.readline(), inStream_2.readline(),
                       inStream_2.readline(), inStream_2.readline()]
            
            check_1 = (jplus in Block_1[1])
            check_2 = (jplus in Block_2[1])
            
            Len += 1
            if check_1 or check_2:
                Count += 1
    else:
        while True:
            line = inStream_1.readline()
            
            try:
                assert line[0] == '@'
            except IndexError:
                break
            
            Block_1 = [line, inStream_1.readline(), inStream_1.readline(),
                       inStream_1.readline()]
            
            Block_2 = [inStream_2.readline(), inStream_2.readline(),
                       inStream_2.readline(), inStream_2.readline()]
            
            check_1 = (jplus in Block_1[1]) or (jminus in Block_1[1])
            check_2 = (jplus in Block_2[1]) or (jminus in Block_2[1])
            
            Len += 1
            if check_1 or check_2:
                Count += 1
    
    sleep()
    
    return Len, Count

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

def uncompressSRA(sra, fastqF, bamF, hdf5F):
    
    if not commandExists('fastq-dump'):
        raise ValueError('Please install fastq-dump first!')
    
    pread = subprocess.Popen(['fastq-dump', sra, "-Z", "--split-files"],
                              stdout = subprocess.PIPE, bufsize = -1)
    
    inStream = pread.stdout
    
    outname = os.path.split(sra)[1].replace('.sra', '') + '_{0}.fastq.gz'
    outf1 = os.path.join(fastqF,outname).format(1)
    outf2 = os.path.join(fastqF,outname).format(2)
    outb1 = os.path.join(bamF,os.path.split(outf1)[1].replace('.fastq.gz','.bam'))
    outb2 = os.path.join(bamF,os.path.split(outf2)[1].replace('.fastq.gz','.bam'))
    hdf5 = os.path.join(hdf5F,os.path.split(sra)[1].replace('.sra','.hdf5'))
    outProc1 = gzipWriter(outf1)
    outProc2 = gzipWriter(outf2)
    outStream1 = outProc1.stdin
    outStream2 = outProc2.stdin
    
    count = 0
    while True:
        line = inStream.readline()
        
        try:
            assert line[0] == '@'
        except AssertionError:
            raise IOError('Invalid fastq file')
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
    
    yield outf1, outb1, outf2, outb2, hdf5

def linkRawFASTQ(fq1, fq2, bamF, hdf5F):
    
    if fq1.endswith('.fastq'):
        outb1 = os.path.join(bamF,os.path.split(fq1)[1].replace('.fastq','.bam'))
        hdf5 = os.path.join(hdf5F,os.path.split(fq1)[1].replace('_1.fastq','.hdf5'))
    else:
        outb1 = os.path.join(bamF,os.path.split(fq1)[1].replace('.fastq.gz','.bam'))
        hdf5 = os.path.join(hdf5F,os.path.split(fq1)[1].replace('_1.fastq.gz','.hdf5'))
    if fq2.endswith('.fastq'):
        outb2 = os.path.join(bamF,os.path.split(fq2)[1].replace('.fastq','.bam'))
    else:
        outb2 = os.path.join(bamF,os.path.split(fq2)[1].replace('.fastq.gz','.bam'))
    
    yield fq1, outb1, fq2, outb2, hdf5

def splitSRA(sra, chunkF, bamF, hdf5F, splitBy=4000000):
    
    if not commandExists('fastq-dump'):
        raise ValueError('Please install fastq-dump first!')

    outname = os.path.split(sra)[1].replace('.sra', '') + '_chunk{0}_{1}.fastq.gz'
    pread = subprocess.Popen(['fastq-dump', sra, "-Z", "--split-files"],
                              stdout = subprocess.PIPE, bufsize = -1)
    inStream = pread.stdout

    halted = False
    for counter in xrange(1000000):
        outf1 = os.path.join(chunkF,outname).format(counter,1)
        outf2 = os.path.join(chunkF,outname).format(counter,2)
        outb1 = os.path.join(bamF,os.path.split(outf1)[1].replace('.fastq.gz','.bam'))
        outb2 = os.path.join(bamF,os.path.split(outf2)[1].replace('.fastq.gz','.bam'))
        hdf5 = os.path.join(hdf5F,os.path.split(outf1)[1].replace('_1.fastq.gz','.hdf5'))
        outProc1 = gzipWriter(outf1)
        outProc2 = gzipWriter(outf2)
        outStream1 = outProc1.stdin
        outStream2 = outProc2.stdin

        for j in xrange(splitBy):

            line = inStream.readline()

            try:
                assert line[0] == '@'
            except AssertionError:
                raise IOError('Invalid fastq file')
            except IndexError:
                halted = True
                break

            fastq_entry = (line, inStream.readline(),
                           inStream.readline(), inStream.readline())

            outStream1.writelines(fastq_entry)
            outStream2.writelines((inStream.readline(), inStream.readline(),
                       inStream.readline(), inStream.readline()))

        outProc1.communicate()
        outProc2.communicate()
        
        yield outf1, outb1, outf2, outb2, hdf5
        
        if halted:
            break

def splitFASTQ(fq1, fq2, chunkF, bamF, hdf5F, splitBy=4000000):
    
    p1 = os.path.split(fq1)[1].split('.')[0].split('_')
    outname1 = '_'.join(p1[:-1]) + '_chunk{0}_{1}.fastq.gz'
    p2 = os.path.split(fq2)[1].split('.')[0].split('_')
    outname2 = '_'.join(p2[:-1]) + '_chunk{0}_{1}.fastq.gz'
    outname = '_'.join(p2[:-1]) + '_chunk{0}.hdf5'
    
    if fq1.endswith('.fastq.gz'):
        pread_1 = subprocess.Popen(['gunzip', fq1, '-c'],
                                    stdout = subprocess.PIPE, bufsize = -1)
    else:
        pread_1 = subprocess.Popen(['cat', fq1],
                                    stdout = subprocess.PIPE, bufsize = -1)
    if fq2.endswith('.fastq.gz'):
        pread_2 = subprocess.Popen(['gunzip', fq2, '-c'],
                                    stdout = subprocess.PIPE, bufsize = -1)
    else:
        pread_2 = subprocess.Popen(['cat', fq2],
                                    stdout = subprocess.PIPE, bufsize = -1)
                                  
    inStream_1 = pread_1.stdout
    inStream_2 = pread_2.stdout
    
    halted = False
    for counter in xrange(1000000):
        outf1 = os.path.join(chunkF,outname1).format(counter,p1[-1])
        outb1 = os.path.join(bamF,os.path.split(outf1)[1].replace('.fastq.gz','.bam'))
        outf2 = os.path.join(chunkF,outname2).format(counter,p2[-1])
        outb2 = os.path.join(bamF,os.path.split(outf2)[1].replace('.fastq.gz','.bam'))
        hdf5 = os.path.join(hdf5F,outname).format(counter)
        outProc1 = gzipWriter(outf1)
        outProc2 = gzipWriter(outf2)
        outStream1 = outProc1.stdin
        outStream2 = outProc2.stdin

        for j in xrange(splitBy):

            line = inStream_1.readline()

            try:
                assert line[0] == '@'
            except AssertionError:
                raise IOError('Invalid fastq file')
            except IndexError:
                halted = True
                break


            fastq_entry_1 = (line, inStream_1.readline(), inStream_1.readline(),
                             inStream_1.readline())
            fastq_entry_2 = (inStream_2.readline(), inStream_2.readline(),
                             inStream_2.readline(), inStream_2.readline())

            outStream1.writelines(fastq_entry_1)
            outStream2.writelines(fastq_entry_2)

        outProc1.communicate()
        outProc2.communicate()
        
        yield outf1, outb1, outf2, outb2, hdf5
        
        if halted:
            break
    
# Convert Matrix to Scipy Sparse Matrix
def toSparse(source, csr = False):
    """
    Convert intra-chromosomal contact matrices to sparse ones.
    
    Parameters
    ----------
    source : str
         Hdf5 file name.
    
    idx2label : dict
        A dictionary for conversion between zero-based indices and
        string chromosome labels.
    
    csr : bool
        Whether to use CSR (Compressed Row Storage) format or not.
    
    """
    import zipfile, tempfile
    from numpy.lib.format import write_array
    from scipy import sparse
    
    lib = h5dict(source, mode = 'r')
    
    ## Uniform numpy-structured-array format
    itype = np.dtype({'names':['bin1', 'bin2', 'IF'],
                          'formats':[np.int, np.int, np.float]})
    
    ## Create a Zip file in NPZ case
    if not csr:
        output = source.replace('.hm', '-sparse.npz')
    else:
        output = source.replace('.hm', '-csrsparse.npz')
    
    Zip = zipfile.ZipFile(output, mode = 'w', allowZip64 = True)
    fd, tmpfile = tempfile.mkstemp(suffix = '-numpy.npy')
    os.close(fd)
    
    log.log(21, 'Sparse Matrices will be saved to %s', output)
    log.log(21, 'Only intra-chromosomal matrices will be taken into account')
    log.log(21, 'Coverting ...')
    
    count = 0
    
    for i in lib:
        if (i != 'resolution') and (i != 'genomeInformation') and (len(set(i.split())) == 1):
            # Used for the dict-like key
            key = lib['genomeInformation']['idx2label'][int(i.split()[0])]
            
            log.log(21, 'Chromosome %s ...', key)
            # 2D-Matrix
            H = lib[i]
            
            if not csr:
                # Triangle Array
                Triu = np.triu(H)
                # Sparse Matrix in Memory
                x, y = np.nonzero(Triu)
                values = Triu[x, y]
                temp = np.zeros(values.size, dtype = itype)
                temp['bin1'] = x
                temp['bin2'] = y
                temp['IF'] = values
            else:
                temp = sparse.triu(H, format = 'csr')
            
            fname = key + '.npy'
            fid = open(tmpfile, 'wb')
            try:
                write_array(fid, np.asanyarray(temp))
                fid.close()
                fid = None
                Zip.write(tmpfile, arcname = fname)
            finally:
                if fid:
                    fid.close()
                    
            log.log(21, 'Done!')
            
            count += 1
    
    if count == 0:
        log.warning('Empty source file!')
            
    # Other information
    for i in ['resolution', 'genomeInformation']:
        fname = '.'.join([i, 'npy'])
        fid = open(tmpfile, 'wb')
        try:
            write_array(fid, np.asanyarray(lib[i]))
            fid.close()
            fid = None
            Zip.write(tmpfile, arcname = fname)
        finally:
            if fid:
                fid.close()
    
    os.remove(tmpfile)
    
    Zip.close()