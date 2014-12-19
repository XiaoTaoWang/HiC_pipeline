#!/usr/bin/env python

# Created on Tue Dec 16 10:22:41 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

## Required Modules
import os, sys, argparse, logging, logging.handlers
import numpy as np

def getargs():
    ## Construct an ArgumentParser object for command-line arguments
    parser = argparse.ArgumentParser(description = '''This software is based on hiclib
                                    (https://bitbucket.org/mirnylab/hiclib), a comprehensive
                                    Python package for Hi-C data analysis.''',
                                    formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    # Version
    parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s 0.1.0',
                        help = 'Print version number and exit')
    
    ## Common Arguments
    parser.add_argument('-p', '--dataFolder', default = '.',
                        help = '''Root directory of original data. We recommend placing sequencing
                        and genome data here.''')
    parser.add_argument('-g', '--genomeName',
                        help = '''Genome folder name. This folder must be placed under dataFolder.
                        Genome sequences should be stored chromosome by chromosome in FASTA format.
                        If gap file is not contained, we will generate a dummy one.''')
    parser.add_argument('-C', '--chroms', nargs = '*', default = ['#', 'X'],
                       help = '''Which chromosomes will be involved. Specially, "#" stands for
                       chromosomes with numerical labels. "--chroms" with zero argument will
                       generate an empty list, in which case all chromosome data will be loaded.''')
    parser.add_argument('-T', '--template', default = 'chr%s.fa',
                        help = '''Template of FASTA file names''')
    parser.add_argument('-G', '--gapFile', default = 'gap.txt', help = '''Gap file name.''')
    
    ## Sub-commands
    subparser = parser.add_subparsers(title = 'sub-commands',
                                      description = '''Read pair mapping, filtering, binning
                                      and iterative correction are contained. You can perform
                                      each stage of the analysis separately, or streamline the
                                      pipeline using "pileup" subcommand.''',
                                      dest = 'subcommand')
    ## Iterative Mapping
    iterM = subparser.add_parser('mapping',
                                 help = '''Map raw pair-end sequencing data to a supplied
                                 genome. Both SRA and FASTQ format are admissible.''',
                                 description = '''An iterative mapping schema is used. The
                                 minimum length is always 25, then the step will be calculated
                                 automatically based on the sequence length. The bowtie2 mapping
                                 software and a fastq-dump tool from SRA toolkit are required.
                                 At least, you should specify --fastqDir, --genomeName,
                                 --bowtiePath, --dataFolder and --metadata yourself.''',
                                 epilog = '''After this command, a BAM folder containing BAM
                                 files for each side of Hi-C molecules and a HDF5 folder containing
                                 hdf5 (dict-like structure format) files for library of matched
                                 Hi-C reads are created under current working directory.''')
    iterM.add_argument('-f', '--fastqDir', help = 'Sequencing data folder. Relative path to dataFolder')
    iterM.add_argument('-F', '--Format', default = 'SRA', choices = ['SRA', 'FASTQ'],
                       help = 'Format of the sequencing data.')
    iterM.add_argument('-m', '-metadata', default = 'datasets.tsv',
                       help = '''Metadata file describing each SRA file. You should place
                       it under current working directory. Four columns are required: prefix
                       of SRA file name, cell line name, biological replicate label, and
                       restriction enzyme name. An example file is distributed along with
                       this software, please check it.''')
    iterM.add_argument('-b', '--bowtiePath', help = 'Path to bowtie2 executable program file.')
    iterM.add_argument('-t', '--threads', type = int, default = 4, help = 'Number bowtie2 threads.')
    iterM.add_argument('-i', '--bowtieIndex',
                       help = '''Path to the bowtie2 genome index. Since the index consists of
                       several files with the different suffices (e.g., hg19.1.bt2, hg19.2.bt.2),
                       provide only the common part. For example, if your genome data hg19.fa
                       and corresponding index files are stored in ~/data/hg19, you need to
                       specify --bowtieIndex as this "--bowtieIndex ~/data/hg19/hg19". When not
                       specified, we will generate one under the genome folder.''')
    iterM.add_argument('--cache', default = '/tmp',
                       help = ''''Set the cache folder. Absolute path is needed.''')
    iterM.set_defaults(func = mapping)
    
    ## Merge files from the same experiment
    multiF = subparser.add_parser('merge',
                                  help = '''Merge files corresponding to the same experiment
                                  together.''',
                                  description = '''This command is useful when you want to merge
                                  several HDF5 files belonging to the same experiment. (Metadata
                                  file is used again)''',
                                  epilog = '''A folder with one or more merged hdf5 files are
                                  generated under current working directory after this command
                                  is called.''')
    multiF.add_argument('-h', '--HDF5',
                        help = '''Path to the folder with hdf5 files which are generated by
                        mapping command.''')
    multiF.add_argument('-m', '-metadata', default = 'datasets.tsv',
                        help = '''Metadata file describing each SRA file. You should place
                        it under current working directory. Four columns are required: prefix
                        of SRA file name, cell line name, biological replicate label, and
                        restriction enzyme name. An example file is distributed along with
                        this software, please check it.''')
    multiF.add_argument('-l', '--level', type = int, default = 2, choices = [1, 2],
                        help = '''Set merging level. 1: hdf5 files from the same biological
                        replicate will be merged, 2: hdf5 files from the same cell line will be
                        merged.''')
    multiF.set_defaults(func = merge)
    
    ## Fragment-level filtering
    removeNoise = subparser.add_parser('filtering',
                                       help = '''Filtering at the level of aligned read pairs
                                       and restriction fragments.''',
                                       description = '''PCR duplications, self-ligation products,
                                       unligated "dangling end" products, random breaks, too
                                       large and too small fragments, and fragments with high
                                       cis-to-trans ratio are all taken into account.''',
                                       epilog = '''A folder with corresponding filtered hdf5
                                       files are generated under current working directory after
                                       calling this command.''')
    removeNoise.add_argument('-u', '--mergedDir',
                             help = '''Path to the merged HDF5 files generated by merge command.
                             If the path points to one file, filtering will only be performed
                             on that file. If the path points to a folder, filtering will be
                             performed on all files contained in that folder.''')
    removeNoise.add_argument('--duplicates', action = 'store_false',
                             help = '''Remove read pairs resulting from PCR amplification.''')
    removeNoise.add_argument('--startNearRsite', action = 'store_true',
                             help = '''Remove reads that start within 5 bp near a restriction site.''')
    removeNoise.add_argument('--RandomBreaks', action = 'store_false',
                             help = '''Remove "random breaks" in which corresponding fragments
                             did not arise from normal restriction digestion.''')
    removeNoise.add_argument('--extremeFragments', action = 'store_false',
                             help = '''Remove too large and too small fragments.''')
    removeNoise.add_argument('--cistotrans', action = 'store_false',
                             help = '''Remove top 0.5% of fragments with the greatest number of
                             reads.''')
    removeNoise.set_defaults(func = filtering)
    
    ## Binning
    binReads = subparser.add_parser('binning',
                                    help = '''Bin filtered reads at certain resolution.''',
                                    description = '''For varying resolutions, three modes are
                                    provided, just choose a proper one.''',
                                    epilog = '''After calling this command, a folder with
                                    created HeatMaps (in HDF5 format) is created under current
                                    working directory.''')
    binReads.add_argument('-f', '--filteredDir',
                          help = '''Path to the filtered HDF5 files generated by filtering
                          command. If the path points to one file, we will only create a HeatMap
                          for that file. If the path points to a folder, we will construct a
                          HeatMap for each file in that folder.''')
    binReads.add_argument('-M', '--mode', default = 'wholeGenome',
                          choices = ['wholeGenome', 'byChromosome', 'withOverlaps'],
                          help = '''Memory usage: withOverlaps > byChromosome > wholeGenome.
                          Resolution capacity (take human genome for example):
                          withOverlaps (10kb) > byChromosome (40kb) > wholeGenome (200kb).''')
    binReads.add_argument('-R', '--resolution', type = int, default = 200000,
                          help = 'Resolution of a heatmap. Unit: bp')
    binReads.set_defaults(func = binning)
    
    ## Iterative Correction
    iterC = subparser.add_parser('correcting',
                                 help = '''Perform iterative corrections on the original HeatMap.''',
                                 description = '''Two modes are provided for different resolutions.
                                 The program will choose a better one for you according to the data
                                 format.''',
                                 epilog = '''After calling this command, a folder with corrected
                                 HeatMaps (in HDF5 format) is created under current working
                                 directory.''')
    iterC.add_argument('-H', '--HeatMap',
                       help = '''Path to the HeatMap files generated by binning command. If the
                       path points to one file, we only correct for that HeatMap. If the path
                       points to a folder, we will perform iterative corrections for all HeaMaps
                       in that folder.''')
    iterC.set_defaults(func = correcting)
                    
    ## Pile Up
    streamline = subparser.add_parser('pileup',
                                      parents = [iterM],
                                      help = '''Perform the entire analysis from sequencing
                                      data to corrected HeatMaps.''',
                                      description = '''A more convenient but less flexible
                                      command for Hi-C data processing.''')
    streamline.set_defaults(func = pileup)
    
    ## Parse the command-line arguments
    commands = sys.argv[1:]
    if not commands:
        commands.append('-h')
    args = parser.parse_args(commands)
    
    return args, commands

def run():
    ## Necessary Modules
    from mirnylib import genome, h5dict
    # Parse Arguments
    args, commands = getargs()
    # Improve the performance if you don't want to run it
    if commands[0] not in ['-h', '-v', '--help', '--version']:
        ## Validity of arguments
        dataLocation = os.path.abspath(os.path.expanduser(args.dataFolder))
        if not os.path.exists(dataLocation):
            logging.error('There is no folder named %s on your system!',
                          dataLocation)
            sys.exit(1)
        genomeFolder = os.path.join(dataLocation, args.genomeName)
        if not os.path.exists(genomeFolder):
            logging.error('%s can not be found at %s', args.genomeName,
                          dataLocation)
            sys.exit(1)
    
        ## Generate a dummy gap file under genome folder if there's no one yet
        gapFile = os.path.join(genomeFolder, args.gapFile)
        if not os.path.exists(gapFile):
            logging.info('No gap file can be found at %s, generating a dummy one ...',
                         genomeFolder)
            tempfile = open(gapFile, 'w')
            tempfile.write('0\tNA1000\t0\t0\t0\tN\t0\tcentromere\tno\n')
            tempfile.flush()
            tempfile.close()
            logging.info('Done!')
        
        # Python Genome Object
        genome_db = genome.Genome(genomeFolder, readChrms = args.chroms)
    

def mapping(args, commands):
    ## Import necessary modules
    import atexit
    import hiclib.mapping as iterM
    
    # A Local Function
    def cleanFile(filename):
        if os.path.exists(filename):
            os.remove(filename)
    
    # Construct bowtie2 genome index
    def buildIndex(genomeFolder):
        """
        Build bowtie2 index files under the provided genome folder.
        
        """
        import glob
        fastaNames = [os.path.join(genomeFolder, i)
                      for i in glob.glob(os.path.join(
                      genomeFolder, args.template % ('*',)))]
        wholeGenome = os.path.join(genomeFolder,
                                   '.'.join([args.genomeName, 'fa']))
        if not os.path.exists(wholeGenome):
            os.system('cat ' + ' '.join(fastaNames) + ' > ' + wholeGenome)
        bowtieIndex = os.path.join(genomeFolder, args.genomeName)
        buildCmd = ['bowtie2-build', '--quiet', wholeGenome, bowtieIndex]
        os.system(' '.join(buildCmd))
        
        return bowtieIndex
    
    def calculateStep(length, minlen, approxStep=10, maxSteps=4):
        """
        Returns minimum length and step based on the length of sequence and
        proposed minimum length.
        """
        actualDif = length - minlen
        if actualDif < approxStep * 0.6:
            return length, 100

        numIter = np.array(np.around(actualDif / float(approxStep)), dtype=int)
        if numIter == 0:
            numIter = 1
        if numIter > maxSteps:
            numIter = maxSteps
        actualStep = actualDif / numIter

        minlen = length - actualStep * numIter

        return minlen, actualStep
    
    ## Validity of arguments
    bowtiePath = os.path.abspath(os.path.expanduser(args.bowtiePath))
    if not os.path.exists(bowtiePath):
        logging.error('Bowtie2 can not be found at %s', bowtiePath)
        sys.exit(1)
    fastqDir = os.path.join(dataLocation, args.fastqDir)
    if not os.path.exists(fastqDir):
        logging.error('%s should be placed under %s', args.fastqDir, dataLocation)
        sys.exit(1)
    mFile = args.metadata
    if not os.path.exists(mFile):
        logging.error('Metadata file %s can not be found at current working directory!',
                      mFile)
        sys.exit(1)
    cache = os.path.abspath(os.path.expanduser(args.cache))
    if not os.path.exists(cache):
        logging.warning('%s does not exist on your system, trying to create one',
                        cache)
        os.makedirs(cache)
    
    ## Construct bowtie2 genome index if there's no one yet
    if '--bowtieIndex' in commands:
        bowtieIndex = os.path.abspath(os.path.expanduser(args.bowtieIndex))
    else:
        logging.info('Index path is not provided, generating under the'
                     ' genome folder ...')
        bowtieIndex = buildIndex(genomeFolder)
        logging.info('Done!')
    
    ## Sequencing Data Format
    Format = args.Format.lower()
    sraNames = [os.path.join(fastqDir, i) for i in glob.glob(os.path.join(
                fastqDir, '%s.%s' % ('*', Format)))]
    ## Sequencing Length
    lengths = os.path.join(fastqDir, 'lengths')
    if not os.path.exists(lengths):
        os.mkdir(lengths)
    if Format == 'sra':
        Set = set([os.path.basename(i)[:-4] for i in sraNames])
        for i in sraNames:
            calLength = ['fastq-dump', '-Z', i, '|', 'head', '-n', '2',
                         '|', 'tail', '-n', '1', '|', 'wc', '-c', '>',
                         os.path.join(lengths, os.path.basename(i)[:-4])]
            os.system(' '.join(calLength))
    else:
        Set = set([os.path.basename(i)[:-8] for i in sraNames])
        for i in Set:
            leftSide = os.path.join(fastqDir, i + '_1.fastq')
            calLength = ['head', '-n', '2', leftSide, '|', 'tail', '-n', '1',
                         'wc', '-c', '>', os.path.join(lengths, i)]
            os.system(' '.join(calLength))
    
    ## Output Folders
    bamFolder = 'bams-%s' % args.genomeName
    hdf5F = 'hdf5-%s' % args.genomeName
    args.HDF5 = hdf5F # To communicate with next processing step (merge)
    if not os.path.exists(bamFolder):
        os.mkdir(bamFolder)
    if not os.path.exists(hdf5F):
        os.mkdir(hdf5F)
    
    # Read Metadata
    metadata = [l.rstrip().split() for l in open(mFile)]
    database = dict([(i[0], i[-1]) for i in metadata])
    for i in 3 * sorted(list(Set)):
        # Parameters used in iterative mapping
        lengthFile = os.path.join(lengths, i)
        if Format == 'sra':
            length = (int(open(lengthFile).readlines()[0]) - 1) / 2
        else:
            length = int(open(lengthFile).readlines()[0])
        minlen, step = calculateStep(length, 25)
        
        finalFile = os.path.join(hdf5F, '%s.hdf5' % i)
        lockFile = os.path.join(hdf5F, '%s.lock' % i)
        
        if os.path.exists(finalFile) and not os.path.exists(lockFile):
            logging.info('% is there, skipping', finalFile)
            continue
        
        if os.path.exists(lockFile):
            logging.info('Someone is working on %s', finalFile)
            continue
        
        lock = open(lockFile, 'w')
        lock.close()
        
        atexit.register(cleanFile, lockFile)
        cleanup = ['rm', '-rf', os.path.join(bamFolder, '%s*' % i)]
        os.system(' '.join(cleanup))
        
        ## Iterative Mapping
        # Common Parameters
        Parameters = {'bowtie_path': bowtiePath, 'bowtie_index_path': bowtieIndex,
                      'min_seq_len': minlen, 'len_step': step, 'nthreads': args.threads,
                      'temp_dir': cache, 'bowtie_flags': '--very-sensitive'}
        if Format == 'sra':
            sourceFile = os.path.join(fastqDir, i + '.sra')
            # The First Side
            iterM.iterative_mapping(fastq_path = sourceFile,
                                    out_sam_path = '%s/%s_1.bam' % (bamFolder, i),
                                    seq_start = 0,
                                    seq_end = length,
                                    bash_reader = 'fastq-dump -Z',
                                    **Parameters)
            # The Second Side
            iterM.iterative_mapping(fastq_path = sourceFile,
                                    out_sam_path = '%s/%s_2.bam' % (bamFolder, i),
                                    seq_start = length,
                                    seq_end = 2 * length,
                                    bash_reader = 'fastq-dump -Z',
                                    **Parameters)
        else:
            iterM.iterative_mapping(fastq_path = os.path.join(fastqDir, i + '_1.fastq'),
                                    out_sam_path = '%s/%s_1.bam' % (bamFolder, i),
                                    **Parameters)
            iterM.iterative_mapping(fastq_path = os.path.join(fastqDir, i + '_2.fastq'),
                                    out_sam_path = '%s/%s_2.bam' % (bamFolder, i),
                                    **Parameters)
        
        ## Parse the mapped sequences into a Python data structure
        ## Assign the ultra-sonic fragments to restriction fragments
        lib = h5dict.h5dict(finalFile)
        iterM.parse_sam(sam_basename1 = '%s/%s_1.bam' % (bamFolder, i),
                        sam_basename2 = '%s/%s_2.bam' % (bamFolder, i),
                        out_dict = lib,
                        genome_db = genome_db,
                        enzyme_name = database[i],
                        save_seqs = False)
        
        os.remove(lockFile)

def merge(args):
    # Required Modules
    from hiclib.fragmentHiC import HiCdataset
    
    ## Validity of arguments
    Sources = os.path.abspath(os.path.expanduser(args.HDF5))
    if not os.path.exists(Sources):
        logging.error('There is no folder named %s on your system!', Sources)
        sys.exit(1)
    mFile = args.metadata
    if not os.path.exists(mFile):
        logging.error('%s can not be found under current working directory!', mFile)
        sys.exit(1)
    # Output Folder
    mergedFolder = 'merged-%s' % args.genomeName
    if not os.path.exists(mergedFolder):
        os.mkdir(mergedFolder)
    args.mergedDir = mergedFolder # To communicate with next processing step (filtering)
    ## Now, start merging
    metadata = [l.rstrip().split() for l in open(mFile)]
    ## Hierarchical merging structures
    bioReps = set((i[1], i[3], i[2]) for i in metadata)
    cellLines = set((i[1], i[3]) for i in metadata)
    ## The First level, biological replicates
    queueL1 = []
    for rep in bioReps:
        filenames = [os.path.join(Sources, '%s.hdf5' % i[0]) for i in metadata
                    if ((i[1], i[3], i[2]) == rep)]
        outfile = os.path.join(mergedFolder, '%s-%s-%s-merged.hdf5' % rep)
        enzyme = rep[1]
        queueL1.append((filenames, outfile, enzyme))
    for member in queueL1:
        ## Parsing individual files, no filtering processes are applied.
        lanePools = []
        for source in member[0]:
            parseName = os.path.join(mergedFolder, '%s-parsed.hdf5' % os.path.basename(source).replace('.hdf5', ''))
            parseF = HiCdataset(filename = parseName,
                                genome = genome_db,
                                mode = 'w')
            parseF.parseInputData(source, enzymeToFillRsites = member[-1], noFiltering = True)
            lanePools.append(parseName)
        ## Merge files altogether
        fragments = HiCdataset(filename = member[1], genome = genome_db, mode = 'w')
        fragments.merge(lanePools)
        # Clean up parsed individual files
        for delFile in lanePools:
            os.remove(delFile)
    
    if args.level == 2:
        ## The Second level, cell lines, optional
        queueL2 = []
        for cell in cellLines:
            filenames = [os.path.join(mergedFolder, '%s-%s-%s-merged.hdf5' % i) for i in bioReps
                         if ((i[1], i[3]) == cell)]
            outfile = os.path.join(mergedFolder, '%s-%s-allReps-merged.hdf5' % cell)
            queueL2.append((filenames, outfile))
        
        for member in queueL2:
            fragments = HiCdataset(filename = member[1], genome = genome_db, mode = 'w')
            fragments.merge(member[0])

def filtering(args):
    ## Necessary Modules
    from hiclib.fragmentHiC import HiCdataset
    import glob
    
    def core(filename):
        # Read data from unfiltered library
        lib = HiCdataset(filename, genome = genome_db, mode = 'r')
        ## Retrive dangling ends and estimate library size of sequencing
        DSmask = (lib.chrms1 >= 0) * (lib.chrms2 >= 0)
        sameFragMask = (lib.fragids1 == lib.fragids2) * DSmask
        cutDifs = lib.cuts2[sameFragMask] > lib.cuts1[sameFragMask]
        s1 = lib.strands1[sameFragMask]
        s2 = lib.strands2[sameFragMask]
        SSDE = (s1 != s2)
        SS = SSDE * (cutDifs == s2)
        Dangling = SSDE & (~SS)
        dist = lib.evaluate("a = - cuts1 * (2 * strands1 -1) - "
                            "cuts2 * (2 * strands2 - 1)",
                            ["cuts1", "cuts2", "strands1", "strands2"])
        Dangling_L = dist[sameFragMask][Dangling]
        library_L = int(np.ceil((np.percentile(Dangling_L, 95))))
        ## Create a new library for filtering
        filteredF = os.path.join(filteredFolder, os.path.basename(filename).replace('merged', 'filtered'))
        fragments = HiCdataset(filteredF, genome = genome_db,
                               maximumMoleculeLength = library_L,
                               mode = 'w')
        fragments.parseInputData(filename, removeSS = True)
        ## Additional Filtering
        if args.duplicates:
            fragments.filterDuplicates()
        if args.startNearRsite:
            fragments.filterRsiteStart(offset=5)
        if args.RandomBreaks:
            t_N = fragments.N
            fragments.maskFilter((fragments.dists1 + fragments.dists2) <= library_L)
            fragments.metadata["330_removeRandomBreaks"] = t_N - fragments.N
        if args.extremeFragments:
            fragments.filterLarge()
        if args.cistotrans:
            fragments.filterExtreme(cutH=0.005, cutL=0)
        
    ## Validity of arguments
    Sources = os.path.abspath(os.path.expanduser(args.mergedDir))
    if not os.path.exists(Sources):
        logging.error('%s does not exists on your system!', Sources)
        sys.exit(1)
    
    # Output Folder
    filteredFolder = 'filtered-%s' % args.genomeName
    if not os.path.exists(filteredFolder):
        os.mkdir(filteredFolder)
    args.filteredDir = filteredFolder # To communicate with next processing step (binning)
    
    ## Two cases: a directory or a single file
    if os.path.isdir(Sources):
        queue = [os.path.join(Sources, i) for i in glob.glob(os.path.join(Sources, '*-merged.hdf5'))]
        if len(queue) == 0:
            logging.error('No proper files can be found at %s!', Sources)
            sys.exit(1)
        else:
            for f in queue:
                core(f)
    else:
        core(Sources)

def binning(args):
    ## Necessary Modules
    from hiclib.fragmentHiC import HiCdataset
    import glob
    
    ## Validity of arguments
    Sources = os.path.abspath(os.path.expanduser(args.filteredDir))
    if not os.path.exists(Sources):
        logging.error('%s does not exists on your system!', Sources)
        sys.exit(1)
    
    # Output Dir
    hFolder = 'Heatmaps-%s' % args.genomeName
    if not os.path.exists(hFolder):
        os.mkdir(hFolder)
    # To communicate with next processing step (correcting)
    args.HeatMap = hFolder
    
    ## Generate HeatMaps
    if os.path.isdir(Sources):
        queue = [os.path.join(Sources, i) for i in glob.glob(os.path.join(Sources, '*-filtered.hdf5'))]
        if len(queue) == 0:
            logging.error('No proper files can be found at %s!', Sources)
            sys.exit(1)
    else:
        queue = [Sources]
    
    # Appropriate Units
    unit, denominator = ('K', 1000) if (args.resolution / 1000 < 1000) else ('M', 1000000)
    nLabel = str(args.resolution / denominator) + unit
    for f in queue:
        hFile = os.path.join(hFolder, os.path.basename(f).replace('.hdf5', '-%s.hm' % nLabel))
        fragments = HiCdataset(f, genome = genome_db, mode = 'r')
        ## Different Modes
        if args.mode == 'wholeGenome':
            fragments.saveHeatmap(hFile, resolution = args.resolution)
        if args.mode == 'byChromosome':
            fragments.saveByChromosomeHeatmap(hFile, resolution = args.resolution)
        if args.mode == 'withOverlaps':
            fragments.saveHiResHeatmapWithOverlaps(hFile, resolution = args.resolution)

def correcting(args):
    ## Two modes
    def Lcore(filename):
        # Necessary Modules
        from hiclib import binnedData
        # Output file
        cFile = os.path.join(cFolder, os.path.basename(filename).replace('.hm', '_c.hm'))
        # Create a binnedData object, load the data.
        BD = binnedData.binnedData(resolution, genome_db)
        name = '-'.join(os.path.basename(filename).split('-')[:3])
        BD.simpleLoad(filename, name)
        ## Perform ICE
        # Remove the contacts between loci located within the same bin.
        BD.removeDiagonal()
        # Remove bins with less than half of a bin sequenced.
        BD.removeBySequencedCount(0.5)
        # Remove 0.5% bins with the lowest number of records
        BD.removePoorRegions(cutoff = 0.5, coverage = True)
        BD.removePoorRegions(cutoff = 0.5, coverage = False)
        # Truncate top 0.05% of inter-chromosomal counts (possibly, PCR blowouts).
        BD.truncTrans(high = 0.0005)
        # Perform iterative correction.
        BD.iterativeCorrectWithoutSS()
        # Save the iteratively corrected heatmap.
        BD.export(name, cFile)

    def Hcore(filename):
        # Necessary Modules
        from hiclib.highResBinnedData import HiResHiC
        # Output file
        cFile = os.path.join(cFolder, os.path.basename(filename).replace('.hm', '_c.hm'))
        # Intermediate file
        iFile = os.path.join(cFolder, os.path.basename(filename).replace('.hm', '_i.hm'))
        ## Perform ICE
        BD = HiResHiC(genome_db, resolution, iFile, mode='w')
        BD.loadData(dictLike = filename)
        BD.removeDiagonal()
        ## Remove x% poor regions ...
        temp = BD.getMarginals(normalizeForIC = False)
        allMarg = np.concatenate(temp)
        percent = (1 - (allMarg >= 100).sum() / (allMarg > 0).sum()) * 100
        BD.removePoorRegions(percent)
        BD.iterativeCorrection(1e-2)
        BD.export(cFile)
        # Remove intermediate files
        os.remove(iFile)
        
    ## Validity of arguments
    Sources = os.path.abspath(os.path.expanduser(args.HeatMap))
    if not os.path.exists(Sources):
        logging.error('%s does not exists on your system!', Sources)
        sys.exit(1)
    
    ## Output Dir
    cFolder = 'Corrected-%s' % args.genomeName
    if not os.path.exists(cFolder):
        os.mkdir(cFolder)
    
    ## Corrections start
    if os.path.isdir(Sources):
        queue = [os.path.join(Sources, i) for i in glob.glob(os.path.join(Sources, '*.hm'))]
        if len(queue) == 0:
            logging.error('No proper files can be found at %s!', Sources)
            sys.exit(1)
    else:
        queue = [Sources]
    
    for f in queue:
        # Raw Data
        raw = h5dict.h5dict(f, mode = 'r')
        Keys = raw.keys()
        if 'heatmap' in Keys: # Low resolution case
            resolution = int(raw['resolution'])
            Lcore(f)
        else: # High resolution case
            rlabel = os.path.basename(f).split('-')[-1].replace('.hm', '')
            ends = rlabel[-1]
            resolution = (int(rlabel[:-1]) * 1000) if (ends == 'K') else (int(rlabel[:-1]) * 1000000)
            Hcore(f)
            
def pileup(args):
    mapping(args)
    args.level = 2
    merge(args)
    filtering(args)
    args.mode = 'withOverlaps'
    args.resolution = 10000
    binning(args)
    correcting(args)
    

if __name__ == '__main__':
    run()