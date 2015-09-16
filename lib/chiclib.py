# Created on Mon May 25 16:57:52 2015

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

import logging, os, time, re
import matplotlib
matplotlib.use('Agg')
import numpy as np
from mirnylib.genome import Genome
from hiclib.fragmentHiC import HiCdataset
from hiclib.highResBinnedData import HiResHiC
from mirnylib.numutils import uniqueIndex, fillDiagonal, fasterBooleanIndexing
from mirnylib.h5dict import h5dict

## Plot Settings
import matplotlib.pyplot as plt
# Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['axes.labelsize'] = 13
matplotlib.rcParams['xtick.labelsize'] = 13
matplotlib.rcParams['ytick.labelsize'] = 13
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['ytick.minor.size'] = 5
matplotlib.rcParams['xtick.major.pad'] = 6
matplotlib.rcParams['ytick.major.pad'] = 6
matplotlib.rcParams['xtick.minor.pad'] = 6
matplotlib.rcParams['ytick.minor.pad'] = 6

colorPool = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C']

log = logging.getLogger(__name__)

class myGenome(Genome):
    
    def _extractChrmLabel(self, fastaName):
        # First assume a whole filename as input (e.g. 'chr01.fa')
        _, fastaName = os.path.split(fastaName)
        regexp = self.chrmFileTemplate % ('(.*)')
        search_results = re.search(regexp, fastaName)
        # If not, assume that only the name is supplied as input (e.g. 'chr01')
        if search_results is None:
            regexp = self.chrmFileTemplate.split('.')[0] % ('(.*)')
            search_results = re.search(regexp, fastaName)

        if search_results is None:
            raise Exception(
                'The filename {} does not match the template {}.'.format(
                fastaName, self.chrmFileTemplate))

        chrm_label = search_results.group(1)

        return chrm_label

# A Customized HiCdataset Class
class cHiCdataset(HiCdataset):
    
    def __init__(self, filename, genome, maximumMoleculeLength = 500,
                 inMemory = False, mode = "a"):
        
        self.vectors = {
            # chromosomes for each read.
            "chrms1": "int8", "chrms2": "int8",

            "mids1": "int32", "mids2": "int32",
            # midpoint of a fragment, determined as "(start+end)/2"

            "fraglens1": "int32", "fraglens2": "int32",
            # fragment lengthes

            "distances": "int32",
            # distance between fragments. If -1, different chromosomes.
            # If -2, different arms.

            "fragids1": "int64", "fragids2": "int64",
            # IDs of fragments. fragIDmult * chromosome + location
            # distance to rsite
            "dists1": "int32", "dists2": "int32",
            # precise location of cut-site
            "cuts1": "int32", "cuts2": "int32",
            "strands1": "bool", "strands2": "bool",
            }
        self.metadata = {}

        #-------Initialization of the genome and parameters-----
        self.mode = mode
        self.genome = genome

        self.chromosomeCount = self.genome.chrmCount
        self.fragIDmult = self.genome.fragIDmult  # used for building heatmaps

        self.maximumMoleculeLength = maximumMoleculeLength

        self.filename = os.path.abspath(os.path.expanduser(filename))  # File to save the data
        self.chunksize = 5000000
        # Chunk size for h5dict operation, external sorting, etc.

        self.inMemory = inMemory

        self.h5dict = h5dict(self.filename, mode = mode, in_memory = inMemory)
        
        if 'chrms1' in self.h5dict.keys():
            chrms1 = self.chrms1
            self.DSnum = self.N = len(chrms1)
    
    def parseInputData(self, dictLike, **kwargs):
        
        import numexpr
        
        if not os.path.exists(dictLike):
            raise IOError('File not found: %s' % dictLike)
        
        dictLike = h5dict(dictLike, 'r')
        
        self.chrms1 = dictLike['chrms1']
        self.chrms2 = dictLike['chrms2']
        self.cuts1 = dictLike['cuts1']
        self.cuts2 = dictLike['cuts2']
        self.strands1 = dictLike['strands1']
        self.strands2 = dictLike['strands2']
        self.dists1 = np.abs(dictLike['rsites1'] - self.cuts1)
        self.dists2 = np.abs(dictLike['rsites2'] - self.cuts2)
        self.mids1 = (dictLike['uprsites1'] + dictLike['downrsites1']) / 2
        self.mids2 = (dictLike['uprsites2'] + dictLike['downrsites2']) / 2
        self.fraglens1 = np.abs(
            (dictLike['uprsites1'] - dictLike['downrsites1']))
        self.fraglens2 = np.abs(
            (dictLike['uprsites2'] - dictLike['downrsites2']))
        self.fragids1 = self.mids1 + np.array(self.chrms1,
                                              dtype='int64') * self.fragIDmult
        self.fragids2 = self.mids2 + np.array(self.chrms2,
                                              dtype='int64') * self.fragIDmult
        
        distances = np.abs(self.mids1 - self.mids2)
        distances[self.chrms1 != self.chrms2] = -1
        self.distances = distances  # Distances between restriction fragments
        del distances
        
        self.N = len(self.chrms1)

        try:
            dictLike['misc']['genome']['idx2label']
            self.updateGenome(self.genome,
                              oldGenome = dictLike["misc"]["genome"]["idx2label"])
        except KeyError:
            assumedGenome = myGenome(self.genome.genomePath)
            self.updateGenome(self.genome, oldGenome = assumedGenome)

        # Discard dangling ends and self-circles
        DSmask = (self.chrms1 >= 0) * (self.chrms2 >= 0)
        self.metadata['100_NormalPairs'] = DSmask.sum()

        sameFragMask = self.evaluate("a = (fragids1 == fragids2)",
                     ["fragids1", "fragids2"]) * DSmask

        cutDifs = self.cuts2[sameFragMask] > self.cuts1[sameFragMask]
        s1 = self.strands1[sameFragMask]
        s2 = self.strands2[sameFragMask]
        SSDE = (s1 != s2)
        SS = SSDE * (cutDifs == s2)
        SS_N = SS.sum()
        SSDE_N = SSDE.sum()
        sameFrag_N = sameFragMask.sum()
        self.metadata['120_SameFragmentReads'] = sameFrag_N
        self.metadata['122_SelfLigationReads'] = SS_N
        self.metadata['124_DanglingReads'] = SSDE_N - SS_N
        self.metadata['126_UnknownMechanism'] = sameFrag_N - SSDE_N
        
        mask = DSmask * (-sameFragMask)

        del DSmask, sameFragMask
        
        noSameFrag = mask.sum()
        
        # distance between sites facing each other
        dist = self.evaluate("a = numexpr.evaluate('- cuts1 * (2 * strands1 -1) - "
                             "cuts2 * (2 * strands2 - 1)')",
                             ["cuts1", "cuts2", "strands1", "strands2"],
                             constants={"numexpr":numexpr})

        readsMolecules = self.evaluate(
            "a = numexpr.evaluate('(chrms1 == chrms2) & (strands1 != strands2) &  (dist >=0) &"
            " (dist <= maximumMoleculeLength)')",
            internalVariables=["chrms1", "chrms2", "strands1", "strands2"],
            externalVariables={"dist": dist},
            constants={"maximumMoleculeLength": self.maximumMoleculeLength, "numexpr": numexpr})

        mask *= (readsMolecules == False)
        extraDE = mask.sum()
        self.metadata['210_ExtraDanglingReads'] = -extraDE + noSameFrag
        if mask.sum() == 0:
            raise Exception('No reads left after filtering. Please, check the input data')

        del dist, readsMolecules
        
        self.maskFilter(mask)
    
    def updateGenome(self, newGenome, oldGenome = 'current'):

        assert isinstance(newGenome, myGenome)
        
        newN = newGenome.chrmCount
        if oldGenome == "current":
            oldGenome = self.genome
        upgrade = newGenome.upgradeMatrix(oldGenome)
        if isinstance(oldGenome, myGenome):
            if oldGenome.hasEnzyme():
                newGenome.setEnzyme(oldGenome.enzymeName)
            oldGenome = oldGenome.idx2label
        
        chrms1 = np.array(self.chrms1, int)
        chrms2 = np.array(self.chrms2, int)

        if upgrade is not None:
            upgrade[upgrade == -1] = 9999  # to tell old SS reads from new SS reads

            chrms1 = upgrade[chrms1]
            self.chrms1 = chrms1
            del chrms1

            chrms2 = upgrade[chrms2]
            self.chrms2 = chrms2
            
        mask = ((self.chrms1 < newN) * (self.chrms2 < newN))
        self.genome = newGenome
        self.maskFilter(mask)
    
    def filterDuplicates(self):

        Nds = self.N

        # an array to determine unique rows. Eats 16 bytes per DS record
        dups = np.zeros((Nds, 2), dtype="int64", order="C")

        dups[:, 0] = self.chrms1
        dups[:, 0] *= self.fragIDmult
        dups[:, 0] += self.cuts1
        dups[:, 1] = self.chrms2
        dups[:, 1] *= self.fragIDmult
        dups[:, 1] += self.cuts2
        dups.sort(axis=1)
        dups.shape = (Nds * 2)
        strings = dups.view("|S16")
            # Converting two indices to a single string to run unique
        uids = uniqueIndex(strings)
        del strings, dups
        stay = np.zeros(Nds, bool)
        stay[uids] = True  # indexes of unique DS elements
        del uids
        
        self.metadata["310_DuplicatedRemoved"] = len(stay) - stay.sum()
        self.maskFilter(stay)
    
    def filterRsiteStart(self, offset=5):

        expression = "mask = (np.abs(dists1 - fraglens1) >= offset) * "\
        "((np.abs(dists2 - fraglens2) >= offset) )"
        
        mask = self.evaluate(expression,
                             internalVariables=["dists1", "fraglens1",
                                                "dists2", "fraglens2"],
                             constants={"offset": offset, "np": np},
                             outVariable=("mask", np.zeros(self.N, bool)))
                             
        self.metadata["320_StartNearRsiteReads"] = len(mask) - mask.sum()
        self.maskFilter(mask)
    
    def maskFilter(self, mask):
        
        # Uses 16 bytes per read
        length = 0
        ms = mask.sum()
        assert mask.dtype == np.bool
        self.N = ms
        for name in self.vectors:
            data = self._getData(name)
            ld = len(data)
            if length == 0:
                length = ld
            else:
                if ld != length:
                    self.delete()
            newdata = fasterBooleanIndexing(data, mask, outLen=ms,
                                        bounds=False)  # see mirnylib.numutils
            del data
            self._setData(name, newdata)
            del newdata
        del mask
        
        self.h5dict["metadata"] = self.metadata
    
    def merge(self, filenames):

        h5dicts = [h5dict(i, mode = 'r') for i in filenames]
        
        if all(["metadata" in i for i in h5dicts]):
            metadatas = [mydict["metadata"] for mydict in h5dicts]
            # print metadatas
            newMetadata = metadatas.pop()
            for oldData in metadatas:
                for key, value in oldData.items():
                    if (key in newMetadata):
                        newMetadata[key] += value
                    else:
                        log.warning('The key %s can not be found in some files',
                                    key)
            self.metadata = newMetadata
            self.h5dict["metadata"] = self.metadata

        for name in self.vectors.keys():
            res = []
            IfIn = [(name in mydict.keys()) for mydict in h5dicts]
            if not all(IfIn):
                continue
            for mydict in h5dicts:
                res.append(mydict[name])
            res = np.concatenate(res)
            self.N = len(res)
            self.DSnum = self.N
            self._setData(name, res)
            self.h5dict.flush()
            time.sleep(0.2)  # allow buffers to flush
        
        Types = ['LeftType', 'RightType', 'InnerType', 'OuterType']
        check = all([(i in j) for i in Types for j in h5dicts])
        if check:
            LeftType = np.zeros(50, dtype = int)
            RightType = np.zeros(50, dtype = int)
            InnerType = np.zeros(50, dtype = int)
            OuterType = np.zeros(50, dtype = int)
            for mydict in h5dicts:
                LeftType += mydict['LeftType']
                RightType += mydict['RightType']
                InnerType += mydict['InnerType']
                OuterType += mydict['OuterType']
        
            self.h5dict['LeftType'] = LeftType
            self.h5dict['RightType'] = RightType
            self.h5dict['InnerType'] = InnerType
            self.h5dict['OuterType'] = OuterType
        
    
    def printMetadata(self, saveTo):
        
        self.metadata = self.h5dict['metadata']
        
        Total = self.metadata['000_SequencedReads']
        Ureads = self.metadata['010_UniqueMappedReads']
        ligSeq = self.metadata['020_LigationCounts']
        selfLig = self.metadata['122_SelfLigationReads']
        longrange = self.metadata['412_IntraLongRangeReads(>=20Kb)']
        contacts = self.metadata['400_TotalContacts']
        
        Uratio = float(Ureads) / Total
        Lratio = float(ligSeq) / Total
        Fratio = float(selfLig) / Total
        longRatio = float(longrange) / contacts
        
        with open(saveTo, 'w') as myfile:
            for i in sorted(self.metadata):
                if (i[2] != '0'):
                    myfile.write('\t\t')
                elif (i[1] != '0') and (i[2] == '0'):
                    myfile.write('\t')
                myfile.write(str(i))
                myfile.write(':   ')
                myfile.write(str(self.metadata[i]))
                myfile.write('\n')
            myfile.write('\nCritical Indicators:\n')
            myfile.write('Unique-Mapping Ratio = %d / %d = %.4f\n' % (Ureads, Total, Uratio))
            myfile.write('Ligation-Junction Ratio = %d / %d = %.4f\n' % (ligSeq, Total, Lratio))
            myfile.write('Intra-Fragment Ratio = %d / %d = %.4f\n' % (selfLig, Total, Fratio))
            myfile.write('Long-Range Ratio = %d / %d = %.4f\n' % (longrange, contacts, longRatio))
                
    def saveHeatmap(self, filename, resolution, countDiagonalReads = 'Once'):

        try:
            os.remove(filename)
        except:
            pass

        tosave = h5dict(path = filename, mode = 'w')
        
        heatmap = self.buildAllHeatmap(resolution, countDiagonalReads)

        tosave['heatmap'] = heatmap
        
        del heatmap
        
        chromosomeStarts = np.array(self.genome.chrmStartsBinCont)
        numBins = self.genome.numBins
            
        tosave['resolution'] = resolution
        tosave['genomeBinNum'] = numBins
        tosave['genomeIdxToLabel'] = self.genome.idx2label
        tosave['chromosomeStarts'] = chromosomeStarts
    
    def saveByChromosomeHeatmap(self, filename, resolution = 40000,
                                includeTrans = False,
                                countDiagonalReads = "Once"):
        """
        Saves chromosome by chromosome heatmaps to h5dict.
        
        This method is not as memory demanding as saving all x all heatmap.

        Keys of the h5dict are of the format ["1 1"], where chromosomes are
        zero-based, and there is one space between numbers.

        Parameters
        ----------
        filename : str
            Filename of the h5dict with the output
            
        resolution : int
            Resolution to save heatmaps
            
        includeTrans : bool, optional
            Build inter-chromosomal heatmaps (default: False)
            
        countDiagonalReads : "once" or "twice"
            How many times to count reads in the diagonal bin

        """
        if countDiagonalReads.lower() not in ["once", "twice"]:
            raise ValueError("Bad value for countDiagonalReads")
            
        self.genome.setResolution(resolution)
        
        pos1 = self.evaluate("a = np.array(mids1 / {res}, dtype = 'int32')"
                             .format(res=resolution), "mids1")
        pos2 = self.evaluate("a = np.array(mids2 / {res}, dtype = 'int32')"
                             .format(res=resolution), "mids2")
                             
        chr1 = self.chrms1
        chr2 = self.chrms2
        
        # DS = self.DS  # 13 bytes per read up to now, 16 total
        mydict = h5dict(filename)

        for chrom in xrange(self.genome.chrmCount):
            if includeTrans == True:
                mask = ((chr1 == chrom) + (chr2 == chrom))
            else:
                mask = ((chr1 == chrom) * (chr2 == chrom))
            # Located chromosomes and positions of chromosomes
            c1, c2, p1, p2 = chr1[mask], chr2[mask], pos1[mask], pos2[mask]
            if includeTrans == True:
                # moving different chromosomes to c2
                # c1 == chrom now
                mask = (c2 == chrom) * (c1 != chrom)
                c1[mask], c2[mask], p1[mask], p2[mask] = c2[mask].copy(), c1[
                    mask].copy(), p2[mask].copy(), p1[mask].copy()
                del c1  # ignore c1
                args = np.argsort(c2)
                c2 = c2[args]
                p1 = p1[args]
                p2 = p2[args]

            for chrom2 in xrange(chrom, self.genome.chrmCount):
                if (includeTrans == False) and (chrom2 != chrom):
                    continue
                start = np.searchsorted(c2, chrom2, "left")
                end = np.searchsorted(c2, chrom2, "right")
                cur1 = p1[start:end]
                cur2 = p2[start:end]
                label = np.asarray(cur1, "int64")
                label *= self.genome.chrmLensBin[chrom2]
                label += cur2
                maxLabel = self.genome.chrmLensBin[chrom] * \
                           self.genome.chrmLensBin[chrom2]
                counts = np.bincount(label, minlength = maxLabel)
                assert len(counts) == maxLabel
                mymap = counts.reshape((self.genome.chrmLensBin[chrom], -1))
                if chrom == chrom2:
                    mymap = mymap + mymap.T
                    if countDiagonalReads.lower() == "once":
                        fillDiagonal(mymap, np.diag(mymap).copy() / 2)
                mydict["%d %d" % (chrom, chrom2)] = mymap
        
        mydict['resolution'] = resolution

        return
    
    def buildAllHeatmap(self, resolution, countDiagonalReads = 'Once'):
        
        # 8 bytes per record + heatmap
        self.genome.setResolution(resolution)
        numBins = self.genome.numBins
        label = self.genome.chrmStartsBinCont[self.chrms1]
        label = np.asarray(label, dtype="int64")
        label += self.mids1 / resolution
        label *= numBins
        label += self.genome.chrmStartsBinCont[self.chrms2]
        label += self.mids2 / resolution
        
        counts = np.bincount(label, minlength=numBins ** 2)
        if len(counts) > numBins ** 2:
            raise StandardError("\nHeatMap exceed length of the genome!")

        counts.shape = (numBins, numBins)
        for i in xrange(len(counts)):
            counts[i, i:] += counts[i:, i]
            counts[i:, i] = counts[i, i:]
        if countDiagonalReads.lower() == "once":
            diag = np.diag(counts)
            fillDiagonal(counts, diag / 2)
        elif countDiagonalReads.lower() == "twice":
            pass
        else:
            raise ValueError("Bad value for countDiagonalReads")
            
        return counts
    
    def typePlot(self, outfile, dpi = 500):
        
        Keys = ['LeftType', 'RightType', 'InnerType', 'OuterType']
        if not all([(i in self.h5dict for i in Keys)]):
            raise StandardError
        
        LeftType = self.h5dict['LeftType'][:25]
        RightType = self.h5dict['RightType'][:25]
        InnerType = self.h5dict['InnerType'][:25]
        OuterType = self.h5dict['OuterType'][:25]
        
        Total = LeftType + RightType + InnerType + OuterType
        Total = Total.astype(np.float)
        
        LeftRatio = LeftType / Total
        RightRatio = RightType / Total
        InnerRatio = InnerType / Total
        OuterRatio = OuterType / Total
        
        fig = plt.figure(figsize = (15, 9))
        ax = fig.add_subplot(111)
        lines = []
        labels = ['Left Type', 'Right Type', 'Inner Type', 'Outer Type']
        x = np.arange(1, LeftRatio.size + 1)
        idx = 0
        for y in [LeftRatio, RightRatio, InnerRatio, OuterRatio]:
            L = ax.plot(x, y, color = colorPool[idx], linewidth = 2)
            lines.extend(L)
            idx += 1
        
        ax.set_xlabel('Genomic Separation, KB')
        ax.set_ylabel('Type Ratio')
        
        ax.set_ylim((0, 1))
        ax.set_title('Read Pair Type Statistics')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        ax.legend(lines, labels, frameon = False, fontsize = 11, labelspacing = 1,
                  handletextpad = 1, borderpad = 1, markerscale = 1, numpoints = 1,
                  ncol = 2, loc = 'upper right')
        
        plt.savefig(outfile, dpi = dpi)
        plt.close()
        
        
    def __setattr__(self, x, value):
        
        if x == 'vectors':
            return object.__setattr__(self, x, value)

        if x in self.vectors.keys():
            self._setData(x, value)
        else:
            return object.__setattr__(self, x, value)
    
    def __getattribute__(self, x):
        
        if x == 'vectors':
            return object.__getattribute__(self, x)

        if x in self.vectors.keys():
            a = self._getData(x)
            return a
        else:
            return object.__getattribute__(self, x)
    
    def _getData(self, name):
        
        return self.h5dict[name]
    
    def _setData(self, name, data):
        
        dtype = np.dtype(self.vectors[name])
        data = np.asarray(data, dtype=dtype)
        self.h5dict[name] = data

class cHiResHiC(HiResHiC):
    
    def export(self, filename, mode = 'cis'):
        mydict = h5dict(filename)
        if mode == 'cis':
            for i in self.cisKeys:
                data = self.data[i].getData()
                mydict["%d %d" % i] = data
        else:
            for i in self.allKeys:
                data = self.data[i].getData()
                mydict["%d %d" % i] = data
        mydict["resolution"] = self.resolution