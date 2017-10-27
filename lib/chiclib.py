# Created on Mon May 25 16:57:52 2015

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

import logging, os, time, re, tempfile, gc
from textwrap import dedent
import numexpr, matplotlib
matplotlib.use('Agg')
import numpy as np
from mirnylib.genome import Genome
from hiclib.fragmentHiC import HiCdataset
from hiclib.hicShared import mydtype, mydtypeSorter, searchsorted, h5dictBinarySearch
from hiclib.binnedData import binnedData
from mirnylib import numutils
from mirnylib.numutils import uniqueIndex, fillDiagonal, externalMergeSort, completeIC
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

colorPool = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#348ABD', '#A60628']

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
    
    def __init__(self, filename, genome, maximumMoleculeLength=500,
                 mode="a", tmpFolder="/tmp", dictToStoreIDs="h5dict",
                 compression="gzip", compression_opts=3):
        
        self.vectors = {
            # chromosomes for each read.
            "chrms1": "int16", "chrms2": "int16",
            #strand to which the read maps
            "strands1": "bool", "strands2": "bool",
            #Start of the read ("ultrasonic" cut site)
            "cuts1": "int32", "cuts2": "int32"}
        
        self.vectors2 = {
            # fragment lengthes
            "fraglens1": "int32", "fraglens2": "int32",
             # fragid as defined in the manual
            "fragids1": "int64", "fragids2": "int64",
            # midpoint of a fragment, determined as "(start+end)/2"
            "mids1": "int32", "mids2": "int32",
            # distance from a cut site to the restriction fragment
            "dists1": "int32", "dists2": "int32",
             # distance between fragments. If -1, different chromosomes.
             # If -2, different arms.
            "distances": "int32",
            }
        self.vectors3 = {
            # absolute ID of a restriction fragment (0,1,2,... numFragments)
            "rfragAbsIdxs1": "int32", "rfragAbsIdxs2": "int32"}
        
        if dictToStoreIDs == "dict":
            self.rfragIDDict = {}
        elif dictToStoreIDs == "h5dict":
            self.rfragIDDict = h5dict()
        else:
            self.rfragIDDict = dictToStoreIDs
            
        self.metadata = {}
        self.tmpDir = tmpFolder
        if not os.path.exists(self.tmpDir):
            os.makedirs(self.tmpDir)

        #-------Initialization of the genome and parameters-----
        self.mode = mode
        self.genome = genome

        self.chromosomeCount = self.genome.chrmCount
        self.fragIDmult = self.genome.fragIDmult  # used for building heatmaps
        self.rFragIDs = self.genome.rfragMidIds
        self.rFragLens = np.concatenate(self.genome.rfragLens)
        self.rFragMids = np.concatenate(self.genome.rfragMids)
        self.rsites = self.genome.rsiteIds
        # to speed up searchsorted we use positive-only numbers
        self.rsitesPositive = self.rsites + 2*self.fragIDmult

        self.maximumMoleculeLength = maximumMoleculeLength
        
        # File to save the data
        self.filename = os.path.abspath(os.path.expanduser(filename))
        # Chunk size for h5dict operation, external sorting, etc.
        self.chunksize = 10000000

        self.h5dict = h5dict(self.filename, mode=mode, in_memory=False)
        self.h5dict.setCompression(compression, compression_opts)
        
        if "chrms1" in list(self.h5dict.keys()):
            self.N = len(self.h5dict.get_dataset("chrms1"))
        if "metadata" in self.h5dict:
            self.metadata = self.h5dict["metadata"]
    
    def parseInputData(self, dictLike, **kwargs):
        
        if not os.path.exists(dictLike):
            raise IOError('File not found: %s' % dictLike)
        
        dictLike = h5dict(dictLike, 'r')
        
        self.chrms1 = dictLike['chrms1']
        self.chrms2 = dictLike['chrms2']
        self.cuts1 = dictLike['cuts1']
        self.cuts2 = dictLike['cuts2']
        self.strands1 = dictLike['strands1']
        self.strands2 = dictLike['strands2']
        
        self.N = len(self.chrms1)
        
        dictLike['misc']['genome']['idx2label']
        self.updateGenome(self.genome,
                          oldGenome = dictLike["misc"]["genome"]["idx2label"])
        
        DSmask = (self.chrms1 >= 0) * (self.chrms2 >= 0)
        self.metadata['100_DoubleUniqueMapped'] = DSmask.sum()
        
        # Discard dangling ends and self-circles
        sameFragMask = self.evaluate("a = (rfragAbsIdxs1 == rfragAbsIdxs2)", 
                                     ["rfragAbsIdxs1", "rfragAbsIdxs2"]) * DSmask
        cutDifs = self.cuts2[sameFragMask] > self.cuts1[sameFragMask]
        s1 = self.strands1[sameFragMask]
        s2 = self.strands2[sameFragMask]
        SSDE = (s1 != s2)
        SS = SSDE * (cutDifs == s2)
        DM = SSDE & np.logical_not(SS)
        SS_N = SS.sum()
        DM_N = DM.sum()
        SSDE_N = SSDE.sum()
        sameFrag_N = sameFragMask.sum()
        self.metadata['120_SameFragmentReads'] = sameFrag_N
        self.metadata['122_SelfLigationReads'] = SS_N
        self.metadata['124_DanglingReads'] = DM_N
        self.metadata['126_UnknownMechanism'] = sameFrag_N - SSDE_N
        mask = DSmask * (~sameFragMask)
        noSameFrag = mask.sum()
        
        del DSmask
        
        ## Dissect dangling ends ...
        Ddists1 = self.fraglens1 - self.dists1
        Ddists2 = self.fraglens2 - self.dists2
        extD1 = Ddists1[sameFragMask][DM]
        extD2 = Ddists2[sameFragMask][DM]
        extLen = self.fraglens1[sameFragMask][DM].astype(float)
        extD = np.r_['0,2', extD1, extD2].min(0)
        
        del Ddists1, Ddists2, extD1, extD2
        
        # distance between sites facing each other
        dist = self.evaluate("a = numexpr.evaluate('- cuts1 * (2 * strands1 -1) - "
                             "cuts2 * (2 * strands2 - 1)')",
                             ["cuts1", "cuts2", "strands1", "strands2"],
                             constants={"numexpr":numexpr})
        extSpace = dist[sameFragMask][DM]
        
        self.h5dict['_DanglingDetials'] = {'extLen':extLen,
                                           'extD':extD,
                                           'extSpace':extSpace}
        
        del sameFragMask

        readsMolecules = self.evaluate(
            "a = numexpr.evaluate('(chrms1==chrms2) & (strands1!=strands2) & (dist>=0) &"
            " (dist<=maximumMoleculeLength)')",
            internalVariables=["chrms1", "chrms2", "strands1", "strands2"],
            externalVariables={"dist": dist},
            constants={"maximumMoleculeLength": self.maximumMoleculeLength, "numexpr": numexpr})
        
        mask *= (readsMolecules==False)
        extraDE = mask.sum()
        self.metadata['210_ExtraDanglingReads'] = -extraDE + noSameFrag
        if mask.sum() == 0:
            raise Exception('No reads left after filtering. Please check the input data')

        del dist, readsMolecules
        
        self.maskFilter(mask)
    
    def updateGenome(self, newGenome, oldGenome='current'):

        assert isinstance(newGenome, myGenome)
        
        newN = newGenome.chrmCount
        if oldGenome == "current":
            oldGenome = self.genome
        upgrade = newGenome.upgradeMatrix(oldGenome)
        if isinstance(oldGenome, myGenome):
            if (oldGenome.hasEnzyme()==True) and (newGenome.hasEnzyme()==False):
                newGenome.setEnzyme(oldGenome.enzymeName)
        
        chrms1 = np.array(self.chrms1, int)
        chrms2 = np.array(self.chrms2, int)

        if upgrade is not None:
            upgrade[upgrade == -1] = 9999  # to tell old SS reads from new SS reads

            chrms1 = upgrade[chrms1]
            self.chrms1 = chrms1
            del chrms1

            chrms2 = upgrade[chrms2]
            self.chrms2 = chrms2
            del chrms2
            
        mask = ((self.chrms1 < newN) * (self.chrms2 < newN))
        self.genome = newGenome
        self.maskFilter(mask)
        
    def make_tempfile(self):
        
        tl = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
        kw = {'suffix':tl, 'dir':self.tmpDir}
        fd, tmpfil = tempfile.mkstemp(**kw)
        os.close(fd)
        
        return tmpfil
        
    def filterDuplicates(self, chunkSize=50000000):
        
        if self.N > 200000000:
            mode = "hdd"
        else:
            mode = "ram"
        
        if mode == "ram":
            dups = np.zeros((self.N, 2), dtype="int64", order="C")
            dups[:, 0] = self.chrms1
            dups[:, 0] *= self.fragIDmult
            dups[:, 0] += self.cuts1
            dups[:, 1] = self.chrms2
            dups[:, 1] *= self.fragIDmult
            dups[:, 1] += self.cuts2
            dups.sort(axis=1)
            dups.shape = (self.N * 2)
            strings = dups.view("|S16")
            # Converting two indices to a single string to run unique
            uids = uniqueIndex(strings)
            del strings, dups
            stay = np.zeros(self.N, bool)
            stay[uids] = True  # indexes of unique DS elements
            del uids
        else:
            tmpfil = self.make_tempfile()
            a = h5dict(tmpfil)
            a.add_empty_dataset("duplicates", (self.N,), dtype="|S24")
            a.add_empty_dataset("temp", (self.N,), dtype="|S24")
            dset = a.get_dataset("duplicates")
            tempdset = a.get_dataset("temp")
            code = dedent("""
            tmp = np.array(chrms1, dtype=np.int64) * fragIDmult + cuts1
            tmp2 = np.array(chrms2, dtype=np.int64) * fragIDmult + cuts2
            newarray = np.zeros((len(tmp),3), dtype = np.int64)
            newarray[:,0] = tmp
            newarray[:,1] = tmp2
            newarray[:,:2].sort(axis=1)
            newarray[:,2] = np.arange(start, end, dtype=np.int64)
            newarray.shape = (3*len(tmp))
            a = np.array(newarray.view("|S24"))
            assert len(a) == len(chrms1)
            """)
            self.evaluate(code, ["chrms1", "cuts1", "chrms2", "cuts2"],
                          constants={"np":np, "fragIDmult":self.fragIDmult},
                          outVariable=("a", dset))
            stay = np.zeros(self.N, bool)
            numutils.externalMergeSort(dset, tempdset, chunkSize=chunkSize)
            bins = list(range(0, self.N-1000, self.chunksize)) + [self.N-1]
            for start, end in zip(bins[:-1], bins[1:]):
                curset = dset[start:end+1]
                curset = curset.view(dtype=np.int64)
                curset.shape = (len(curset)//3, 3)
                unique = (curset[:-1,0] != curset[1:,0]) + (curset[:-1,1] != curset[1:,1])
                stay[curset[:,2][unique]] = True
                if end == self.N-1:
                    stay[curset[-1,2]] = True
            del a
            os.remove(tmpfil)
            
        self.metadata["310_DuplicatedRemoved"] = len(stay) - stay.sum()
        self.maskFilter(stay)
    
    def maskFilter(self, mask):
        
        # Uses 16 bytes per read
        for i in list(self.rfragIDDict.keys()):
            del self.rfragIDDict[i]
        
        length = 0
        ms = mask.sum()
        assert mask.dtype == np.bool
        self.N = ms
        if ms == 0:
            return

        for name in self.vectors:
            data = self._getData(name)
            ld = len(data)
            if length == 0:
                length = ld
            else:
                if ld != length:
                    self.delete()
            newdata = data[mask]
            del data
            self._setData(name, newdata)
            del newdata
        del mask
        
        self.h5dict['metadata'] = self.metadata
    
    def _sortData(self):

        if not hasattr(self, "dataSorted"):
            tmpfil = self.make_tempfile()
            mydict = h5dict(tmpfil, 'w')
            data = mydict.add_empty_dataset("sortedData", (self.N,), mydtype)
            tmp = mydict.add_empty_dataset("trash", (self.N,), mydtype)
            code = dedent("""
            a = np.empty(len(chrms1), dtype = mydtype)
            mask = (chrms1 > chrms2) | ( (chrms1 == chrms2) & (cuts1 > cuts2))

            chrms2[mask],chrms1[mask] = chrms1[mask].copy(), chrms2[mask].copy()
            cuts1[mask],cuts2[mask] = cuts2[mask].copy(), cuts1[mask].copy()
            strands1[mask],strands2[mask] = strands2[mask].copy(),strands1[mask].copy()

            a["chrms1"] = chrms1
            a["pos1"] = cuts1
            a["chrms2"] = chrms2
            a["pos2"] = cuts2
            a["strands1"] = strands1
            a["strands2"] = strands2
            """)
            self.evaluate(expression=code, internalVariables=["chrms1","chrms2","cuts1","cuts2","strands1","strands2"],
                          constants={"np":np, "mydtype":mydtype}, outVariable=("a", data))

            externalMergeSort(data, tmp, sorter=mydtypeSorter, searchsorted=searchsorted,
                              chunkSize=max(150000000, self.chunksize))
            sdata = mydict.get_dataset("sortedData")

            c1 = self.h5dict.get_dataset("chrms1")
            c2 = self.h5dict.get_dataset("chrms2")
            p1 = self.h5dict.get_dataset("cuts1")
            p2 = self.h5dict.get_dataset("cuts2")
            s1 = self.h5dict.get_dataset("strands1")
            s2 = self.h5dict.get_dataset("strands2")

            for start, end in self._getChunks():
                data = sdata[start:end]
                c1[start:end] = data["chrms1"]
                c2[start:end] = data["chrms2"]
                p1[start:end] = data["pos1"]
                p2[start:end] = data["pos2"]
                s1[start:end] = data["strands1"]
                s2[start:end] = data["strands2"]
            self.dataSorted = True
            del mydict
            os.remove(tmpfil)
            gc.collect()
    
    def _getChunks(self, chunkSize="default"):
        
        if chunkSize == "default":
            chunkSize = self.chunksize
        if chunkSize > 0.5 * self.N:
            return [(0, self.N)]
        points = list(range(0,self.N-chunkSize//2,chunkSize)) + [self.N]
        
        return list(zip(points[:-1], points[1:]))
    
    def merge(self, filenames):

        h5dicts = [h5dict(i, mode='r') for i in filenames]
        
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
        
        self.N = sum([len(i.get_dataset("strands1")) for i in h5dicts])

        for name in self.vectors.keys():
            if name in self.h5dict:
                del self.h5dict[name]
            self.h5dict.add_empty_dataset(name, (self.N,), self.vectors[name])
            dset = self.h5dict.get_dataset(name)
            position = 0
            for mydict in h5dicts:
                cur = mydict[name]
                dset[position:position+len(cur)] = cur
                position += len(cur)
            self.h5dict.flush()
            time.sleep(0.2)  # allow buffers to flush
        
        self._sortData()
        
        stats = ['extLen', 'extD', 'extSpace']
        pool = {}
        check = all([('_DanglingDetials' in h) for h in h5dicts])
        if check:
            for name in stats:
                res = []
                for mydict in h5dicts:
                    res.append(mydict['_DanglingDetials'][name])
                res = np.concatenate(res)
                pool[name] = res
            self.h5dict['_DanglingDetials'] = pool
        
        Types = ['LeftType', 'RightType', 'InnerType', 'OuterType']
        pool = {}
        check = all([('_DirectionTypeStats' in h) for h in h5dicts])
        if check:
            for Type in Types:
                tmp = np.zeros(50, dtype=int)
                for mydict in h5dicts:
                    tmp += mydict['_DirectionTypeStats'][Type]
                pool[Type] = tmp
            self.h5dict['_DirectionTypeStats'] = pool
        
        tempDict = {}
        for mydict in h5dicts:
            if 'genomeInformation' in mydict:
                tempDict.update(mydict['genomeInformation'])
        self.h5dict['genomeInformation'] = tempDict
        
    
    def printMetadata(self, saveTo):
        
        self.metadata = self.h5dict['metadata']
        
        longrange = self.metadata['412_IntraLongRangeReads(>=20Kb)']
        contacts = self.metadata['400_TotalContacts']
        longRatio = float(longrange) / contacts
                         
        fromMapping = ['000_SequencedReads', '010_UniqueMappedReads',
                       '020_LigationCounts']
        
        check = all([m in self.metadata for m in fromMapping])
        if check:
            Total = self.metadata['000_SequencedReads']
            Ureads = self.metadata['100_DoubleUniqueMapped']
            ligSeq = self.metadata['030_LigationCounts']
            selfLig = self.metadata['122_SelfLigationReads']
            dangling = self.metadata['124_DanglingReads']
            Uratio = float(Ureads) / Total
            Lratio = float(ligSeq) / Total
            Fratio = float(selfLig) / Total
            Dratio = float(dangling) / Total
            usage = float(contacts) / Total
        ##-----------------------------------------------------------------
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
            if check:
                myfile.write('Double Unique Mapped Ratio = %d / %d = %.4f\n' % (Ureads, Total, Uratio))
                myfile.write('Ligation-Junction Ratio = %d / %d = %.4f\n' % (ligSeq, Total, Lratio))
                myfile.write('Self-Ligation Ratio = %d / %d = %.4f\n' % (selfLig, Total, Fratio))
                myfile.write('Dangling-Reads Ratio = %d / %d = %.4f\n' % (dangling, Total, Dratio))
            myfile.write('Long-Range Ratio = %d / %d = %.4f\n' % (longrange, contacts, longRatio))
            if check:
                myfile.write('Data Usage = %d / %d = %.4f\n' % (contacts, Total, usage))
                
    def saveHeatmap(self, filename, resolution, gInfo):

        try:
            os.remove(filename)
        except:
            pass

        tosave = h5dict(path=filename, mode='w')
        
        heatmap = self.buildAllHeatmap(resolution)

        tosave['heatmap'] = heatmap
        
        del heatmap
        
        chromosomeStarts = np.array(self.genome.chrmStartsBinCont)
            
        tosave['resolution'] = resolution
        tosave['chromosomeStarts'] = chromosomeStarts
        tosave['genomeInformation'] = gInfo
    
    def saveByChromosomeHeatmap(self, filename, resolution, gInfo,
                                includeTrans=False):
        
        self.genome.setResolution(resolution)

        mydict = h5dict(filename)

        for chrom in range(self.genome.chrmCount):
            c1 = self.h5dict.get_dataset("chrms1")
            p1 = self.h5dict.get_dataset("cuts1")
            low = h5dictBinarySearch(c1, p1, (chrom,-1), "left")
            high = h5dictBinarySearch(c1, p1, (chrom,999999999), "right")

            chr1 = self._getVector("chrms1", low, high)
            chr2 = self._getVector("chrms2", low, high)
            pos1 = np.array(self._getVector("mids1",low,high)//resolution,
                            dtype=np.int32)
            pos2 = np.array(self._getVector("mids2",low,high)//resolution,
                            dtype=np.int32)

            assert (chr1==chrom).all()  # getting sure that bincount worked

            args = np.argsort(chr2)
            chr2 = chr2[args]
            pos1 = pos1[args]
            pos2 = pos2[args]

            for chrom2 in range(chrom, self.genome.chrmCount):
                if (includeTrans == False) and (chrom2 != chrom):
                    continue
                start = np.searchsorted(chr2, chrom2, "left")
                end = np.searchsorted(chr2, chrom2, "right")
                cur1 = pos1[start:end]
                cur2 = pos2[start:end]
                label = np.array(cur1, "int64")
                label *= self.genome.chrmLensBin[chrom2]
                label += cur2
                maxLabel = self.genome.chrmLensBin[chrom] * \
                           self.genome.chrmLensBin[chrom2]
                counts = np.bincount(label, minlength=maxLabel)
                mymap = counts.reshape((self.genome.chrmLensBin[chrom], -1))
                if chrom==chrom2:
                    mymap = mymap + mymap.T
                    fillDiagonal(mymap, np.diag(mymap).copy()/2)
                mydict["%d %d" % (chrom, chrom2)] = mymap
        
        mydict['resolution'] = resolution
        mydict['genomeInformation'] = gInfo

        return
    
    def buildAllHeatmap(self, resolution):
        
        for start, end in self._getChunks(30000000):
            # 8 bytes per record + heatmap
            self.genome.setResolution(resolution)
            numBins = self.genome.numBins
            label = self.genome.chrmStartsBinCont[self._getVector("chrms1",start,end)]
            label = np.asarray(label,dtype="int64")
            label += (self._getVector("mids1",start,end)//resolution).astype(np.int64)
            label *= numBins
            label += self.genome.chrmStartsBinCont[self._getVector("chrms2",start,end)]
            label += (self._getVector("mids2",start,end)//resolution).astype(np.int64)
            counts = np.bincount(label, minlength=numBins**2)
            if len(counts) > numBins ** 2:
                raise StandardError("\nHeatMap exceed length of the genome!")

            counts.shape = (numBins, numBins)
            try:
                heatmap += counts  # @UndefinedVariable
            except:
                heatmap = counts

        for i in range(len(heatmap)):
            heatmap[i,i:] += heatmap[i:,i]
            heatmap[i:,i] = heatmap[i,i:]
        diag = np.diag(heatmap)
        fillDiagonal(heatmap, diag/2)
        
        return heatmap
    
    def typePlot(self, outfile, dpi = 500):
        
        LT = self.h5dict['_DirectionTypeStats']['LeftType'][:25]
        RT = self.h5dict['_DirectionTypeStats']['RightType'][:25]
        IT = self.h5dict['_DirectionTypeStats']['InnerType'][:25]
        OT = self.h5dict['_DirectionTypeStats']['OuterType'][:25]
        
        Total = LT + RT + IT + OT
        Total = Total.astype(np.float)
        
        LR = LT / Total
        RR = RT / Total
        IR = IT / Total
        OR = OT / Total
        
        fig = plt.figure(figsize = (15, 9))
        ax = fig.add_subplot(111)
        lines = []
        labels = ['Left Type', 'Right Type', 'Inner Type', 'Outer Type']
        x = np.arange(1, LR.size + 1)
        idx = 0
        for y in [LR, RR, IR, OR]:
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
    
    def dangStats(self, prefix, offset = 2, dpi = 500):
        
        extD = self.h5dict['_DanglingDetials']['extD']
        extLen = self.h5dict['_DanglingDetials']['extLen']
        extSpace = self.h5dict['_DanglingDetials']['extSpace']
        
        m = int(np.percentile(extD[extD<0],0.5))
        extLen = extLen[extD>=m]
        extSpace = extSpace[extD>=m]
        extD = extD[extD>=m]
        if offset < -m:
            offset = -m
        
        extD += offset
        extR = extD / extLen
        
        lenlow = np.percentile(extLen, 0.1)
        lenhigh = np.percentile(extLen, 99.9)
        spacelow = np.percentile(extSpace, 0.1)
        spacehigh = np.percentile(extSpace, 99.9)
        
        lenM = (extLen >= lenlow) & (extLen <= lenhigh)
        spaceM = (extSpace >= spacelow) & (extSpace <= spacehigh)
        
        fig = plt.figure(figsize = (15, 9))
        ax = fig.add_subplot(111)
        ax.hist(extR[lenM], bins = 15, color = colorPool[-2])
        ax.set_title('Start Site of Dandling Ends Relative to Fragments (Ratio)')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        plt.savefig('-'.join([prefix, 'danglingStart']) + '.png', dpi = dpi)
        plt.close()
        
        fig = plt.figure(figsize = (15, 9))
        ax = fig.add_subplot(111)
        ax.hist(extSpace[spaceM], bins = 15, color = colorPool[-1])
        ax.set_title('Estimated Library Size Distribution (bp)')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        plt.savefig('-'.join([prefix, 'librarySize']) + '.png', dpi = dpi)
        plt.close()

class cBinnedData(binnedData):
    
    def export(self, name, outFilename):

        if not name in self.dataDict:
            raise ValueError("No data {name}".format(name=name))
            
        toexport = {}
        toexport["heatmap"] = self.dataDict[name]
        toexport["resolution"] = self.resolution
        toexport["chromosomeStarts"] = self.chromosomeStarts
        myh5dict = h5dict(outFilename, mode="w")
        myh5dict.update(toexport)

class HiResHiC(object):
    
    def __init__(self, genome, resolution, raw):
        
        self.genome = genome
        self.resolution = resolution
        self.rawdata = h5dict(raw, 'r')
        self.cisKeys = ['{0} {0}'.format(i) for i in self.genome.idx2label]
    
    def iterativeCorrection(self, outname):
        
        mydict = h5dict(outname)
        for key in self.cisKeys:
            bychr = self.rawdata[key]
            corrected = completeIC(bychr, returnBias=False)
            mydict[key] = corrected
        mydict["resolution"] = self.resolution