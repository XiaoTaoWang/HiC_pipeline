# Created on Tue Dec 23 21:15:19 2014

# Author: XiaoTao Wang
# Organization: HuaZhong Agricultural University

import logging
import numpy as np
from hiclib.fragmentHiC import HiCdataset
from mirnylib.numutils import uniqueIndex

log = logging.getLogger(__name__)

# A customized HiCdataset class, which makes filtering processes more flexible
class cHiCdataset(HiCdataset):
    
    def parseInputData(self, dictLike, commandArgs, **kwargs):
        '''
        Added Parameters
        ----------------
        commandArgs : NameSpace
            A NameSpace object defined by argparse.            
        '''
        ## Necessary Modules
        import numexpr
        
        # Simply load merged data
        log.log(21, 'Loading data ...')
        self.merge([dictLike])
        log.log(21, 'Done!')
        
        log.log(21, 'Basic statistics on your data:')
        
        log.log(21, 'Total Reads ...')
        # Total Reads
        self.trackLen = len(self.chrms1)
        log.log(21, self.trackLen)
        
        self.metadata["100_TotalReads"] = self.trackLen
        self.metadata["152_removedUnusedChromosomes"] = self.trackLen - self.N
        self.metadata["150_ReadsWithoutUnusedChromosomes"] = self.N
        
        log.log(21, 'Total DS Reads ...')
        DSmask = (self.chrms1 >= 0) * (self.chrms2 >= 0)
        self.metadata["200_totalDSReads"] = DSmask.sum()
        log.log(21, self.metadata["200_totalDSReads"])
        
        self.metadata["201_DS+SS"] = len(DSmask)
        self.metadata["202_SSReadsRemoved"] = len(DSmask) - DSmask.sum()
        
        mask = DSmask
        
        log.log(21, 'Determining Hi-C library size from dangling ends ...')
        
        ## Information based on restriction fragments
        sameFragMask = self.evaluate("a = (fragids1 == fragids2)",
                                     ["fragids1", "fragids2"]) * DSmask
        cutDifs = self.cuts2[sameFragMask] > self.cuts1[sameFragMask]
        s1 = self.strands1[sameFragMask]
        s2 = self.strands2[sameFragMask]
        SSDE = (s1 != s2)
        SS = SSDE * (cutDifs == s2)
        Dangling = SSDE & (~SS)
        SS_N = SS.sum()
        SSDE_N = SSDE.sum()
        sameFrag_N = sameFragMask.sum()
        
        dist = self.evaluate("a = - cuts1 * (2 * strands1 -1) - "
                             "cuts2 * (2 * strands2 - 1)",
                             ["cuts1", "cuts2", "strands1", "strands2"])
        Dangling_L = dist[sameFragMask][Dangling]
        library_L = int(np.ceil((np.percentile(Dangling_L, 95))))
        self.maximumMoleculeLength = library_L
        
        log.log(21, library_L)
        
        readsMolecules = self.evaluate(
            "a = numexpr.evaluate('(chrms1 == chrms2) & (strands1 != strands2) &  (dist >=0) &"
            " (dist <= maximumMoleculeLength)')",
            internalVariables=["chrms1", "chrms2", "strands1", "strands2"],
            externalVariables={"dist": dist},
            constants={"maximumMoleculeLength": self.maximumMoleculeLength, "numexpr": numexpr})
        
        if commandArgs.sameFragments:
            log.log(21, 'Removing read pairs located in the same restriction fragments ...')            
            mask *= (-sameFragMask)
            noSameFrag = mask.sum()
            self.metadata["210_sameFragmentReadsRemoved"] = sameFrag_N
            self.metadata["212_Self-Circles"] = SS_N
            self.metadata["214_DandlingEnds"] = SSDE_N - SS_N
            self.metadata["216_error"] = sameFrag_N - SSDE_N
            mask *= (readsMolecules == False)
            extraDE = mask.sum()
            self.metadata["220_extraDandlingEndsRemoved"] = -extraDE + noSameFrag
            log.log(21, 'Done!')
            log.log(21, '%s reads are remained', extraDE)
            
        if commandArgs.RandomBreaks:
            log.log(21, 'Removing "Random Breaks" ...')
            
            ini_N = extraDE
            mask *= ((self.dists1 + self.dists2) <= library_L)
            rb_N = ini_N - mask.sum()
            self.metadata["330_removeRandomBreaks"] = rb_N
            
            log.log(21, 'Done!')
            log.log(21, '%s reads are remained', (ini_N - rb_N))
        
        if mask.sum() == 0:
            raise Exception(
                'No reads left after filtering. Please, check the input data')
            
        del DSmask, sameFragMask
        del dist, readsMolecules
        
        self.metadata["300_ValidPairs"] = self.N
        
        self.maskFilter(mask)
    
    def filterDuplicates(self):
        
        log.log(21, 'Filtering duplicates in DS reads ...')
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
        uflen = len(self.ufragments)
        Remained_N = stay.sum()
        self.metadata["320_duplicatesRemoved"] = len(stay) - Remained_N
        self.maskFilter(stay)
        assert len(self.ufragments) == uflen  # self-check
        
        log.log(21, 'Done!')
        log.log(21, '%s reads are remained', Remained_N)
    
    def filterRsiteStart(self, offset = 5):
        """
        Removes reads that start within x bp near rsite

        Parameters
        ----------

        offset : int
            Number of bp to exclude next to rsite, not including offset

        """
        log.log(21, 'Filtering reads starting within 5 bp near the restriction site ...')

        expression = "mask = (np.abs(dists1 - fraglens1) >= offset) * "\
        "((np.abs(dists2 - fraglens2) >= offset) )"
        mask = self.evaluate(expression,
                             internalVariables=["dists1", "fraglens1",
                                                "dists2", "fraglens2"],
                             constants={"offset": offset, "np": np},
                             outVariable=("mask", np.zeros(self.N, bool)))
        Remained_N = mask.sum()
        self.metadata["310_startNearRsiteRemoved"] = len(mask) - Remained_N
        self.maskFilter(mask)
        
        log.log(21, 'Done!')
        log.log(21, '%s reads are remained', Remained_N)
        
    def filterLarge(self, cutlarge = 100000, cutsmall = 100):
        """
        Removes very large and small fragments.

        Parameters
        ----------
        cutlarge : int
            Remove fragments larger than it
        cutsmall : int
            Remove fragments smaller than it
        """
        self._buildFragments()
        
        log.log(21, 'Removing too large and too small fragments ...')
        
        p = (self.ufragmentlen < (cutlarge)) * (self.ufragmentlen > cutsmall)
        N1 = self.N
        self.fragmentFilter(self.ufragments[p])
        N2 = self.N
        self.metadata["340_removedLargeSmallFragments"] = N1 - N2
        self._dumpMetadata()
        
        log.log(21, 'Done!')
        log.log(21, '%s reads are remained', N2)
    
    def filterExtreme(self, cutH = 0.005, cutL = 0):
        """
        Removes fragments with most and/or least # counts

        Parameters
        ----------
        cutH : float, 0<=cutH < 1, optional
            Fraction of the most-counts fragments to be removed
            
        cutL : float, 0<=cutL<1, optional
            Fraction of the least-counts fragments to be removed
        """
        self._buildFragments()
        
        log.log(21, 'Removing the top 0.5% fragments with the greatest number of reads ...')
        
        s = self.fragmentSum()
        ss = np.sort(s)

        valueL, valueH = np.percentile(ss, [100. * cutL, 100 * (1. - cutH)])
        news = (s >= valueL) * (s <= valueH)
        N1 = self.N
        self.fragmentFilter(self.ufragments[news])
        self.metadata["350_removedFromExtremeFragments"] = N1 - self.N
        self._dumpMetadata()

        log.log(21, 'Done!')
        log.log(21, '%s reads are remained', self.N)
    
    def maskFilter(self, mask):
        """
        Use numpy's internal mask mechanism instead.

        Parameters
        ----------
        mask : array of bools
            Indexes of reads to keep
            
        """
        # Uses 16 bytes per read
        length = 0
        ms = mask.sum()
        
        assert mask.dtype == np.bool
        
        self.N = ms
        self.DSnum = self.N
        
        if hasattr(self, "ufragments"):
            del self.ufragmentlen, self.ufragments
            
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
        
        self.rebuildFragments()
