# Created on Thu Sep 13 21:44:24 2018
# Author: XiaoTao Wang

import os, time
import numpy as np
import matplotlib, pickle, glob
matplotlib.use('Agg')
import matplotlib.pyplot as plt

## Plot Settings
import matplotlib.pyplot as plt
# Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['axes.labelsize'] = 15
matplotlib.rcParams['xtick.labelsize'] = 15
matplotlib.rcParams['ytick.labelsize'] = 15
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['ytick.minor.size'] = 5
matplotlib.rcParams['xtick.major.pad'] = 6
matplotlib.rcParams['ytick.major.pad'] = 6
matplotlib.rcParams['xtick.minor.pad'] = 6
matplotlib.rcParams['ytick.minor.pad'] = 6

colorPool = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#348ABD', '#A60628']

def properU(pos):
    
    m_part = int(pos) // 1000000 # Megabase
    k_part = (int(pos) % 1000000) // 1000 # Kilobase
    b_part = int(pos) % 1000 # base pair
    
    if (m_part > 0) and (k_part > 0) and (b_part > 0):
        return ''.join([str(m_part), 'M', str(k_part), 'K', str(b_part), 'bp'])
    elif (m_part > 0) and (k_part > 0) and (b_part == 0):
        return ''.join([str(m_part), 'M', str(k_part), 'K'])
    elif (m_part > 0) and (k_part == 0) and (b_part > 0):
        return ''.join([str(m_part), 'M', str(b_part), 'bp'])
    elif (m_part > 0) and (k_part == 0) and (b_part == 0):
        return ''.join([str(m_part), 'M'])
    elif (m_part == 0) and (k_part > 0) and (b_part > 0):
        return ''.join([str(k_part), 'K', str(b_part), 'bp'])
    elif (m_part == 0) and (k_part > 0) and (b_part == 0):
        return ''.join([str(k_part), 'K'])
    elif (m_part == 0) and (k_part == 0) and (b_part > 0):
        return ''.join([str(b_part), 'bp'])

def printStats(stats, saveTo):

    keys = [k for k in stats if (not k in ['libsize','danglingStart']) and (not k.startswith('dist_freq'))]

    longrange = stats['412_IntraLongRangeReads(>=20Kb)']
    contacts = stats['400_TotalContacts']
    
    Total = stats['000_SequencedReads']
    DUreads = stats['010_DoubleSideMappedReads']
    Uratio = float(DUreads) / Total
    usage = float(contacts) / Total
    longRatio = float(longrange) / contacts
    if '120_SameFragmentReads' in stats:
        selfLig = stats['122_SelfLigationReads']
        dangling = stats['124_DanglingReads']
        Fratio = float(selfLig) / Total
        Dratio = float(dangling) / Total

    # Tree-like print
    with open(saveTo, 'w') as myfile:
        for i in sorted(keys):
            if (i[2] != '0'):
                myfile.write('\t\t')
            elif (i[1] != '0') and (i[2] == '0'):
                myfile.write('\t')
            myfile.write(str(i))
            myfile.write(':  ')
            myfile.write(str(stats[i]))
            myfile.write('\n')
        myfile.write('\nCritical Indicators:\n')
        myfile.write('Double Unique Mapped Ratio = {0} / {1} = {2:.4f}\n'.format(DUreads, Total, Uratio))
        if '120_SameFragmentReads' in stats:
            myfile.write('Self-Ligation Ratio = {0} / {1} = {2:.4f}\n'.format(selfLig, Total, Fratio))
            myfile.write('Dangling-Reads Ratio = {0} / {1} = {2:.4f}\n'.format(dangling, Total, Dratio))
        myfile.write('Long-Range Ratio = {0} / {1} = {2:.4f}\n'.format(longrange, contacts, longRatio))
        myfile.write('Data Usage = {0} / {1} = {2:.4f}\n'.format(contacts, Total, usage))


def typePlot(stats, outfile, dpi = 300):

    from collections import defaultdict

    LT = {}; RT = {}; IT = {}; OT = {}
    Total = defaultdict(int)
    for k in stats:
        if not k.startswith('dist_freq'):
            continue
        _, dist, t = k.split('/')
        if not '-' in dist:
            continue
        dist = int(dist.split('-')[1])
        if (dist > 10000000) and (dist <= 1):
            continue
        if t == '++':
            RT[dist] = stats[k]
        elif t == '--':
            LT[dist] = stats[k]
        elif t == '+-':
            IT[dist] = stats[k]
        else:
            OT[dist] = stats[k]

        Total[dist] += stats[k]
    
    lt = []; rt = []; it = []; ot = []
    for k in sorted(LT):
        if Total[k] > 0:
            lt.append(LT[k]/Total[k])
            rt.append(RT[k]/Total[k])
            it.append(IT[k]/Total[k])
            ot.append(OT[k]/Total[k])
        else:
            lt.append(np.nan)
            rt.append(np.nan)
            it.append(np.nan)
            ot.append(np.nan)
    lt = np.r_[lt]; rt = np.r_[rt]; it = np.r_[it]; ot = np.r_[ot]

    x = np.arange(0.25, 8, 0.125)
    xticks = list(range(1, 8))
    xticklabels = ['10bp', '100bp', '1K', '10K', '100K', '1M', '10M']

    fig = plt.figure(figsize = (10, 7))
    ax = fig.add_subplot(111)
    lines = []
    labels = ['Left Type', 'Right Type', 'Inner Type', 'Outer Type']
    idx = 0
    for y in [lt, rt, it, ot]:
        L = ax.plot(x, y[:x.size], color = colorPool[idx], linewidth = 2)
        lines.extend(L)
        idx += 1
    
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('Genomic Separation', fontsize=18)
    ax.set_ylabel('Type Ratio', fontsize=18)
    
    ax.set_ylim((0, 1))
    ax.set_title('Read Pair Type Statistics', fontsize=20)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
        
    ax.legend(lines, labels, frameon = False, fontsize = 11, labelspacing = 1,
            handletextpad = 1, borderpad = 1, markerscale = 1, numpoints = 1,
            ncol = 2, loc = 'upper right')
        
    plt.savefig(outfile, dpi = dpi, bbox_inches='tight')
    plt.close()
    
def plot_libsize(stats, outplot, dpi = 300):

    libsize = stats['libsize']

    low = np.percentile(libsize, 0.2)
    high = np.percentile(libsize, 99.8)
    libsize = libsize[(libsize<=high) & (libsize>=low)]
    fig = plt.figure(figsize = (10, 7))
    ax = fig.add_subplot(111)
    ax.hist(libsize, bins = 15, color = colorPool[-1])
    ax.set_title('Estimated Library Size Distribution (bp)', fontsize=20)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.savefig(outplot, dpi = dpi, bbox_inches='tight')
    plt.close()

def plot_piechart(stats, outplot, dpi = 300):

    total = stats['000_SequencedReads']
    double = stats['010_DoubleSideMappedReads']
    unmapped = total - double
    pcr = stats['130_DuplicateRemoved']
    intra = stats['410_IntraChromosomalReads']
    inter = stats['420_InterChromosomalReads']
    if '120_SameFragmentReads' in stats:
        samefrag = stats['120_SameFragmentReads']
        counts = [unmapped, pcr, samefrag, intra, inter][::-1]
        colors = ['#6391cf', '#cc4f57', '#a1ca54', '#8f6bb5', '#66bfd1'][::-1]
        labels = ['Unmapped/Single mapped', 'PCR duplicates', 'Self-ligation/Dangling reads',
                  'Intra contacts', 'Inter contacts'][::-1]
        ratios = [c/sum(counts) for c in counts]
    else:
        counts = [unmapped, pcr, intra, inter][::-1]
        colors = ['#6391cf', '#cc4f57', '#8f6bb5', '#66bfd1'][::-1]
        labels = ['Unmapped/Single mapped', 'PCR duplicates', 'Intra contacts', 'Inter contacts'][::-1]
        ratios = [c/sum(counts) for c in counts]

    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    patches, texts, autotexts = ax.pie(counts, colors = colors, autopct='%1.0f%%',
                                   startangle=90, labeldistance=1.03, pctdistance=0.6)
    for w in patches:
        w.set_linewidth(1.5)
        w.set_edgecolor('w')
    
    for i, t in enumerate(autotexts):
        if ratios[i] > 0.05:
            t.set_fontsize(14)
        else:
            t.set_text('')
    
    plt.axis('equal')
    plt.tight_layout()
    ax.legend(patches[::-1], labels[::-1], frameon=False, bbox_to_anchor=(1,0.2,0.08,0.6), ncol=1, fontsize=13, loc=3)

    plt.savefig(outplot, bbox_inches='tight', dpi = dpi)
    plt.close()


def plot_dangling_details(stats, outplot, dpi = 300):

    danglingStart = stats['danglingStart']

    danglingStart = danglingStart[danglingStart<=0.5]
    fig = plt.figure(figsize = (15, 9))
    ax = fig.add_subplot(111)
    ax.hist(danglingStart, bins = 15, color = colorPool[-2])
    ax.set_title('Start Site of Dandling Ends Relative to Fragments (Ratio)')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.savefig(outplot, dpi = dpi)
    plt.close()


def outStatsCache(stats_pool, outpre):

    lockfil = outpre + '.lock'
    while os.path.exists(lockfil):
        time.sleep(3)
    
    lock = open(lockfil, 'wb') # aquire the file lock
    lock.close()

    pattern = outpre + '*'
    allfile = glob.glob(pattern)
    allfile.remove(lockfil)
    if not len(allfile):
        cache = outpre + '.1'
    else:
        counts = [int(f.split('.')[-1]) for f in allfile]
        counts.sort()
        suffix = str(counts[-1]+1)
        cache = '.'.join([outpre, suffix])
    
    with open(cache, 'wb') as out:
        pickle.dump(stats_pool, out)
    
    os.remove(lockfil) # release the lock


def loadStats(cache_pre):

    pattern = cache_pre + '*'
    allfile = glob.glob(pattern)
    if not len(allfile):
        raise Exception('Stats cache {0} cannot be found, exit'.format(cache_pre))
    stats_pool = {}
    for f in allfile:
        with open(f, 'rb') as source:
            tmp = pickle.load(source)
        
        for key in tmp:
            if not key in stats_pool:
                stats_pool[key] = tmp[key]
            else:
                if len(tmp[key]) > len(stats_pool[key]):
                    stats_pool[key] = tmp[key]
    
    return stats_pool

def checkKeys(stats_pool, keys):

    return all([(k in stats_pool) for k in keys])

def update_stats_pool(stats_pool, keys, cache_pre):

    pattern = cache_pre + '*'
    allfile = glob.glob(pattern)
    for f in allfile:
        with open(f, 'rb') as source:
            tmp = pickle.load(source)
        for k in keys:
            if k in tmp:
                stats_pool[k] = tmp[k]
                if os.path.exists(f):
                    os.remove(f) # clean duplicate records

    return stats_pool
