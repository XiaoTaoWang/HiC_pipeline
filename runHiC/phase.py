import sys, pairtools
from bisect import bisect_left
import numpy as np

if pairtools.__version__.startswith('0'):
    from pairtools import _fileio, _headerops, _pairsam_format
else:
    from pairtools.lib import fileio as _fileio
    from pairtools.lib import headerops as _headerops
    from pairtools.lib import pairsam_format as _pairsam_format

def parse_phased_SNPs(infil):

    # four columns: chrom, 0-based pos, hap1, hap2
    hap1 = {}
    hap2 = {}
    SNPs = {}
    with open(infil, 'r') as source:
        for line in source:
            chrom, pos, b1, b2 = line.rstrip().split()
            pos = int(pos)
            if not chrom in SNPs:
                SNPs[chrom] = []
            SNPs[chrom].append(pos)
            hap1[(chrom, pos)] = b1.upper()
            hap2[(chrom, pos)] = b2.upper()
    
    for chrom in SNPs:
        SNPs[chrom] = sorted(SNPs[chrom])
    
    return SNPs, hap1, hap2

def locate_bisect(sorted_list, start, end):

    si = bisect_left(sorted_list, start)
    ei = bisect_left(sorted_list, end)

    return sorted_list[si:ei]

def parse_CIGAR(cigar_str):

    cigar_tuples = []
    if not len(cigar_str):
        return cigar_tuples
    
    if not cigar_str[0].isdigit():
        return cigar_tuples
    
    length = ''
    operation = ''
    for s in cigar_str:
        if s.isdigit():
            length += s
        else:
            operation = s
            cigar_tuples.append((int(length), operation))
            length = ''
            operation = ''
    
    return cigar_tuples

def calculate_positions(cigar_str, ref_pos):

    ref_intervals = []
    read_intervals = []
    read_pos = 0
    cigar = parse_CIGAR(cigar_str)
    for length, operation in cigar:
        if operation in ['M', '=', 'X']:
            ref_intervals.append((ref_pos, ref_pos + length))
            read_intervals.append((read_pos, read_pos + length))
            ref_pos += length
            read_pos += length
        elif operation == 'I':
            read_pos += length
        elif operation in ['D', 'N']:
            ref_pos += length
        elif operation == 'S':
            read_pos += length
        elif operation == 'H':
            continue
        elif operation == 'P':
            continue
    
    return ref_intervals, read_intervals

def parse_sam(record):

    SAM_SEP = _pairsam_format.SAM_SEP
    INTER_SAM_SEP = _pairsam_format.INTER_SAM_SEP
    
    alignments = record.split(INTER_SAM_SEP)
    collect = []
    for align in alignments:
        tmp = align.split(SAM_SEP)
        D = {
            'chrom': tmp[2],
            'pos': int(tmp[3]) - 1,
            'cigar': tmp[5],
            'seq': tmp[9]
        }
        collect.append(D)
    
    return collect

def find_match(align, pair, maximum_dist):

    # maximum_dist --> max_molecule_size of runHiC
    dists = []
    for c, p in pair:
        if c != align['chrom']:
            dists.append(np.inf)
        else:
            dists.append(abs(p - align['pos']))
    
    if min(dists) > maximum_dist:
        return -1
    else:
        return np.argmin(dists)

def phase_pairs(line, SNPs, hap1, hap2, include_readid=False, include_sam=False,
                maximum_dist=750):

    SEP = _pairsam_format.PAIRSAM_SEP
    parse = line.rstrip().split(SEP)
    readID, c1, p1, c2, p2, strand1, strand2, pair_type = parse[:8]
    p1, p2 = int(p1), int(p2)
    sam1 = parse[8]
    sam2 = parse[9]
    mapq1 = parse[10]
    mapq2 = parse[11]
    if not include_readid:
        readID = '.'
    
    collect = [] # store all alignments of the read pair
    collect.extend(parse_sam(sam1))
    collect.extend(parse_sam(sam2))

    # the values are either 0 or 1
    # each value represents one-base-pair support
    # from either haplotype 1 or haplotype 2
    phase_support_1 = []
    phase_support_2 = []
    for align in collect:
        if not align['chrom'] in SNPs:
            continue

        mi = find_match(align, [(c1, p1), (c2, p2)], maximum_dist)
        if mi == -1:
            continue

        ref_intervals, read_intervals = calculate_positions(align['cigar'], align['pos'])
        for ref_iv, read_iv in zip(ref_intervals, read_intervals):
            ref_start, ref_end = ref_iv
            read_start, read_end = read_iv
            snp_hits = locate_bisect(SNPs[align['chrom']], ref_start, ref_end)
            if not len(snp_hits):
                continue

            for p in snp_hits:
                hap1_ = hap1[(align['chrom'], p)]
                hap2_ = hap2[(align['chrom'], p)]
                rp = read_start + (p - ref_start)
                obs = align['seq'][rp].upper()
                if obs == hap1_:
                    if mi == 0:
                        phase_support_1.append(0)
                    else:
                        phase_support_2.append(0)
                elif obs == hap2_:
                    if mi == 0:
                        phase_support_1.append(1)
                    else:
                        phase_support_2.append(1)
    
    if not len(phase_support_1):
        phase1 = '.'
    else:
        count_0 = phase_support_1.count(0)
        count_1 = phase_support_1.count(1)
        if count_0 == count_1:
            phase1 = '.'
        elif count_0 > count_1:
            phase1 = '0'
        else:
            phase1 = '1'

    if not len(phase_support_2):
        phase2 = '.'
    else:
        count_0 = phase_support_2.count(0)
        count_1 = phase_support_2.count(1)
        if count_0 == count_1:
            phase2 = '.'
        elif count_0 > count_1:
            phase2 = '0'
        else:
            phase2 = '1'
    
    if not include_sam:
        cols = [readID, c1, str(p1), c2, str(p2), strand1, strand2,
                pair_type, phase1, phase2, mapq1, mapq2]
    else:
        cols = [readID, c1, str(p1), c2, str(p2), strand1, strand2,
                pair_type, sam1, sam2, mapq1, mapq2, phase1, phase2]
    
    return cols


def phase_pipeline(pairs_path, output, snp_fil, nproc_in=3, nproc_out=8,
                   include_readid=False, include_sam=False, maximum_dist=750):

    SNPs, hap1, hap2 = parse_phased_SNPs(snp_fil)
    instream = _fileio.auto_open(pairs_path, mode='r', nproc=nproc_in)
    outstream = _fileio.auto_open(output, mode='w', nproc=nproc_out)

    # Parse the input stream:
    header, body_stream = _headerops.get_header(instream)

    # modify the header
    columns = header.pop().split()
    if not include_sam:
        sam1_col = columns.index('sam1')
        sam2_col = columns.index('sam2')
        columns[sam1_col] = 'phase1'
        columns[sam2_col] = 'phase2'
    else:
        columns.extend(['phase1', 'phase2'])
    header.append(' '.join(columns))

    # write out the updated header
    outstream.writelines((l+'\n' for l in header))
    outstream.flush()

    # iterate each line of the input pairsam file
    for line in body_stream:
        cols = phase_pairs(line, SNPs, hap1, hap2, include_readid=include_readid,
                           include_sam=include_sam, maximum_dist=maximum_dist)
        outstream.write('\t'.join(cols) + '\n')

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()