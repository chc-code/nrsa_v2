#! /usr/bin/env python3

"""
functions used for NRSA package
"""
import os
import sys
import re
import time
import json, pickle
import gzip
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
import fisher
import bisect
pw_code = os.path.dirname(os.path.realpath(__file__))


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr
def red(s, p=False):
    s_new = f'\u001b[31;1m{s}\u001b[0m'
    if p:
        print(s_new)
    return s_new
def green(s, p=False):
    s_new = f'\u001b[32;1m{s}\u001b[0m'
    if p:
        print(s_new)
    return s_new

def getlogger(fn_log=None, logger_name=None, nocolor=False, verbose=False):
    import logging
    logger_name = logger_name or "main"
    
    try:
        logger = logging.getLogger(logger_name)
    except:
        logger = logging.getLogger('terminal')

    class CustomFormatter(logging.Formatter):
    
        def __init__(self, nocolor=False):
            self.nocolor = nocolor
        colors = {
            'black': '\u001b[30;20m',
            'red': '\u001b[31;20m',
            'r': '\u001b[31;20m',
            'bold_red': '\u001b[31;1m',
            'rb': '\u001b[31;1m',
            'green': '\u001b[32;20m',
            'g': '\u001b[32;20m',
            'gb': '\u001b[32;1m',
            'yellow': '\u001b[33;20m',
            'blue': '\u001b[34;20m',
            'b': '\u001b[34;20m',
            'purple': '\u001b[35;1m',
            'p': '\u001b[35;1m',
            'grey': '\u001b[38;20m',
        }
        FORMATS = {
            logging.WARNING: colors['purple'],
            logging.ERROR: colors['bold_red'],
            logging.CRITICAL: colors['bold_red'],
        }
    
        def format(self, record):
            format_str = "%(asctime)s  %(levelname)-6s %(funcName)-20s  line: %(lineno)-5s  %(message)s"
            reset = "\u001b[0m"
            log_fmt = None
            
            record.msg = str(record.msg)
            if self.nocolor:
                pass
            elif '@' in record.msg[:10]:
                try:
                    icolor, tmp = record.msg.split('@', 1)
                    log_fmt = self.colors.get(icolor)
                    if log_fmt:
                        record.msg = tmp
                except:
                    raise
                    pass
            else:
                log_fmt = self.FORMATS.get(record.levelno)
            if log_fmt:
                record.msg = log_fmt + record.msg + reset
            formatter = logging.Formatter(format_str, datefmt='%Y-%m-%d %H:%M:%S')
            return formatter.format(record)
    
    logger.setLevel('DEBUG')
    handler_names = {_.name for _ in logger.handlers}
    if 'console' not in handler_names:
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(CustomFormatter(nocolor=nocolor))
        console.setLevel('DEBUG' if verbose else 'INFO')
        console.name = 'console'
        logger.addHandler(console)

    if fn_log and 'file' not in handler_names:
        fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
        fh_file.setLevel('DEBUG')
        fh_file.setFormatter(CustomFormatter())
        fh_file.name = 'file'
        logger.addHandler(fh_file)
    return logger
logger = getlogger('NRSA.run.log', 'NRSA')

bases_set = set('ATGCatgc')
def refine_chr(chr_):
    chr_ = chr_.lower()
    for _ in ['chromosome', 'chrom', 'chr']:
        chr_ = chr_.replace(_, '')
    return {'23': 'x', '24': 'y', '25': 'mt', 'm': 'mt'}.get(chr_) or chr_

def build_idx_for_fa(fn_fa):
    # >chr1
    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    res = {'note': 'the line_width does not include the end newline symbol, it is the count of the regular char only'}
    line_width = {}
    # line_width,  the most frequent sequence in each line. the width does not include the end newline symbol
    # other keys = each chr, removed the chr prefix,  the value is like [s_byte, length]
    if '.gz' in fn_fa:
        logger.error(f'Please use decompressed fasta file')
        sys.exit(1)
        return None

    lb = os.path.basename(fn_fa).rsplit('.', 1)[0]
    fn_fa = os.path.abspath(fn_fa)
    
    # check if have write permission to the folder of the fa file
    if not os.access(os.path.dirname(fn_fa), os.W_OK):
        home = os.path.expanduser("~")
        fa_out_dir = f'{home}/.nrsa'
        os.makedirs(fa_out_dir, exist_ok=True)
    else:
        fa_out_dir = os.path.dirname(fn_fa)
        
    fn_idx = f'{fa_out_dir}/{lb}.idx.pkl'
    if os.path.exists(fn_idx):
        logger.debug(f'loading fasta index file: {fn_idx}')
        with open(fn_idx, 'rb') as f:
            res = pickle.load(f)
            return res

    logger.info(f'building index for {fn_fa}')
    chr_ = None
    with open(fn_fa) as f:
        while True:
            i = f.readline()
            if not i:
                break
            if i[0] == '>':
                chr_ = refine_chr(i[1:-1].lower())

                s_pos = f.tell()
                res[chr_] = [s_pos, 0]
            elif chr_:
                len1 = len(i) - 1
                res[chr_][1] += len1
                line_width.setdefault(len1, 0)
                line_width[len1] += 1
    
    tmp = sorted(line_width.items(), key=lambda _: _[1], reverse=True)
    word_wrap, freq = tmp[0]
    logger.debug(f'word wrap = {word_wrap}, lines width = {tmp}')
    res['line_width'] = word_wrap
    with open(fn_idx, 'wb') as o:
        pickle.dump(res, o)
    return res

def get_seqence_from_fa(idx, fh, chr_, start, end):
    """
    idx is from loading of the json {fn_fa}.idx.json file, fh is the opened file handle of the fasta file
    already converted to uppercase
    """
    chr_ = refine_chr(chr_)
    line_width = idx['line_width']
    len1 = end - start + 1
    if chr_ not in idx:
        logger.error(f'invalid chromosome {chr_}, valid = {sorted(idx.keys())}')
        return 1
    byte_0, chr_total_len = idx[chr_]
    if end > chr_total_len:
        logger.error(f'end position out of range! chr{chr_} : total length = {chr_total_len} ')
        return 1
    lines_start = start // line_width
    lines_end = end // line_width
    lines_span = lines_end - lines_start
    start_byte = start + byte_0 -1 + lines_start
    fh.seek(start_byte)
    return fh.read(len1 + lines_span).replace('\n', '').upper()

def build_idx_for_bed(fn_bed, fn_idx, step=100):
    """
    built the index for the bed file, the index file is in json format
    key1 = chr, key2 = position (by default, step size = 1000), value is the byte position in the bed file
    the step size is used to reduce the number of lines to be read, so that the index file is smaller
    for 2G bed files, the index file size is around 9M for chunksize= 1000, and around 50M for chunksize of 100
    the overall run time will be reduced by 70% for chunksize = 100 compared to chunksize = 1000
    stratify by strand, for plus strand, the end position is used to get the chunk number, however, because the end position is not sorted, so file byte position of a larger chunk number can be smaller than a smaller chunk number.
    This will lead to problem, e.g. we searched for chunk100, and start from the file byte postion, and we expect we can find all the variants in chunk100 and chunk101 and after, but in fact, some reads in chunk101 or chunk102 can be before this file position.
    The solution is, after building this index, we'll do adjustment to the plus strand index, and for a specific chunk, the byte position is set to be the minimum of the byte position of the chunk and the  chunks after
    """
    lb = get_lb(fn_bed)
    res = {'step': step, 'fn_lb': lb}
    chunks = 0
        
    with open(fn_bed) as f:
        while True:
            i = f.readline()
            if not i:
                break
            line = i[:-1].split('\t')
            try:
                # 3 = id, 4 = score, 5 = strand
                chr_, start, end, strand = line[0], line[1], line[2], line[5]
                # chr_, start, strand = line[0], line[1], line[5]
            except:
                logger.error(f'{fn_bed}: invalid line format: {line}, position = {f.tell()}')
                sys.exit(1)
            ires = res.setdefault(chr_, {}).setdefault(strand, {})
            read_end = int(start) + 1 if strand == '-' else int(end)
            chunk = read_end // step   # the chunk number
            # chunk = (int(start) + 1) // step
            if chunk not in ires:
                pos_line_end = f.tell()
                # the above pos is the start of the next line, so need to get the position of the current line
                pos_line_start = pos_line_end - len(i)
                chunks += 1
                ires[chunk] = pos_line_start
    res['chunks'] = chunks
    # update the chromosome name
    res = {refine_chr(k): v for k, v in res.items()}

    with open(fn_idx, 'wb') as o:
        pickle.dump(res, o)
    logger.info(f'Buildingn index done for {os.path.basename(fn_bed)}, chunks = {chunks}')
    return res

def check_is_sorted(fn_bed):
    """
    check if the bed file is sorted or not
    """
    chr_prev = None
    chr_set = set()
    with open(fn_bed) as f:
        last = 0
        for i in f:
            chr_, start = i.split('\t', 2)[:2]
            start = int(start)
            if chr_ != chr_prev:
                if chr_ in chr_set:
                    # the chr_ has been seen before
                    return False
                chr_set.add(chr_)
                chr_prev = chr_
            elif start < last:
                return False
            last = start
    return True


def get_ref(organism, fa_in=None, gtf=None):
    # hg19.chrom.sizes
    # hg19_rCRS.fa
    # hg19-tss-tts.txt
    # RefSeq-hg19-exon-formatted.gtf
    # RefSeq-hg19-tss-tts.txt
    # RefSeq-hg19-tss.txt
    # UCSC-hg19-exon-uniq-reformated-filter.gtf
    # UCSC-hg19-filter-tss.txt
    pw_code = os.path.dirname(os.path.realpath(__file__))
    pw_ref = f'{pw_code}/ref'
    
    ref = {
        'fa': f'{pw_code}/fa/{organism}/{organism}.fa',

        'gtf': f'{pw_ref}/{organism}/RefSeq-{organism}.gtf', # use the new gtf file
        'tss': f'{pw_ref}/{organism}/RefSeq-{organism}-tss.txt',
        'tss_tts': f'{pw_ref}/{organism}/RefSeq-{organism}-tss-tts.txt',
        # 'tss_filter': f'{pw_ref}/{organism}/UCSC-{organism}-filter-tss.txt',
    }
    
    if fa_in:
        ref['fa'] = fa_in
    if gtf:
        ref['gtf'] = gtf
        
    err = 0
    not_found = []
    for k in list(ref):
        fn = ref[k]
        if not os.path.exists(fn):
            if os.path.exists(f'{fn}.gz'):
                ref[k] = f'{fn}.gz'
            else:
                not_found.append(fn)
                err = 1
    if err:
        logger.error(f"the following files are not found:")
        print('\n\t'.join(not_found))
        return None

    return ref

def bench(s, lb, d):
    now = time.time()
    d[lb] += now - s
    return now

def get_FDR_per_sample(fn_peak, fn_fdr, header_fdr):
    """
    get the pvalue and FDR based on the peak file, retain only 1 transcript per gene, (keep the transcript with the highest pp count)
    """
    df_peak = pd.read_csv(fn_peak, sep='\t')
    # collapse genes
    # Transcript	Gene	ppc	ppm	ppd	pps	gbc	gbm	gbd	pauseIndex
    # NM_000014.5	A2M	40	50	0.8	12-:0	43	47522	0.0009	884.13023
    df_peak = df_peak.groupby('Gene').apply(lambda x: x.loc[x['ppc'].idxmax()], )
    


def get_peak(count_file_info, chr_, strand, start, end, s_gene, gene_seq, window_size, step_size, only_mapped=False, time_cost=None, already_parsed_windows=None, get_summit_flag=False):
    """
    use the new method, i.e. in the preprocessing step, get the mapped reads count for each 10bp window, then directly sum the mapped reads count for the new window
    """
    dummy = [0,0,0,100]  # for test purpose, return this dummy value
    
    # fn_lb, fn_bed_count, idx_d, res_chunk
    # count_file_info = {'fn': fn_bed_count, 'idx_d': idx_d, 'res_chunk': res_chunk, 'fh':  file handle}
    peak_window = [None, None] # [mapped_reads_count, window_start_pos]
    count_file_idx = count_file_info['idx_d']
    count_file_chunk = count_file_info['res_chunk']
    count_file_fh = count_file_info['fh']
    
    count_chunk_size = count_file_chunk['bin_size'] # default = 10bp
    idx_chunk_size = count_file_idx['bin_size'] # default = 50bp
    
    k_whole = f'{chr_}@{strand}@{start}@{end}'
    if already_parsed_windows and k_whole in already_parsed_windows['peak']:
        time_cost['reused']['whole_region'] += 1
        return already_parsed_windows['peak'][k_whole]
    
    for pos in range(start, end - window_size + 2, step_size):
        s_window, e_window = pos, pos + window_size - 1
        # k = f'{chr_}@{strand}@{s_window}@{e_window}'
        # continue
        # if already_parsed_windows and k in already_parsed_windows:
        #     ct_mapped_read = already_parsed_windows[k]
        #     time_cost['reused']['window'] += 1
        # else:
            # continue
            # each window, split into 3 parts
            # part1 = s_window to the first bin boundary
            # part2 = the bins in the window
            # part3 = the last bin boundary to e_window
        
        bin_start = s_window // count_chunk_size
        bin_end = e_window // count_chunk_size
        
        window_both_side = []
        p2_bin_list = []
        ct_mapped_read = 0
        
        p1_s, p1_e = (s_window, (bin_start + 1) * count_chunk_size - 1)
        if bin_end == bin_start:
            window_both_side.append([p1_s, e_window])
            # there is no p2 and p3
        else:
            p3_s, p3_e = ((bin_end) * count_chunk_size, e_window)
            window_both_side.append([p1_s, p1_e])
            window_both_side.append([p3_s, p3_e])
            if bin_end > bin_start + 1:
                p2_bin_list = range(bin_start + 1, bin_end)
        for ps, pe in window_both_side:
            k1 = f'{chr_}@{strand}@{ps}@{pe}'
            if already_parsed_windows and k1 in already_parsed_windows:
                ct_mapped_read += already_parsed_windows[k1]
                time_cost['reused']['side'] += 1
            else:
                ict_side = get_mapped_reads(count_file_fh, count_file_idx, chr_, strand, ps, pe)
                ct_mapped_read += ict_side
                already_parsed_windows[k1] = ict_side
        if p2_bin_list:
            for p2_bin in p2_bin_list:
                ct_mapped_read += count_file_chunk[chr_][strand].get(p2_bin, 0)
            
        if peak_window[0] is None or ct_mapped_read > peak_window[0]:
                peak_window = [ct_mapped_read, pos]

    if not peak_window[0]:
        ires = [0, 0, 0]
        summit_pos = start if strand == '+' else end
        summit_pos_str = f'{chr_}{strand}:{summit_pos}'
        ires.append(summit_pos_str)
        already_parsed_windows['peak'][k_whole] = ires
        return ires
    s_peak_window = peak_window[1]
    e_peak_window = s_peak_window + window_size - 1
    
    if get_summit_flag and peak_window[0] > 0:
        summit_pos = get_summit(count_file_fh, count_file_idx, chr_, strand, s_peak_window, e_peak_window)
    else:
        summit_pos = s_peak_window
    summit_pos_str = f'{chr_}{strand}:{summit_pos}'
    
    

    # use the precalculated invalid bases set, but it is actually slower
    # if invalid_bases_set:
    #     invalid_bases_in_window = len(set(range(s_peak_window, e_peak_window + 1)) & invalid_bases_set)
    # else:
    #     invalid_bases_in_window = 0
    # mappable_sites = window_size - invalid_bases_in_window
    if gene_seq:
        window_seq = gene_seq[s_peak_window - s_gene: e_peak_window - s_gene + 1]
        mappable_sites = window_seq.count('A') + window_seq.count('T') + window_seq.count('G') + window_seq.count('C')
    else:
        mappable_sites = window_size

    ratio = peak_window[0] / mappable_sites if mappable_sites else 0
    ires = [peak_window[0], mappable_sites, ratio, summit_pos_str] # [mapped_reads_count, mappable_sites, ratio, summit_pos]
    already_parsed_windows['peak'][k_whole] = ires

    return ires

def get_mapped_reads(count_file_fh, count_file_idx, chr_, strand, start, end):
    # logger.warning('mofify here')
    return 20

def get_mapped_reads_real(count_file_fh, count_file_idx, chr_, strand, start, end):
    """
    get the coverage of the region, the coverage is the number of **read ends** at each position in this region
    each call will cost around 100us.
    """
    chunk_size = count_file_idx['bin_size']
    start, end = sorted([start, end])
    chunk_s = start // chunk_size
    chunk_e = end // chunk_size
    strand_idx = 2 if strand == '+' else 3 # plus strand is 3rd column, minus strand is 4th column
    idx_bed_chr = count_file_idx[chr_]
    fh_byte_start = None
    
    if chunk_e > chunk_s:
        all_possible_chunks = sorted(range(chunk_s, chunk_e + 1))
        for ichunk in all_possible_chunks:
            if ichunk in idx_bed_chr:
                byte_pos_tmp = idx_bed_chr[ichunk]
                if fh_byte_start is None or byte_pos_tmp < fh_byte_start:
                    fh_byte_start = byte_pos_tmp
    else:
        if chunk_s in idx_bed_chr:
            fh_byte_start = idx_bed_chr[chunk_s]
    if fh_byte_start is None:
        # no reads in this region
        return 0

    count_file_fh.seek(fh_byte_start)
    n_parsed = 0
    for i in count_file_fh:
        # 1	13033	0	2
        line = i[:-1].split('\t')
        chr_line, pos, ict = line[0], line[1], line[strand_idx]
        pos = int(pos)
        ict = int(ict)
        if start <= pos <= end:
            n_parsed += ict
            continue
        if  pos > end or chr_line != chr_:
            # there will be no more lines matching the region
            break
    return n_parsed

def get_summit(count_file_fh, count_file_idx, chr_, strand, start, end):
    """
    get the positionn with the hightest count in the region
    """
    chunk_size = count_file_idx['bin_size']
    start, end = sorted([start, end])
    chunk_s = start // chunk_size
    chunk_e = end // chunk_size
    strand_idx = 2 if strand == '+' else 3 # plus strand is 3rd column, minus strand is 4th column
    idx_bed_chr = count_file_idx[chr_]

    all_possible_chunks = sorted(range(chunk_s, chunk_e + 1))
    fh_byte_start = None

    if chunk_e > chunk_s:
        all_possible_chunks = sorted(range(chunk_s, chunk_e + 1))
        for ichunk in all_possible_chunks:
            if ichunk in idx_bed_chr:
                byte_pos_tmp = idx_bed_chr[ichunk]
                if fh_byte_start is None or byte_pos_tmp < fh_byte_start:
                    fh_byte_start = byte_pos_tmp
    else:
        if chunk_s in idx_bed_chr:
            fh_byte_start = idx_bed_chr[chunk_s]
    if fh_byte_start is None:
        # no reads in this region
        return 0

    count_file_fh.seek(fh_byte_start)
    
    summit = [0, 0]
    for i in count_file_fh:
        # 1	13033	0	2
        line = i[:-1].split('\t')
        chr_line, pos, ict = line[0], line[1], line[strand_idx]
        pos = int(pos)
        ict = int(ict)
        if start <= pos <= end:
            if summit[0] < ict:
                summit = [ict, pos]
            continue

        if  pos > end or chr_line != chr_:
            # there will be no more lines matching the region
            break
    return summit[1]


def draw_box_plot(n_gene_cols, pw_out, out_name, n_rep1, n_rep2=None, gn_list=None):
    """
    gn_list, the genes to be plotted, if None, all genes will be plotted
    pw_out: folder name for pause_PROseq.pl outputs
    n_rep1 and n_rep2:  the replacates number for condition1 and conditon2
    n_rep2 is optional
    """
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    n_rep2 = n_rep2 or 0
    gn_list = set(gn_list) if gn_list else None
    # check the pause_PROseq.pl output files
    for i in ['known_gene', 'intermediate']:
        ipw = f'{pw_out}/{i}'
        if not os.path.exists(ipw):
            logger.error(f'output folder of pause_PROseq.py not found: {i}')
            return 1

    # read pp_gb results
    fn_pp_gb = f'{pw_out}/known_gene/normalized_pp_gb.txt'
    fn_pindex = f'{pw_out}/known_gene/pindex.txt'
    fno = f'{pw_out}/intermediate/boxplot_data.txt'
    
    header_str = []
    out_str = {} # key = gn, v = list combined from pp_gb and pindex
    
    n_total_sam = n_rep1 + n_rep2
    
    # parse pp_gb results
    with open(fn_pp_gb) as f:
        header_in = f.readline().strip().split('\t')
        if header_in[0] != 'Transcript':
            logger.error(f'invalid header in {fn_pp_gb}')
            return 1
        idx_gene_cols = list(range(n_gene_cols))
        idx_list = []
        idx_list += idx_gene_cols
        idx_list += [i for i, col in enumerate(header_in) if col[:4] == 'ppc_']
        idx_list += [i for i, col in enumerate(header_in) if col[:4] == 'ppd_']
        header_str = [header_in[i] for i in idx_list]

        for i in f:
            line = i.strip().split('\t')
            gn = line[0]
            if gn_list is None or gn in gn_list:
                k = f'{gn}@{line[1]}'
                ires = [line[i] for i in idx_list]
                out_str[k] = ires

    # parse pindex results
    with open(fn_pindex) as f:
        header_in = f.readline().strip().split('\t')
        # transcript, gene, sam1-pindex, sam1-pvalue, sam1-FDR, sam2-pindex, sam2-pvalue, sam2-FDR...
        if header_in[0] != 'Transcript':
            logger.error(f'invalid header in {fn_pindex}')
            return 1
        idx_list = [i for i, col in enumerate(header_in) if col.endswith('pindex')]
        header_str += [header_in[i] for i in idx_list]
        for i in f:
            line = i.strip().split('\t')
            gn = line[0]
            if gn_list is None or gn in gn_list:
                k = f'{gn}@{line[1]}'
                if k not in out_str:
                    logger.warning(f'{gn} in pindex not found in pp_gb')
                    continue
                ires = out_str[k]
                ires += [line[i] for i in idx_list]

    # write to file
    with open(fno, 'w') as o:
        o.write('\t'.join(header_str) + '\n')
        for k in sorted(out_str):
            o.write('\t'.join(out_str[k]) + '\n')
    
    # draw the box plot, adapted from boxplot.R
    # R CMD --args outdir=\"$inter_dir\" pname=\"$name\" custom=1 rep1=$rep1 rep2=$rep2
    df_box = pd.read_csv(fno, sep='\t')
    headers = df_box.columns
    sam_list = [name.split('ppc_')[1] for name in headers[n_gene_cols: n_gene_cols + n_total_sam]]
    cols_pp_density = [_ for _ in df_box.columns if _.startswith('ppd_')]
    cols_gb_density = [_ for _ in df_box.columns if _.startswith('ppc_')]
    cols_pindex = [_ for _ in df_box.columns if _.endswith('pindex')]
    
    pp_density = df_box.iloc[:, cols_pp_density]
    gb_density = df_box.iloc[:, cols_gb_density]
    pindex =     df_box.iloc[:, cols_pindex]

    pp_density.columns = sam_list
    gb_density.columns = sam_list
    pindex.columns = sam_list
    
    # plot
    fn_plot_pp = f'{pw_out}/known_gene/{out_name}_pp_density.pdf'
    fn_plot_gb = f'{pw_out}/known_gene/{out_name}_gb_density.pdf'
    fn_plot_pindex = f'{pw_out}/known_gene/{out_name}_pausing_index.pdf'
    
    def plot_task(fn, ylabel, data):
        with PdfPages(fn) as pdf:
            plt.figure(figsize=(2, 4))
            
            # Split the data into two groups based on n_rep1
            data_rep1 = data.iloc[:, :n_rep1]
            data_rep2 = data.iloc[:, n_rep1:]
            
            # Plot the boxplot for the first n_rep1 samples (blue color)
            plt.boxplot(data_rep1 * 1000, positions=range(1, n_rep1 + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='blue'), capprops=dict(color='blue'), whiskerprops=dict(color='blue'), medianprops=dict(color='blue'))
            
            # Plot the boxplot for the rest of the samples (red color)
            if n_rep2:
                plt.boxplot(data_rep2 * 1000, positions=range(n_rep1 + 1, n_total_sam + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='red'), capprops=dict(color='red'), whiskerprops=dict(color='red'), medianprops=dict(color='red'))
            
            plt.ylabel(ylabel)
            plt.xticks(range(1, n_total_sam + 1), sam_list, rotation='vertical')
            plt.tight_layout()
            pdf.savefig()
            
    plot_task(fn_plot_pp, 'Reads per Kb (RPK)', pp_density)
    plot_task(fn_plot_gb, 'Reads per Kb (RPK)', gb_density)
    plot_task(fn_plot_pindex, 'Pausing Index', pindex)

    def plot_hist(n_sam, n_prev_col, ppd, gbd, condition_sn):

        plt.figure(figsize=(12, 6 * (n_sam - 1)))
        for k in range(2, n_sam + 1):
            
            xlabel = f"log2({sam_list[0]}/{sam_list[k]})"
            plt.subplot((n_sam - 1), 2, k - 1)
            plt.hist(np.log2(ppd[:, n_prev_col] / ppd[:, k + n_prev_col]), bins=100)
            plt.title(f"ppd: rep1vs.rep{k}")
            plt.xlabel(xlabel)

            plt.subplot((n_sam - 1), 2, k)
            plt.hist(np.log2(gbd[:, n_prev_col] / gbd[:, k + n_prev_col]), bins=100)
            plt.title(f"gbd: rep1vs.rep{k}")
            plt.xlabel(xlabel)

        plt.savefig(f"{pw_out}/known_gene/Reps-condition{condition_sn}.tif")
        plt.close()

    if n_rep1 > 1:
        plot_hist(n_rep1, 0, pp_density, gb_density, 1)
    if n_rep2 > 1:
        plot_hist(n_rep2, n_rep1, pp_density, gb_density, 2)
    

def draw_heatmap_pindex(n_gene_cols, pw_out):
    """
    plot for the pindex change
    """
    from matplotlib import pyplot as plt
    
    fn_pindex_change = f'{pw_out}/known_gene/pindex_change.txt'
    if not os.path.exists(fn_pindex_change):
        logger.error(f'{fn_pindex_change} not found')
        return 1
    
    
    data_plot = pd.read_csv(fn_pindex_change, sep='\t')
    # Transcript	Gene	log2fc	pvalue	FDR
    idx_log2fc = n_gene_cols
    idx_pvalue = n_gene_cols + 1
    idx_fdr = n_gene_cols + 2

    cols = data_plot.columns
    # filter out the genes with FDR > 0.05
    data_plot = data_plot[data_plot.iloc[:, idx_fdr] < 0.05]
    # sort by log2FC
    col_logfc = cols[idx_log2fc]
    data_plot = data_plot.sort_values(by=cols[idx_log2fc])

    color_min = min(data_plot.iloc[:, idx_log2fc])
    color_max = max(data_plot.iloc[:, idx_log2fc])

    colors = list(mcolors.ColorConverter().to_rgb(color) for color in ["green", "yellow", "red"])
    my_palette = mcolors.LinearSegmentedColormap.from_list("my_palette", colors, N=209)

    fig, ax = plt.subplots(figsize=(2, 6))
    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=my_palette), ax=ax, orientation='horizontal')
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels([color_min, 0, color_max])

    img = ax.imshow(x.iloc[:, 2].values.reshape(-1, 1), cmap=my_palette, aspect='auto')

    plt.savefig(f'{pw_out}/known_gene/pindex_change.pdf', bbox_inches='tight')
    plt.close()


def draw_heatmap_pp_change(n_gene_cols, pw_out, fn_glist, fls_ctrl, fls_case, ref_fls, region_size=5000, bin_size=200, outname='heatmap'):
    """
    draw heatmap for pp change
    fn_glist = gene list for plot (the IDs as listed in the first column of pp_change.txt)
    fls_ctrl, fls_case, foramt is like [fn_bed, idx_bed, fn_lb]
    region_size, upstream and downstream distance relative to TSS for plotting PRO-seq signal (bp, default: 5000),should can be divided by bin size
    bin_size (bp, default: 200),should can be divided by 2
    """
    from matplotlib import pyplot as plt
    fls_case = fls_case or []
    
    ppchange = {}
    fn_pp_change = f'{pw_out}/known_gene/pp_change.txt'
    fn_tss = ref_fls['tss']
    fno = f'{pw_out}/intermediate/genes_for_heatmap.pdf'
    fn_split_bin_bed = f'{pw_out}/intermediate/region_bed.tmp'
    fn_data_bed = f'{pw_out}/intermediate/data_bed.tmp'
    fn_tss_padding = f'{pw_out}/intermediate/infile_bed.tmp' # expand the fn_tss with upstream and downstream region by region_size

    if bin_size % 2 or region_size % bin_size:
        logger.error(f'bin_size should be even and region_size should be divisible by bin_size')
        return 1

    bin_number = region_size * 2 / bin_size + 1  # *2 because of the upstream and downstream

    with open(fn_glist) as f:
        gn_list = {i.strip() for i in f}
    
    # read pp_change results
    # Transcript	Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
    with open(fn_pp_change) as f:
        header = f.readline().strip().split('\t')
        idx_transcript = header.index('Transcript')
        idx_logfc = header.index('log2FoldChange')
        for i in f:
            line = i.strip().split('\t')
            transcript_id = line[idx_transcript]
            if transcript_id in gn_list:
                logfc = 'NA' if line[idx_logfc] == 'NA' else float(line[idx_logfc])
                ppchange[transcript_id] = logfc
    
    # tss file
    # chr1	11874	11874	NR_046018.2	+
    tss_res = []
    tss_with_padding = []
    transcript_id_list = []
    strand_info = {}
    with open(fn_tss) as f:
        for i in f:
            line = i.strip().split('\t')
            chr_, pos, transcript_id, strand = [line[_] for _ in [0, 1, 3, 4]]
            pos_int = int(pos)
            chr_ = refine_chr(chr_)
            if transcript_id in ppchange and ppchange[transcript_id] != 'NA':
                logfc = ppchange[transcript_id]
                row = [transcript_id, chr_, pos, strand, logfc]
                transcript_id_list.append(transcript_id)
                strand_info[transcript_id] = strand

                left_boundary = pos_int - region_size - bin_size/2
                right_boundary = pos_int + region_size + bin_size/2
                row_tmp = [chr_, left_boundary, right_boundary, transcript_id, strand]

                tss_res.append(row)
                tss_with_padding.append(row_tmp)
    if len(tss_res) == 0:
        logger.error(f'no gene in tss file match with the gene list, please check the gene names')
        return 1
    
    # write tss_with_beding bed and data_bed
    with open(fn_tss_padding, 'w') as o:
        for i in tss_with_padding:
            o.write('\t'.join(map(str, i)) + '\n')
    with open(fn_data_bed, 'w') as o:
        for i in tss_res:
            o.write('\t'.join(map(str, i)) + '\n')

    # use bedtools to split the fn_tss_padding, srcwinnum, use the orignal region_label and bin number
    # after running, the last column is like  origlb_1, origlb_2, ....  _1, and _2 are the bin number
    cmd = f'bedtools makewindows -b {fn_tss_padding} -n {bin_number} -i srcwinnum > {fn_split_bin_bed}'
    retcode = run_shell(cmd)
    if retcode:  # fail to makewindow
        logger.error(f'failed to run bedtools makewindows')
        return 1
    
    # load the normalization factors (nf.txt)
    fn_nf = f'{pw_out}/intermediate/nf.txt'
    factors_dict = {}
    with open(fn_nf) as f:
        f.readline() # header
        for i in f:
            sample, factor = i.strip().split('\t')
            factors_dict[sample] = float(factor)
            
    # main part, get the overlap of the input bed files with the split bin bed
    # coverage_by_strand_flag = '-s'  # by_strand, if used, need to update the reference tss file
    coverage_by_strand_flag = '' # current setting
    for condition, fls in {'control': fls_ctrl, 'case': fls_case}.items():
        fn_count_sum = f'{pw_out}/intermediate/{condition}.count'
        count = {}
        n_sam_condition = len(fls)
        for fn_bed, idx_bed, fn_lb in fls:
            norm_factor = factors_dict[fn_lb]
            fn_coverage_tmp = f'{pw_out}/intermediate/{fn_lb}.coverage_count.tmp'
            # do we need to calculate the coverage by strand??? 
            cmd = f'bedtools coverage -a {fn_split_bin_bed} -b {fn_bed} -counts {coverage_by_strand_flag} > {fn_coverage_tmp}'
            retcode = run_shell(cmd)
            if retcode:
                logger.error(f'failed to run bedtools coverage for {fn_lb}')
                return 1
            # process the coverage file
            # chr1	323932	323952	NR_028322.1_8   10
            with open(fn_coverage_tmp) as f:
                for i in f:
                    _, transcript_id_chunk, ict = i.strip().rsplit('\t', 2)
                    transcript_id, bin_sn = transcript_id_chunk.rsplit('_', 1)
                    bin_sn = int(bin_sn)
                    ict = int(ict)
                    ict_norm = ict * norm_factor
                    count.setdefault(transcript_id, {})[bin_sn] = ict_norm


        with open(fn_count_sum, 'w') as o:
            half_bin_number = (bin_number - 1) // 2
            # header = ['Transcript'] + [f'up_{i}' for i in range(half_bin_number, 0, -1)] + ['tss'] + [f'down_{i}' for i in range(1, half_bin_number + 1)]
            header = [f'up_{i}' for i in range(half_bin_number, 0, -1)] + ['tss'] + [f'down_{i}' for i in range(1, half_bin_number + 1)]
            print('\t'.join(header), file=o)
            for transcript_id, v1 in count.items():
                strand = strand_info[transcript_id]
                bin_count_list = []
                bin_number_order = range(1, bin_number + 1) if strand == '+' else range(bin_number, 0, -1)
                
                for bin_sn in bin_number_order:
                    bin_count_list.append(str(v1[bin_sn]))
                tmp = "\t".join(bin_count_list)
                print(f'{transcript_id}\t{tmp}', file=o)
    
    # heatmap.R  --args file=\"$data_bed\" outdir=\"$inter_dir\" pname=\"$outname\" window=$bin_size region=$region_size
    # fn_data_bed => data_bed.tmp
    # each row is [transcript_id, chr_, pos, strand, logfc]
    # read the count table
    df_case = pd.read_csv(f'{pw_out}/intermediate/case.count', sep='\t')
    df_ctrl = pd.read_csv(f'{pw_out}/intermediate/control.count', sep='\t')
    case_log = np.log2(df_case + 1)
    ctrl_log = np.log2(df_ctrl + 1)
    df_delta = case_log - ctrl_log
    abs_values = np.abs(df_delta.values.flatten())
    cutoff = np.quantile(abs_values, 0.75) # verified, match with R code
    cutoff1 = round(cutoff, 1)
    df_delta = df_delta.clip(lower=-cutoff, upper=cutoff)
    
    # plot
    
    
    os.system(f'rm {pw_out}/intermediate/*.tmp')


def groseq_fisherexact_pausing_for_gene(c1, l1, c2, l2):
    """
    Perform Fisher's exact test to calculate the p-value for the pausing index of each gene. x is a row of the pausing count dataset
    per run, is about 1ms, so it is quite slow., 1-2ms per run
    user new fisher package, the speed is 100x faster (6us per call)
    must make sure the input are all non-zero
    """
    # if c1 == 0 or c2 == 0 or l1 == 0 or l2 == 0:
    #     return 'NA'
    # else:
        # null hypothesis: the pause read count is uniformly distributed in the gene body and the promoter region

    expectC1 = round(l1 * (c1 + c2) / (l1 + l2))
    expectC2 = (c1 + c2 - expectC1)
    # default = two-sided
    return fisher.pvalue(c1, expectC1, c2, expectC2).two_tail

def groseq_fisherexact_comparison_for_gene(tssCountGene1, gbCountGene1, tssCountGene2, gbCountGene2):
    """
    Perform Fisher's exact test to calculate the p-value for the comparison of pausing index between two genes or 2 conditions of the same gene
    """
    c1 = round(float(tssCountGene1))
    c2 = round(float(gbCountGene1))
    c3 = round(float(tssCountGene2))
    c4 = round(float(gbCountGene2))
    
    if c1 == 0 or c2 == 0 or c3 == 0 or c4 == 0:
        return None
    
    # default = two-sided
    return fisher.pvalue(c1, c2, c3, c4).two_tail

def calc_FDR(pvalues, with_na=True):
    """
    Calculate the false discovery rate (FDR) using the Benjamini-Hochberg method
    with_na, will first filter out the invalid pvalues, and then map the pvalues back to the original index
    """
    
    from statsmodels.stats.multitest import multipletests
    if not with_na:
        return multipletests(pvalues, method='fdr_bh')[1]
    
    pvalues_valid = []
    idx_map = {}
    idx_valid_p = -1
    for idx_old, p in enumerate(pvalues):
        if p is None or isinstance(p, str) or np.isnan(p):
            continue
        idx_valid_p += 1
        pvalues_valid.append(p)
        idx_map[idx_old] = idx_valid_p
    
    fdr = multipletests(pvalues_valid, method='fdr_bh')[1]  # return is a tuple, (reject, corrected p-values, alphaSidak, alphaBonf )
    fdr_all = ['NA' if idx not in idx_map else fdr[idx_map[idx]] for idx in range(len(pvalues))]
    return fdr_all


def run_deseq2(n_gene_cols, data, metadata, ref_level, col_group=None, min_reads=10, filter_out_low_expression=False, 
               size_factors_in=None, size_factor_only=False):
    """
    A count matrix of shape ‘number of samples’ x ‘number of genes’, containing read counts (non-negative integers),
    Metadata (or “column” data) of shape ‘number of samples’ x ‘number of variables’, containing sample annotations that will be used to split the data in cohorts.  rows = sample, columns = gene
    index should be samplename
    meta data, first col is the sample name (index), following columns usually condidion and group
    min_reads:  if a gene has less than min_reads in all samples, it will be removed
    ref_level,  the level used as reference in the comparison
    """
    
    
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    from pydeseq2.preprocessing import deseq2_norm_fit
    
    # data filtering
    
    counttable = data.iloc[:, n_gene_cols:]
    idx_gene = 1
    idx_transcript = 0
    idx_gene_cols = list(range(n_gene_cols))
    if n_gene_cols < 2:
        logger.error(f'DESeq analysis is for known genes only, the first 2 columns should be Transcript and Gene')
    gene_col = data.iloc[:, idx_gene] # use Transcript as index
    counttable.index = data.iloc[:, idx_transcript]
    
    counttable = counttable.T
    metadata.index = counttable.index
    
    # remove samples for which annotations are missing and exclude genes with very low levels of expression
    col_condition = 'condition'
    sample_to_keep = list(~metadata[col_condition].isna())
    counttable = counttable.loc[sample_to_keep, :]
    metadata = metadata.loc[sample_to_keep, :]
    
    if filter_out_low_expression:
        genes_to_keep = counttable.columns[counttable.sum(axis=0) >= min_reads]
        counttable = counttable[genes_to_keep]

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counttable,
        metadata=metadata,
        design_factors=col_condition, # compare samples based on condition
        ref_level=[col_condition, ref_level],
        refit_cooks=True,
        inference=inference,
    )
    if size_factors_in is None:  # gb change
        # tested, the size factors are exactly the same as the R deseq2 code
        dds.fit_size_factors()
        size_factors = dds.obsm['size_factors']
    else:
        size_factors = size_factors_in
        dds.obsm["size_factors"] = size_factors

    if size_factor_only:
        return None, size_factors
    
    # if each group only have one sample
    n_sam = counttable.shape[0]
    if n_sam == 2:
        log2fc = np.log2(data.iloc[:, n_gene_cols + 1] / size_factors[1]) / (data.iloc[:, n_gene_cols] / size_factor[0])
        res_df = data.iloc[:, idx_gene_cols]
        res_df['log2fc'] = log2fc
        return res_df, size_factors

    if size_factors_in is None:  # gb change
        dds.deseq2()
    else:
        counts = dds.X
        (dds.logmeans, dds.filtered_genes) = deseq2_norm_fit(counts)
        dds.layers["normed_counts"] = counts / size_factors[:, None]

        dds.fit_genewise_dispersions()
        dds.fit_dispersion_trend()
        dds.fit_dispersion_prior()
        dds.fit_MAP_dispersions()
        dds.fit_LFC()
        dds.calculate_cooks()
        dds.refit()

    # avialble attributes = X: count tdata, obs = design factors, 
    # obsm = sample matrix, keys = design_matrix, size_factors
    # varm = gene matrix, keys = dispersions, LFC

    # get stat
    stat_res = DeseqStats(dds, inference=inference)
    stat_res.summary()
    
    # the logFC is the same as R, but the p-value and padj are slightly different
    res_df = stat_res.results_df
    # here must convert gene_col to list, otherwise, the value will be all NaN, because the index does not match
    res_df.reset_index('Transcript', inplace=True)
    res_df.insert(1, 'Gene', list(gene_col))
    res_df = res_df.sort_values('padj')
    
    return res_df, size_factors


def filter_pp_gb(data, n_gene_cols, rep1, rep2, skip_filtering=False):
    """
    reused by change_pp_gb_with_case and change_pindex
    """
    n_sam = rep1 + rep2
    cols_raw = data.columns
    col_idx = {
        'ppc': {
            'cond1': list(range(n_gene_cols, n_gene_cols + rep1)),
            'cond2': list(range(n_gene_cols + rep1, n_gene_cols + n_sam)),
            'combined': list(range(n_gene_cols, n_gene_cols + n_sam)),
        },
        'gbc': {
            'cond1': list(range(n_gene_cols + n_sam,             n_gene_cols + n_sam + rep1 * 2, 2)),
            'cond2': list(range(n_gene_cols + n_sam + rep1 * 2, n_gene_cols + n_sam * 3, 2)),
            'combined': list(range(n_gene_cols + n_sam, n_gene_cols + n_sam * 3, 2)),
        },
        'gbd': {
            'cond1': list(range(n_gene_cols + n_sam + 1,             n_gene_cols + n_sam + rep1 * 2, 2)),
            'cond2': list(range(n_gene_cols + n_sam + rep1 * 2 + 1,  n_gene_cols + n_sam * 3, 2)),
            'combined': list(range(n_gene_cols + n_sam + 1, n_gene_cols + n_sam * 3, 2)),
        },
    }
    
    # cond1 = ctrl, cond2 = case
    idx_ppc_combined = col_idx['ppc']['combined']
    idx_gbc_combined = col_idx['gbc']['combined']
    idx_gbd_combined = col_idx['gbd']['combined']
    idx_gene_cols = list(range(n_gene_cols))
    
    
    if skip_filtering:
        data_pass = data
        data_drop = pd.DataFrame()
        return col_idx, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop
    
    tmp1 = sum([1 if cols_raw[i].startswith('gbc_') else 0 for i in idx_gbc_combined])
    tmp2 = sum([1 if cols_raw[i].startswith('gbd_') else 0 for i in idx_gbd_combined])
    if tmp1 != n_sam or tmp2 != n_sam:
        logger.error('the number of samples in gbc or gbd is not correct')
        sys.exit(1)
        return None
    
    sam_list = [re.sub('^ppc_', '', _) for _ in data.columns[idx_ppc_combined]]
    
    data['tmp_ppc_pass'] = data.iloc[:, idx_ppc_combined].apply(lambda x: all(x > 0), axis=1)
    
    for sn, idx_gbc_tmp, idx_gbd_tmp in zip(range(n_sam), idx_gbc_combined, idx_gbd_combined):
        sum_col = data.iloc[:, idx_gbc_tmp].sum()
        data[f'tmp_gbd_pass_{sn+1}'] = data.iloc[:, idx_gbd_tmp] * 10_000_000 / sum_col
    cols_gbd_pass = [f'tmp_gbd_pass_{sn+1}' for sn in range(n_sam)]
    data['tmp_gbd_pass'] = data[cols_gbd_pass].apply(lambda x: x.sum()/n_sam > 0.004, axis=1)
    
    
    # already compared with R result, the same
    data_pass = data[data['tmp_ppc_pass'] & data['tmp_gbd_pass']] # data0 in perl code
    data_drop = data.drop(data_pass.index) # data_0 in perl code

    data_pass = data_pass.drop(columns=[_ for _ in data_pass.columns if _[:4] == 'tmp_'])
    data_drop = data_drop.drop(columns=[_ for _ in data_drop.columns if _[:4] == 'tmp_'])

    return col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop

def change_pp_gb_with_case(n_gene_cols, rep1, rep2, data, out_dir, window_size, factor_flag, factor1=None, factor2=None):
    # data = processed, count_pp_gb.txt, already removed the chr, start, end  and strand column, and already collapsed the transcripts with the same gene
    # if factor1, factor2 is not None, it should be a list
    # n_gene_cols = 2 for known gene, and 1 for eRNA
    
    # columns = "Transcript\tGene"  + ppc_[sam_list], gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ....
    if rep2 == 0:
        factor2 = []
    
    n_sam = rep1 + rep2
    col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop = filter_pp_gb(data, n_gene_cols, rep1, rep2)

    data_pass.iloc[:, idx_gene_cols].to_csv(f'{out_dir}/intermediate/active_gene.txt', sep='\t', index=False, header=False)
    
    data_pass_pp = data_pass.iloc[:, idx_gene_cols + [idx_ppc_combined]] # datapp in R
    data_drop_pp = data_drop.iloc[:, idx_gene_cols + [idx_ppc_combined]] # data_pp in R
    
    data_pass_gb = data_pass.iloc[:, idx_gene_cols + [idx_gbc_combined]] # datagb in R
    data_drop_gb = data_drop.iloc[:, idx_gene_cols + [idx_gbc_combined]] # data_gb in R
    
    if rep2 == 0:
        # case only, the meta data is just the dummy value, won't be used for deseq, just for get the sizefactors
        metadata = pd.DataFrame({
            'condition': ['control'] + ['case'] * (rep1 - 1),
        })
    else:
        metadata = pd.DataFrame({
            'condition': ['control'] * rep1 + ['case'] * rep2,
        })
    ref_level = 'control'

    # with normalization factors
    if factor_flag == 1:
        factor1 = np.array(factor1)
        factor2 = np.array(factor2)
        norm_factors = np.concatenate([factor1, factor2])
        size_factors = 1 / norm_factors
    else:  # without normalization factors
        size_factors = None
    
    if rep2 > 0:
        # gb change
        with HiddenPrints():
            res_df, size_factors = run_deseq2(n_gene_cols, data_pass_gb, metadata, ref_level, size_factors_in=size_factors)
            res_df_full = pd.concat([res_df, data_drop_gb.iloc[:, idx_gene_cols]])
            res_df_full.to_csv(f'{out_dir}/known_gene/gb_change.txt', sep='\t', index=False)
            
            # pp change
            res_df, _ = run_deseq2(n_gene_cols, data_pass_pp, metadata, ref_level=ref_level, size_factors_in=size_factors)
            res_df_full = pd.concat([res_df, data_drop_pp.iloc[:, idx_gene_cols]])
            res_df_full.to_csv(f'{out_dir}/known_gene/pp_change.txt', sep='\t', index=False)
    elif rep1 == 1:
        size_factors = 1 # size factor is 1 for the only sample
    elif rep2 == 0:
        # control only
        _, size_factors = run_deseq2(n_gene_cols, data_pass_gb, metadata, ref_level, size_factors_in=size_factors, size_factor_only=True)

    # get the normalized data in other function, because we need the chr start end strand information, and the data here is already processed
    # save normalization factors
    nf = pd.DataFrame({'sample': sam_list, 'nfactor': norm_factors})
    nf.to_csv(f'{out_dir}/intermediate/nf.txt', sep='\t', index=False)

    return size_factors
    

def get_alternative_isoform_across_conditions(fn, out_dir, rep1, rep2):
    """
    get the alternative isoform across conditions
    1. get the genes with multiple isoforms
    2. for transcript, get the ppc sum of conditon1 and condition2,  compare the two ppc_sum. if the direction is different between 2 isoforms, then it is alternative isoform. e.g. in isoform1, ppc_sum_conditoin1 > ppc_sum_condition2, in isoform2, ppc_sum_conditoin1 < ppc_sum_condition2
    fn = count_pp_gb.txt / normalized_pp_gb.txt
    """
    
    if 'normalized' in fn:
        norm_flag = '.norm'
    else:
        norm_flag = ''
    data = pd.read_csv(fn, sep='\t')
    #  columns = "Transcript\tGene\tchr\tstart\tend\tstrand"  + ppc_[sam_list], gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ....
    # only keep the genes with multiple isoforms with different TSS
    data = data.groupby('Gene').filter(lambda x: x['Transcript'].nunique() > 1)
    data = data.groupby('Gene').filter(lambda x: x['start'].nunique() > 1) # 333 rows

    cols_orig = list(data.columns)
    cols_ppc_all = [i for i in cols_orig if i.startswith('ppc_')]

    sam_list = [re.sub('^ppc_', '', _) for _ in cols_ppc_all]
    sam1, sam2 = (sam_list[:rep1], sam_list[rep1:])
    
    cols_ppc1 = [f'ppc_{sam}' for sam in sam1]
    cols_ppc2 = [f'ppc_{sam}' for sam in sam2]
    
    data['ppc_sum1'] = data[cols_ppc1].sum(axis=1)
    data['ppc_sum2'] = data[cols_ppc2].sum(axis=1)
    # exclude transcripts with same ppc_sum1 and ppc_sum2
    data = data.loc[data['ppc_sum1'] != data['ppc_sum2']]
    
    # exclude transcripts with ppc_sum1 and ppc_sum2 are both < 5
    min_ppc_sum = 5
    data = data.loc[(data['ppc_sum1'] >= min_ppc_sum) | (data['ppc_sum2'] >= min_ppc_sum)]
    
    # keep only genes still have multiple isoforms
    data = data.groupby('Gene').filter(lambda x: x['Transcript'].nunique() > 1)
    data = data.groupby('Gene').filter(lambda x: x['start'].nunique() > 1) # 333 rows

    data = data.sort_values(['Gene', 'chr', 'start'])

    # for the comparison column, the value can be 1, 0 or -1
    # if ppc_sum1 > ppc_sum2, then 1, if ppc_sum1 < ppc_sum2, then -1, if ppc_sum1 == ppc_sum2, then 0
    data['comparison'] = data['ppc_sum1'].sub(data['ppc_sum2']).apply(lambda x: 1 if x > 0 else -1)
    # per gene, keep only the genes with different comparison values
    data = data.groupby('Gene').filter(lambda x: x['comparison'].nunique() > 1)
    
    # get the log2FC of ppc_sum1 and ppc_sum2, log2(ppc_sum2/ppc_sum1)
    ratio = data['ppc_sum2'] / (data['ppc_sum1'].replace(0, np.nan))
    data['log2FC'] = np.log2(ratio.replace(0, np.nan))
    
    # drop comparison column
    data = data.drop('comparison', axis=1)
    
    print('dumping result')
    data.to_csv(f'{out_dir}/known_gene/alternative_isoform{norm_flag}.txt', index=False, sep='\t')

def change_pp_gb(n_gene_cols, fn, out_dir, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0):
    """
    adapted from Rscript change_pp_gb.R
    if factor_flag == 1, then will use these normalization factors, otherwise , will use DESeq2 to normalize the data
    fn = count_pp_gb.txt

    """
    # input file = my $out_pp_gb = $inter_dir . "count_pp_gb.txt";
    # columns = "Transcript\tGene\tchr\tstart\tend\tstrand"  + ppc_[sam_list], gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ....
    err = 0
    rep2 = rep2 or 0
    
    if factor1 is not None and not isinstance(factor1, list):
        logger.error('factor1 should be a list')
        err = 1
    if factor2 is not None and not isinstance(factor2, list):
        logger.error('factor2 should be a list')
        err = 1
    if err:
        return 1
    
    data_raw = pd.read_csv(fn, sep='\t')
    data = data_raw.copy()
    data = data.iloc[:, list(range(n_gene_cols)) + list(range(n_gene_cols + 4, len(data.columns)))]
    cols = data.columns
    cols_ppc = [i for i in cols if i.startswith('ppc_')]
    cols_gbc = [i for i in cols if i.startswith('gbc_')]
    cols_gbd = [i for i in cols if i.startswith('gbd_')]
    
    data['ppc_sum'] = data[cols_ppc].sum(axis=1)
    data['gbc_sum'] = data[cols_gbc].sum(axis=1)
    
    ## !!! collapse the data from transcript level to gene level
    def is_max(col):
        max_v = col.max()
        return [1 if i == max_v else 0 for i in col]

    # colapase by ppc_sum
    data['SameGene_ppc_sum_max'] = data.groupby('Gene')['ppc_sum'].transform(is_max)
    data = data.loc[data['SameGene_ppc_sum_max'] == 1]
    
    # if multiple transcripts have the same max ppc_sum, then choose the one with the max gbc_sum
    data['SameGene_gbc_sum_max'] = data.groupby('Gene')['gbc_sum'].transform(is_max)
    data = data.loc[data['SameGene_gbc_sum_max'] == 1]
    
    # leaving the first one for each gene
    data = data.drop_duplicates(subset='Gene', keep='first').drop(columns=['SameGene_ppc_sum_max', 'SameGene_gbc_sum_max', 'ppc_sum', 'gbc_sum'])
    
    n_sam = rep1 + rep2
    
    # (rep1, rep2, data, out_dir, window_size, factor_flag, factor1=None, factor2=None)
    size_factors = change_pp_gb_with_case(n_gene_cols, rep1, rep2, data, out_dir, window_size, factor_flag, factor1, factor2)

    # get the normalized data
    n_extra_cols = 4 # chr, start, end, strand
    n_prev_cols = n_gene_cols + n_extra_cols
    idx_cols, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop = filter_pp_gb(data_raw, n_prev_cols, rep1, rep2, skip_filtering=True)
    
    norm_factors = 1 / size_factors
    ppc_norm = data_raw.iloc[:, idx_ppc_combined] * norm_factors
    ppd_norm = ppc_norm / window_size
    ppd_norm.columns = [f'ppd_{_}' for _ in sam_list]
    
    gbc_norm = data_raw.iloc[:, idx_gbc_combined] * norm_factors
    gbd_norm = data_raw.iloc[:, idx_gbd_combined] * norm_factors
    data_normalized = pd.concat([data_raw.iloc[:, :n_prev_cols], ppc_norm, ppd_norm, gbc_norm, gbd_norm], axis=1)
    data_normalized.to_csv(f'{out_dir}/known_gene/normalized_pp_gb.txt', sep='\t', index=False)


def cmhtest(row, rep1, rep2):
    """
    performs the Mantel-Haenszel test on a certain gene
    the row is like ppc_sam1, ppc_sam2, ppc_sam3, ppc_sam4, gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2, ...
    """
    import statsmodels as sm
    arr = []
    n_sam = rep1 + rep2
    for i in range(rep1):
        pp1 = row[i] # pp condition 1
        gb1 = row[n_sam + (i + 1) * 2 - 2] # gb condition 1
        for j in range(rep2):
            pp2 = row[rep1 + j]
            gb2 = row[n_sam + rep1 * 2 + (j + 1) * 2 - 2]
            # below is the gb2 defined in R code, which is wrong, here just used to test the consistency
            # gb2 = row[(rep2 + rep1 + j + 1) * 2 - 2]
            arr.append([[pp2, gb2], [pp1, gb1]])
    arr = np.array(arr).T
    cmt = sm.stats.contingency_tables.StratifiedTable(tables=arr)
    odds_ratio = cmt.oddsratio_pooled
    test_res = cmt.test_null_odds(correction=True)
    pvalue = test_res.pvalue
    statistic = test_res.statistic
    # print(f'odds_ratio = {odds_ratio}, pvalue = {pvalue}, statistic = {statistic}') 
    # return log2(odd_ratio), pvalue
    return np.log2(odds_ratio), pvalue


def get_pvalue_2_sample(row, idx_ppc1, idx_ppc2, idx_gbc1, idx_gbc2):
    """
    get the p-value for 2 samples
    using groseq_fisherexact_comparison_for_gene(tssCountGene1, gbCountGene1, tssCountGene2, gbCountGene2)
    """
    tssCountGene1 = row[idx_ppc1]
    gbCountGene1 = row[idx_gbc1]
    tssCountGene2 = row[idx_ppc2]
    gbCountGene2 = row[idx_gbc2]
    pvalue = groseq_fisherexact_comparison_for_gene(tssCountGene1, gbCountGene1, tssCountGene2, gbCountGene2)
    if pvalue is None:
        return 'NA'
    return pvalue
     
def change_pindex(fno_prefix, n_gene_cols, fn, out_dir, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0):
    """
    fno_prefix, default is empty for known genes, longeRNA- for longeRNA
    adapted from Rscript change_pindex.R
    input = count_pp_gb.txt  Transcript, Gene, chr, start, end, strand, (ppc_sam1 - ppc_sum2 ...), (gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ...)
    
    The transcripts were not collapsed to gene level
    """
    
    if fno_prefix not in ['', 'longeRNA-']:
        logger.error(f'invalid file prefix for pindex_change, should be empty or "longeRNA-", input="{fno_prefix}"')
        return 1
    
    is_erna = True if fno_prefix == 'longeRNA-' else False
    data = pd.read_csv(fn, sep='\t')
    n_sam = rep1 + rep2
    if is_erna:
        # the data does not contains the chr, start, end, strand columns, and there is no Gene column
        # no filtering needed
        ana_type = 'eRNA'
    else:
        ana_type = 'Known genes'
        n_cols_prev = n_gene_cols + 4
        cols_keep = list(range(n_gene_cols)) + list(range(n_cols_prev, len(data.columns))) # drop the chr, start, end, strand columns
        data = data.iloc[:, cols_keep]

    col_idx, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop = filter_pp_gb(data, n_gene_cols, rep1, rep2, skip_filtering=is_erna)
    
    if n_sam == 2:
        # use fisher exact test
        idx_ppc1 = col_idx['ppc']['cond1']
        idx_ppc2 = col_idx['ppc']['cond2']
        idx_gbc1 = col_idx['gbc']['cond1']
        idx_gbd1 = col_idx['gbd']['cond1']
        pvalues = data_pass.apply(lambda x: get_pvalue_2_sample(x, idx_ppc1, idx_ppc2, idx_gbc1, idx_gbd1), axis=1)
        fdr = calc_FDR(pvalues)
        # odds_ratio = (ppc2/ppc1) / (gbc2/gbc1) = ppc2 * gbc1 / (ppc1 * gbc2)
        log2fc = data_pass.apply(lambda x: np.log2(x[idx_ppc2] * x[idx_gbc1] / (x[idx_ppc1] * x[idx_gbc2])), axis=1)
        
        data_out = data_pass.iloc[:, idx_gene_cols]
        
        data_out['log2fc'] = log2fc
        data_out['pvalue'] = pvalues
        data_out['FDR'] = fdr

    # run the cmhtest, the data should have 3 dimensions
    else:
        # cmhtest(row, rep1, rep2)  # return = log2(odds_ratio), pvalue
        data_out = data_pass.iloc[:, idx_gene_cols]
        data_out[['log2fc', 'pvalue']] = data_pass.apply(lambda x: cmhtest(x, rep1, rep2), axis=1, result_type='expand')
        data_out['FDR'] = calc_FDR(data_out['pvalue'])
    
    data_out = data_out.sort_values('FDR')
    data_out_full = pd.concat([data_out, data_drop.iloc[:, idx_gene_cols]]).fillna('NA')
    data_out_full.to_csv(f'{out_dir}/known_gene/{fno_prefix}pindex_change.txt', sep='\t', index=False)
    logger.info(f'pindex_change done : {ana_type}')


def change_enhancer(fn, out_dir, condition, rep1, rep2, factor1=None, factor2=None, factor_flag=0):
    """
    fn = intermediate/count_enhancer.txt
    
    condition can be 1 or 2,  for 2, it is case-control design, and for 1, it is ctrl only
    """
    n_total_sam = rep1 + rep2
    fno_change = f'{out_dir}/eRNA/Enhancer_change.txt'
    fno_norm = f'{out_dir}/eRNA/normalized_count_enhancer.txt'
    
    # fn input header
    # Enhancer_ID	chr	start	end	sam1, sam2 ....
    
    # ctrl only
    if condition == 1:
        pass
    


def process_gtf(fn_gtf):
    """
    process the gtf file, get the gene start and end position, and other info
    gtf position is 1-idx, full closed
    """

    err = {'no_transcript_id': 0, 'no_gene_name': 0, 'invalid_line_format': 0, 'invalid_strand': 0}
    fn_gtf = os.path.abspath(fn_gtf)
    # check if have write permission to the folder of the gtf file
    if not os.access(os.path.dirname(fn_gtf), os.W_OK):
        home = os.path.expanduser("~")
        gtf_out_dir = f'{home}/.nrsa'
        os.makedirs(gtf_out_dir, exist_ok=True)
    else:
        gtf_out_dir = os.path.dirname(fn_gtf)
    
    if 'gtf' not in fn_gtf.rsplit('.', 2)[-2:]:
        logger.error(f'the gtf file should have the extension of .gtf: {fn_gtf}')
        sys.exit(1)
        return None, err
    
    fn_gtf_lb = os.path.basename(fn_gtf).replace('.gz', '').replace('.gtf', '')
    fn_gtf_pkl = f'{gtf_out_dir}/{fn_gtf_lb}.gtf_info.pkl'
    if os.path.exists(fn_gtf_pkl):
        logger.debug('loading gtf from pickle')
        with open(fn_gtf_pkl, 'rb') as f:
            return pickle.load(f), err
    
    gtf_col_idx = {
        'chr': 0,
        'source': 1,
        'feature': 2,
        'start': 3,
        'end': 4,
        'score': 5,
        'strand': 6,
        'attribute': 8,
    }
    res_raw = {}  # key = transcript_id, v = {'chr': '', 'strand': '', 'gene_id': '', 'gene_name': '', 'start': 0, 'end': 0}
    
    with gzip.open(fn_gtf, 'rt') if fn_gtf.endswith('.gz') else open(fn_gtf) as f:
        for i in f:
            if i[0] != '#' or not i:
                break
        for i in f:
            line = i.strip().split('\t')
            # chr1    hg19_ncbiRefSeq    exon    66999252    66999355    0.000000    +    .    gene_id "NM_001308203.1"; transcript_id "NM_001308203.1"; gene_name "SGIP1";
            # 1       ensembl_havana  gene    11869   14412   .       +       .       gene_id "ENSG00000223972"; gene_version "4"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";
            line_err = None
            try:
                chr_, start, end, strand, attribute_raw = [line[_] for _ in [gtf_col_idx[_] for _ in ['chr', 'start', 'end', 'strand', 'attribute']]]
                start = int(start)
                end = int(end)
            except:
                line_err = 'invalid_line_format'
            if strand not in {'+', '-'}:
                line_err = 'invalid_strand'
            chr_ = refine_chr(chr_)
            if '_' in chr_:
                continue # skip the line if the chr_ contains underscore (e.g. chr1_KI270706v1_random)

            attributes = attribute_raw.rstrip(';').split(';')
            gene_name = transcript_id = None
            
            for att in attributes:
                try:
                    k, v = att.strip().split(' ', 1)
                except:
                    continue
                k = k.strip('"')
                if k == 'transcript_id':
                    transcript_id = v.strip('"')
                elif k == 'gene_name':
                    gene_name = v.strip('"')
            if attribute_raw is None:
                line_err = 'no_transcript_id'
            if gene_name is None:
                line_err = 'gene_name_not_found'
            if line_err:
                err[line_err] += 1
                continue
            
            key_uniq = f'{chr_}@{strand}@{start}@{end}'
            # {'chr': '', 'strand': '', 'gene_name': '', 'start': 0, 'end': 0}
            res_raw.setdefault(gene_name, {}).setdefault(key_uniq, {'transcript': [], 'chr': chr_, 'strand': strand, 'gene_name': gene_name, 'start': start, 'end': end})
            res_raw[gene_name][key_uniq]['transcript'].append(transcript_id)
            
            ires = res_raw[gene_name][key_uniq]
            if chr_ != ires['chr']:
                continue

            if start < ires['start']:
                ires['start'] = start
            if end > ires['end']:
                ires['end'] = end

    # merge transcripts with same start and end position
    res = {}
    n_merged = 0
    for gn, v1 in res_raw.items():
        for v2 in v1.values():
            ntmp = len(v2['transcript'])
            if ntmp > 1:
                n_merged += ntmp - 1
            transcript_id = ';'.join(v2['transcript'])
            v2['transcript'] = transcript_id
            res[transcript_id] = {k: v2[k] for k in ['chr', 'strand', 'gene_name', 'start', 'end']}
    
    
    logger.warning(f'merged {n_merged} transcripts with same start and end position')
    
    with open(fn_gtf_pkl, 'wb') as o:
        pickle.dump(res, o)

    err_total = sum(err.values())
    if err_total:
        logger.info(f'error in parsing gtf file: {err}')
    return res, err



def check_dependency():
    """
    check if the dependency: HOMER and bedtools are installed
    """
    # check HOMER
    err = 0

    retcode = run_shell('which makeTagDirectory')
    if retcode != 0:
        logger.error("HOMER is not installed, please install HOMER first")
        err = 1

    # check bedtools
    retcode = run_shell('which bedtools')
    if retcode != 0:
        logger.error("bedtools is not installed, please install bedtools first")
        err = 1
    
    # check python packages
    required_packages = ['pydeseq2', 'pandas', 'numpy', 'statsmodels', 'scipy', 'matplotlib', 'fisher']
    import importlib.util
    missing_python_packages = []
    for pkg in required_packages:
        if importlib.util.find_spec(pkg) is None:
            missing_python_packages.append(pkg)
            err = 1
    if missing_python_packages:
        logger.error(f"Missing python packages: {missing_python_packages}, please install them first.")
    
    return err


def run_shell(cmd):
    """
    run shell command and check the return code and stdout stderr
    """
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    retcode = p.returncode
    if retcode:
        logger.error(f"Error running command: {cmd}")
        logger.error(f"stdout: {stdout.decode()}")
        logger.error(f"stderr: {stderr.decode()}")
    return retcode

def get_lb(fn):
    return fn.rsplit('/', 1)[-1].rsplit('.', 1)[0]



def pre_count_for_bed(fh_bed, fn_lb, pwout):
    """
    fh_bed can be the real file handle or a subprocess.Popen.stdout object, so do not use seek or tell
    process the bed file, and get the count of the read end regarding strand.
    will export 2 files, 
    1. dict, key = chr, k2 = strand, k3 = chunk_number. value = count of reads with read end in this chunk.  the chunk_number is got from position // chunk_size  (default chunksize = 10bp)
    2. bed file, record the read end chr, read_end position, strand and the count of reads with this read end position
    """
    chunk_size = 10
    idx_bin_size = 50 # for the per position read end count file (final output), will build the index

    fn_out = f'{pwout}/intermediate/bed/{fn_lb}.count.per_base.txt'
    fn_out_idx = f'{pwout}/intermediate/bed/{fn_lb}.count.per_base.idx.pkl'
    fn_count_bin = f'{pwout}/intermediate/bed/{fn_lb}.count.bin_of_{chunk_size}.pkl'
    
    
    if os.path.exists(fn_out) and os.path.getsize(fn_out) > 10 and os.path.exists(fn_out_idx) and os.path.exists(fn_count_bin):
        logger.info(f'bed file already pre-counted: {fn_lb}')
        with open(fn_out_idx, 'rb') as f:
            idx_d = pickle.load(f)
        with open(fn_count_bin, 'rb') as f:
            res_chunk = pickle.load(f)
        return fn_lb, fn_out, idx_d, res_chunk

    if fh_bed is None:
        return None
    res_chunk = {'bin_size': chunk_size}
    res_bed = {}
    idx_d = {'bin_size': idx_bin_size}  # key = chr, k2 = bin_number, value = byte position in the file
    file_out_pos = 0 # current byte position in the output count file

    
    def dump_pos_to_file(chr_, res_bed_chr, fo, file_out_pos):
        dump_pos_l = sorted(res_bed_chr)
        for dp in dump_pos_l:
            ct = res_bed_chr.pop(dp)
            chunk = dp // idx_bin_size
            idx_tmp = idx_d.setdefault(chr_, {})
            line_out = f'{chr_}\t{dp}\t{ct[0]}\t{ct[1]}\n'
            if chunk not in idx_tmp:
                idx_tmp[chunk] = file_out_pos
            fo.write(line_out)
            file_out_pos += len(line_out)
        return file_out_pos

    n = 0
    section_size, unit = 1_000_000, 'M'
    with open(fn_out, 'w') as fo:
        for i in fh_bed:
            # chr1	10511	10569	A00758:60:HNYK7DSXX:4:2461:27181:23171	32	-

            n += 1
            if n % section_size == 0:
                logger.debug(f'Processing {n/section_size:.0f}{unit} lines')
            chr_, s, e, read_id, _, strand = i[:-1].split('\t')
            if strand == '+':
                read_end = int(e)
            else:
                read_end = int(s) + 1

            chunk = read_end // chunk_size
            res_chunk.setdefault(chr_, {}).setdefault(strand, {}).setdefault(chunk, 0)
            res_chunk[chr_][strand][chunk] += 1
            
            res_bed.setdefault(chr_, {}).setdefault(read_end, [0, 0])  # for each read end, the count of + and - strand
            idx_strand = 0 if strand == '+' else 1
            res_bed[chr_][read_end][idx_strand] += 1

        # refine the chr_
        res_bed = {refine_chr(k): v for k, v in res_bed.items()}
        res_chunk = {refine_chr(k): v for k, v in res_chunk.items()}

        for chr_ in sorted(res_bed):
            # dump the res_bed_chr to file
            logger.debug(f'{fn_lb} - dumping the pos list for chr {chr_}')
            res_bed_chr = res_bed[chr_]
            if chr_[:2].lower() == 'gl':
                continue
            file_out_pos = dump_pos_to_file(chr_, res_bed_chr, fo, file_out_pos)
        
    # dump the pickle
    with open(fn_out_idx, 'wb') as f:
        pickle.dump(idx_d, f)
    
    with open(fn_count_bin, 'wb') as f:
        pickle.dump(res_chunk, f)
    return fn_lb, fn_out, idx_d, res_chunk


def process_input(pwout, fls, bed_idx_step=100):
    """
    process input files, convert bam to bed if needed, and build the index for bed files
    return a list of [fn_bed, idx_bed, fn_lb]
    """
    if fls is None or not fls:
        return None
    err = 0
    res = []
    lb_map = {}
    for fn in fls:
        fn_lb = get_lb(fn)
        fn_out = f'{pwout}/intermediate/bed/{fn_lb}.count.per_base.txt'

        if os.path.exists(fn_out) and os.path.getsize(fn_out) > 10:
            tmp = pre_count_for_bed(None, fn_lb, pwout)
            if tmp is not None:
                res.append(tmp)
                continue
        if fn.endswith('.bam'):
            logger.info(f'Converting {fn}')
            cmd = f"""bedtools bamtobed -i {fn}"""
            with Popen(cmd, shell=True, stdout=PIPE, text=True) as p:
                with p.stdout as fh_bed:
                    res.append(pre_count_for_bed(fh_bed, fn_lb, pwout))
        elif fn.endswith('.bed'):
            with open(fn) as fh_bed:
                res.append(pre_count_for_bed(fh_bed, fn_lb, pwout))
        else:
            logger.error(f"Input file '{fn}' should be in bed or bam format")
            err = 1
            continue
    if err:
        return None
    return res


def process_tss_tts(fn):
    """
    example, /Users/files/work/jb/work/NRSA_v2/new/ref/hg19/RefSeq-hg19-tss-tts.txt
    chr1	11874	14409	NR_046018.2	+
    """
    tss_tts = {}
    with open(fn, 'r') as file:
        for line in file:
            temp = line.strip().split('\t')
            tss_tts[temp[3]] = {
                'chr': temp[0],
                'tss': temp[1],
                'tts': temp[2],
                'strand': temp[4]
            }

    return tss_tts

def filter_by_tss_tts(tss_tts_info, line):
    """
    filter the line by tss_tts_info, return True if the line is in the tss_tts_info, otherwise, False
    input line is enhancer out or lerna out
    """
    flag = 0
    temp = in_list[i].split('\t')
    chr = temp[0]
    start = temp[1]
    end = temp[2]
    
    for key in gene:
        if chr != gene[key]['chr']:
            continue
        
        if gene[key]['strand'] == '+':
            if start <= gene[key]['tts'] and (gene[key]['tss'] - end) < filter_tss:
                flag = 1
                break
            if end >= gene[key]['tss'] and (start - gene[key]['tts']) < filter_tts:
                flag = 1
                break
        elif gene[key]['strand'] == '-':
            if start <= gene[key]['tss'] and (gene[key]['tts'] - end) < filter_tts:
                flag = 1
                break
            if end >= gene[key]['tts'] and (start - gene[key]['tss']) < filter_tss:
                flag = 1
                break
    
    if flag == 0:
        filter_outstr += in_list[i] + "\n"

