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
import traceback
# import bisect
from statsmodels.stats.multitest import multipletests
import statsmodels.stats.contingency_tables as contingency_tables
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
    chr_ = chr_.strip().lower()
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

def dummy(idx, fh, chr_, start, end):
    return 'ATCG'

def get_seqence_from_fa(idx, fh, chr_, start, end):
    """
    idx is from loading of the json {fn_fa}.idx.json file, fh is the opened file handle of the fasta file
    already converted to uppercase
    for 20k transcripts, get the whole gene sequence cost about 2 second
    """
    # get_seqence_from_fa.note = 'normal'
    
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


def check_is_sorted(fn_bed):
    """
    check if the bed file is sorted or not, 
    """
    chr_prev = None
    chr_set = set()
    logger.info(f'Checking if bed file is sorted: {fn_bed}')
    with gzip.open(fn_bed, 'rt') if fn_bed.endswith('.gz') else open(fn_bed) as f:
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
                not_found.append(f'{k}: {fn}')
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


def get_peak_method1(count_per_base, chr_, strand, s, e):
    # sum from the per-base pre-count dict
    ct_sum = 0
    strand_idx = 0 if strand == '+' else 1
    res_chr = count_per_base[chr_]
    return sum([res_chr.get(i, [0, 0])[strand_idx] for i in range(s, e + 1)])

def get_peak_method2(count_per_base, count_bin, chr_, strand_idx, s, e, bin_size):
    # compared to get_peak_method1 function
    # test a 1.5kb region to get the mapped reads count
    # method2, combination of per-site count and bin_collected count (bin_size = 10), cost 19us 
    # method1, use the per-site count only, cost 150us
    # so the method2 is 7 times faster than method1
    # for small window, e.g 50bp, method2 is still faster, 5.8us vs 10.8us
    # will also test if will improve when the pre_count bin_size is 20bp
    # for bin_size = 20, speed is 12.5 and 5.4us for 1.5kb and 50bp window, both faster than bin_size = 10
    
    
    # res is outer dict, stores the per_site count
    # bin_size = count_bin['bin_size']
    bin_start = s // bin_size
    bin_end = e // bin_size
    if bin_end - bin_start < 2:
        return sum([count_per_base[chr_].get(i, [0, 0])[strand_idx] for i in range(s, e + 1)])
    # the region is covered by 3 or more bins
    points = []
    left_mod = s % bin_size
    right_mod = e % bin_size
    if left_mod:
        points += range(s, (bin_start + 1) * bin_size)
    if right_mod != bin_size - 1:
        bin_list = range(bin_start + 1, bin_end)
        points += range(bin_end * bin_size, e + 1)
    else:
        bin_list = range(bin_start + 1, bin_end + 1)
    return sum([count_per_base[chr_].get(i, [0, 0])[strand_idx] for i in points]) + sum([count_bin[chr_].get(i, [0, 0])[strand_idx] for i in bin_list])

def get_peak(count_per_base, count_bin, chr_, strand, gene_raw_s, strand_idx, pp_start, pp_end, gb_start, gb_end, gb_len_mappable, gene_seq, window_size, step_size, bin_size, prev_pp_peak):
    # {'chr': '5', 'strand': '-', 'gene_name': 'PFDN1', 'start': 139682626, 'end': 139682689}
    # pro_up = pro_down = 500
    # gb_start_pos = 1000 # count the gb density from 1000bp downstream of the gene start site
    # using the already_parsed_windows will cost more time

    gbc = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, gb_start, gb_end, bin_size=bin_size)
    if gbc == 0:
        gbd = 0
    else:
        gbd = gbc / gb_len_mappable
    # gbc = gbd = 0
    # k_pp_region = f'{chr_}{strand}:{pp_start}'
    
    
    # get all the sites count in pp region
    pp_region_count = [count_per_base[chr_].get(i, [0, 0])[strand_idx] for i in range(pp_start, pp_end + 1)]
    window_ct = [(i, sum(pp_region_count[i: i + window_size])) for i in range(0, pp_end - pp_start - window_size + 2, step_size)]
    peak_window_start, ppc  = max(window_ct, key=lambda _: _[1]) 
    peak_window_ct = [(i, pp_region_count[i]) for i in range(peak_window_start, peak_window_start + window_size)]
    summit_pos, summit_count = max(peak_window_ct, key=lambda _: _[1])
    summit_pos += pp_start
    summit_pos_str = f'{chr_}{strand}:{summit_pos}'
    
    if ppc == 0:
        ppd = 0
    elif gene_seq:
        peak_windows_start_abs = peak_window_start + pp_start - gene_raw_s
        mappable_sites = window_size - gene_seq[peak_windows_start_abs: peak_windows_start_abs + window_size].count('N')
        ppd = ppc / mappable_sites
    else:
        ppd = ppc / window_size

    pp_res = {'ppc': ppc, 'ppd': ppd, 'mappable_sites': window_size, 'summit_pos': summit_pos_str, 'summit_count': summit_count}
    # prev_pp_peak[k_pp_region] = pp_res
    # gbc = gbd = 0
    # return 0, 0, {'ppc': 0, 'ppd': 0, 'mappable_sites': 50, 'summit_pos': 0, 'summit_count': 0}
    return gbc, gbd,  pp_res


def get_mapped_reads(count_file_fh, count_file_idx, chr_, strand, start, end):
    """
    get the coverage of the region, the coverage is the number of **read ends** at each position in this region
    each call will cost around 100us.
    """
    bin_size = count_file_idx['bin_size']
    start, end = sorted([start, end])
    chunk_s = start // bin_size
    chunk_e = end // bin_size
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
    bin_size = count_file_idx['bin_size']
    start, end = sorted([start, end])
    chunk_s = start // bin_size
    chunk_e = end // bin_size
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


def draw_box_plot(n_gene_cols, pwout, out_name, n_rep1, n_rep2=None, gn_list=None):
    """
    gn_list, the genes to be plotted, if None, all genes will be plotted
    pwout: folder name for pause_PROseq.pl outputs
    n_rep1 and n_rep2:  the replacates number for condition1 and conditon2
    n_rep2 is optional
    """
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    n_rep2 = n_rep2 or 0
    gn_list = set(gn_list) if gn_list else None
    # check the pause_PROseq.pl output files
    for i in ['known_gene', 'intermediate']:
        ipw = f'{pwout}/{i}'
        if not os.path.exists(ipw):
            logger.error(f'output folder of pause_PROseq.py not found: {i}')
            return 1

    # read pp_gb results
    fn_pp_gb = f'{pwout}/known_gene/normalized_pp_gb.txt'
    fn_pindex = f'{pwout}/known_gene/pindex.txt'
    fno = f'{pwout}/intermediate/boxplot_data.txt'
    
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
        idx_list += [i for i, col in enumerate(header_in) if col[:4] == 'ppd_']
        idx_list += [i for i, col in enumerate(header_in) if col[:4] == 'gbd_']
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
        header_str += [header_in[_] for _ in idx_list]
        
        for i in f:
            line = i[:-1].split('\t')
            gn = line[0]
            if gn_list is None or gn in gn_list:
                k = f'{gn}@{line[1]}'
                if k not in out_str:
                    logger.warning(f'{gn} in pindex not found in pp_gb')
                    continue
                ires = out_str[k]
                try:
                    vals_pindex = [line[_] for _ in idx_list]
                except:
                    logger.error(line)
                    sys.exit(1)
                if '' in vals_pindex:
                    logger.error(line)
                    sys.exit(1)
                ires += vals_pindex
                

    # write to file
    with open(fno, 'w') as o:
        o.write('\t'.join(header_str) + '\n')
        for k in sorted(out_str):
            o.write('\t'.join(out_str[k]) + '\n')
    
    # draw the box plot, adapted from boxplot.R
    # R CMD --args outdir=\"$inter_dir\" pname=\"$name\" custom=1 rep1=$rep1 rep2=$rep2
    df_box = pd.read_csv(fno, sep='\t')
    cols_pp_density = [_ for _ in df_box.columns if _.startswith('ppd_')]
    cols_gb_density = [_ for _ in df_box.columns if _.startswith('gbd_')]
    cols_pindex = [_ for _ in df_box.columns if _.endswith('pindex')]
    sam_list = [name.split('ppd_')[1] for name in cols_pp_density]
    # logger.info(df_box.head())
    
    pp_density = df_box.loc[:, cols_pp_density]
    gb_density = df_box.loc[:, cols_gb_density]
    pindex =     df_box.loc[:, cols_pindex]
    
    # logger.info(pindex.head(20))

    pp_density.columns = sam_list
    gb_density.columns = sam_list
    pindex.columns = sam_list
    
    # plot
    fn_plot_pp = f'{pwout}/known_gene/{out_name}_pp_density.pdf'
    fn_plot_gb = f'{pwout}/known_gene/{out_name}_gb_density.pdf'
    fn_plot_pindex = f'{pwout}/known_gene/{out_name}_pausing_index.pdf'
    
    def plot_task(fn, ylabel, data, factor=1):
        with PdfPages(fn) as pdf:
            plt.figure(figsize=(2, 4))
            
            # Split the data into two groups based on n_rep1
            data_rep1 = data.iloc[:, :n_rep1]
            data_rep1 = [data_rep1[_].dropna() * factor for _ in data_rep1]
            
            if n_rep2:
                data_rep2 = data.iloc[:, n_rep1:]
                data_rep2 = [data_rep2[_].dropna() * factor for _ in data_rep2]
            
            # previously, data_rep1 and data_rep2 were all x 1000 as the input of box plot data
            
            # Plot the boxplot for the first n_rep1 samples (blue color)
            # plt.boxplot(data_rep1 * 1000, positions=range(1, n_rep1 + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='blue'), capprops=dict(color='blue'), whiskerprops=dict(color='blue'), medianprops=dict(color='blue'))
            plt.boxplot(data_rep1, positions=range(1, n_rep1 + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='blue'), capprops=dict(color='blue'), whiskerprops=dict(color='blue'), medianprops=dict(color='blue'))
            
            # Plot the boxplot for the rest of the samples (red color)
            if n_rep2:
                # plt.boxplot(data_rep2 * 1000, positions=range(n_rep1 + 1, n_total_sam + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='red'), capprops=dict(color='red'), whiskerprops=dict(color='red'), medianprops=dict(color='red'))
                plt.boxplot(data_rep2, positions=range(n_rep1 + 1, n_total_sam + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='red'), capprops=dict(color='red'), whiskerprops=dict(color='red'), medianprops=dict(color='red'))
            
            plt.ylabel(ylabel)
            plt.xticks(range(1, n_total_sam + 1), sam_list, rotation='vertical')
            plt.tight_layout()
            pdf.savefig()
            
    plot_task(fn_plot_pp, 'Reads per Kb (RPK)', pp_density, 1000)
    plot_task(fn_plot_gb, 'Reads per Kb (RPK)', gb_density, 1000)
    plot_task(fn_plot_pindex, 'Pausing Index', pindex)

    def plot_hist(n_sam, n_prev_col, ppd, gbd, condition_sn):

        plt.figure(figsize=(12, 6 * (n_sam - 1)))
        for k in range(1, n_sam):
            col1 = n_prev_col
            col2 = k + n_prev_col
            col1_str = ppd.columns[col1]
            col2_str = ppd.columns[col2]
            xlabel = f"log2({sam_list[0]}/{sam_list[k]})"
            plt.subplot((n_sam - 1), 2, k)
            ppd_tmp = ppd.loc[(ppd.iloc[:, col1] > 0) & (ppd.iloc[:, col2] > 0), [col1_str, col2_str]]
            plt.hist(np.log2(ppd_tmp.iloc[:, 0] / ppd_tmp.iloc[:, 1]), bins=100)
            plt.title(f"ppd: rep1vs.rep{k + 1}")
            plt.xlabel(xlabel)

            plt.subplot((n_sam - 1), 2, k + 1)
            gbd_tmp = gbd.loc[(gbd.iloc[:, col1] > 0) & (gbd.iloc[:, col2] > 0), [col1_str, col2_str]]
            plt.hist(np.log2(gbd_tmp.iloc[:, 0] / gbd_tmp.iloc[:, 1]), bins=100)

            plt.title(f"gbd: rep1vs.rep{k + 1}")
            plt.xlabel(xlabel)

        plt.savefig(f"{pwout}/known_gene/Reps-condition{condition_sn}.tif")
        plt.close()
    
    if n_rep1 > 1:
        plot_hist(n_rep1, 0, pp_density, gb_density, 1)
    if n_rep2 > 1:
        plot_hist(n_rep2, n_rep1, pp_density, gb_density, 2)
    

def draw_heatmap_pindex(n_gene_cols, pwout):
    """
    plot for the pindex change
    """
    from matplotlib import pyplot as plt
    import matplotlib.colors as mcolors
    
    fn_pindex_change = f'{pwout}/known_gene/pindex_change.txt'
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

    img = ax.imshow(data_plot.iloc[:, 2].values.reshape(-1, 1), cmap=my_palette, aspect='auto')

    plt.savefig(f'{pwout}/known_gene/pindex_change.pdf', bbox_inches='tight')
    plt.close()
def get_line_count(fn):
    n = 0
    with gzip.open(fn) if fn.endswith('.gz') else open(fn) as f:
        for i in f:
            n += 1
    return n

def calculate_signal(fn_lb, fn_bed, fn_split_bin_bed, fn_coverage):
    s = time.time()
    cmd = f'bedtools coverage -a {fn_split_bin_bed} -b {fn_bed} -counts -sorted > {fn_coverage}'
    retcode = run_shell(cmd)
    if retcode:
        logger.error(f'failed to run bedtools coverage for {fn_lb}')
        return 1
    # process the coverage file
    # chr1	323932	323952	1   10
    logger.debug(f'bedtools coverage for {fn_lb}, done, time used = {dur:.2f}s')
    count = {}
    with open(fn_coverage) as f:
        for i in f:
            line = i[:-1].split('\t')
            bin_sn, ict = line[3], line[4]
            count.setdefault(bin_sn, 0)
            count[bin_sn] += int(ict)
    return count

def draw_signal(pwprj, fn_enhancer_center, fls_ctrl,fls_case, distance=2000, bin_size=20, signal_type='p'):
    # fn_enhancer_center #enhancer central position file, the output of eRNA.pl (Enhancer_centralposition.txt); required, enhancer central position file, the output of eRNA.pl (Enhancer_center.txt), or the file of chromation location of interest (should contains at least two columns separated by tab or space, e.g "chr11	83078317")
    # example line = chr1	565340	565340	+, no header
    
    # my $distance = 2000;	#up and down relative distance from peak centre;
    # my $bin_size=20; #bin size for smoothing;
    # my @case_bed;	#read alignment files in bed/bam format (case);
    # my @control_bed;	#read alignment files in bed/bam format (control);
    # my $signal="p"; The signal type (-s) should be either 'p' or 'o', p - signal for PROseq was input in pause_PROseq.pl; o - signal for other data such as histone modification (default: p)
    fls_case = fls_case or []
    if signal_type not in {'p', 'o'}:
        logger.error(f'invalid signal type: {signal_type}, valid = "p" or "o"')
        return 1
    pwout = f'{pwprj}/eRNA'
    pwinter = f'{pwout}/intermediate'
    fn_padding = f'{fn_enhancer_center}.padding.tmp'
    fn_split_bin_bed = f'{fn_enhancer_center}.split_regions.bed'

    n_valid_sites = 0
    n_bad_line = 0 # the line which 2nd column is not numeric
    min_dist_before = distance + bin_size/2
    nfactor = {}  # key = fn_lb, v = total lines in the bed file, the line_count should be calculated during the pre-counting step
    res = {}
    half_bin_size = int(bin_size / 2)
    
    fn_bed = fls_ctrl[0][1]
    fn_chr_map = fn_bed.replace('.sorted.bed', '.chr_map.pkl')
    with open(fn_chr_map, 'rb') as f:
        chr_map = pickle.load(f)

    with open(fn_enhancer_center) as f, open(fn_padding, 'w') as o:
        for i in f:
            try:
                line = i[:-1].split('\t')
                pos = int(line[1])
                if pos > min_dist_before:
                    n_valid_sites += 1
                    chr_ = chr_map.get(line[0], line[0])
                    left_boundary = pos_int - region_size - half_bin_size
                    right_boundary = pos_int + region_size + half_bin_size
                    print(f'{chr_}\t{left_boundary}\t{right_boundary}', file=o)
            except:
                n_bad_line += 1
    if n_bad_line:
        logger.warning(f'{n_bad_line} lines in {fn_enhancer_center} are invalid, 2nd column is not numeric')

    # get the line count for each bed file
    for fn_lb, fn_bed in fls_ctrl + fls_case:
        fn_line_count = fn_bed.replace('.sorted.bed', '.line_count.txt')
        if os.path.exists(fn_line_count):
            with open(fn_line_count) as f:
                nfactor[fn_lb] = int(f.readline())
        else:
            line_count = get_line_count(fn)
            nfactor[fn_lb] = line_count
            with open(fn_line_count, 'w') as o:
                o.write(str(line_count))

    # build the regions around the enhancer center for plot
    # no region name, the final output from makewindows is like
    # chr1	20	30	3  the 4th column is the region number
    def makewindows(fn_lb, fn_tss_padding, bin_size, chr_map):
        # rewrite the bedtools makewindows, because the original one does not contain the strand information
        fno = f'{pwout}/intermediate/{fn_lb}.split_regions.bed'
        with open(fn_tss_padding) as f, open(fno, 'w') as o:
            for i in f:
                chr_, s, e, strand = i[:-1].split('\t')
                bin_sn = 0
                chr_orig = chr_map.get(chr_, chr_)
                for pos in range(s, e, bin_size):
                    bin_sn += 1
                    print(f'{chr_orig}\t{pos}\t{pos + bin_size}\t{bin_sn}\t{strand}', file=o)
        return fno
    
    # build the new tss for the first file, and use it  for the rest
    fn_lb = fls_ctrl[0][0]
    fn_chr_map = f'{pw_bed}/{fn_lb}.chr_map.pkl'
    with open(fn_chr_map, 'rb') as f:
        chr_map = pickle.load(f)
    fn_split_bin_bed = makewindows(fn_lb, fn_padding, bin_size, chr_map)
    

    # get the signal results
    for fn_lb, fn_bed in fls_ctrl + fls_case:
        fn_coverage = f'{pwinter}/{fn_lb}.signal.coverage.txt'
        res[fn_lb] = calculate_signal(fn_lb, fn_bed, fn_split_bin_bed, fn_coverage)
    
    # save nfactor for type = 'o'
    if signal_type == 'o':
        with open(f'{pwinter}/nf_o.txt') as o:
            print('sample\tnfactor', file=o)
            for fn_lb, line_count in nfactor.items():
                print(f'{fn_lb}\t{10_000_000/line_count}', file=o)

    # export the signal results
    fn_pro_signal = f'{pwinter}/PROseq_signal.txt'
    with open(fn_pro_signal, 'w') as o:
        print('position\t' + '\t'.join([fn_lb for fn_lb, _ in fls_ctrl + fls_case]), file=o)
        bin_sn = 0
        for bin_pos in range(-distance, distance + 1, bin_size):
            bin_sn += 1
            row = [str(bin_pos)]
            for fn_lb, _ in fls_ctrl + fls_case:
                row.append(f'{res[fn_lb][bin_sn]/n_valid_sites:.4f}')
            print('\t'.join(row), file=o)
    

    # plot the signal
    # signal = 1 if signal_type == 'o' else 0
    # file=fn_pro_signal  outdir=\"$inter_dir\" sample=\"$samplename_str\" pname=\"$name\" rep1=$rep1 rep2=$rep2 signal=1 \' $cmd_file

def draw_heatmap_pp_change(n_gene_cols, pwout, pw_bed, fls_ctrl, fls_case, fn_tss, region_size=5000, bin_size=200, outname='heatmap'):
    """
    draw heatmap for pp change
    fls_ctrl, fls_case, foramt is like [fn_lb, fn_bed]
    region_size, upstream and downstream distance relative to TSS for plotting PRO-seq signal (bp, default: 5000),should can be divided by bin size
    bin_size (bp, default: 200),should can be divided by 2
    """
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    fls_case = fls_case or []
    
    ppchange = {}
    fn_pp_change = f'{pwout}/known_gene/pp_change.txt'
    fn_out_pdf = f'{pwout}/known_gene/heatmap.pdf'
    # fn_data_bed = f'{pwout}/intermediate/data_bed.tmp'
    fn_tss_padding = f'{pwout}/intermediate/infile_bed.tmp' # expand the fn_tss with upstream and downstream region by region_size
    time_bedtools_coverage = 0

    if bin_size % 2 or region_size % bin_size:
        logger.error(f'bin_size should be even and region_size should be divisible by bin_size')
        return 1

    bin_number = int(region_size * 2 / bin_size + 1)  # *2 because of the upstream and downstream

    def to_number(s):
        if not s or s.upper() == 'NA':
            return 'NA'
        return float(s)

    # read pp_change results
    # Transcript	Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
    n_case, n_ctrl = len(fls_case), len(fls_ctrl)
    
    # gene_list is the genes with non-NA log2FoldChange in the pp_change.txt
    with open(fn_pp_change) as f:
        header = f.readline().strip().split('\t')
        logger.debug(header)
        idx_transcript = header.index('Transcript')
        idx_logfc = header.index('log2FoldChange')
        for i in f:
            line = i[:-1].split('\t')
            transcript_id = line[idx_transcript]
            logfc = to_number(line[idx_logfc])
            if logfc != 'NA':
                ppchange[transcript_id] = logfc

    # tss file
    # chr1	11874	11874	NR_046018.2	+
    # tss_res = []
    tss_with_padding_n = 0
    strand_info = {}
    half_bin_size = int(bin_size / 2)
    with open(fn_tss) as f, open(fn_tss_padding, 'w') as o:
        for i in f:
            line = i.strip().split('\t')
            chr_, pos, transcript_id, strand = [line[_] for _ in [0, 1, 3, 4]]
            pos_int = int(pos)
            chr_ = refine_chr(chr_)
            if transcript_id in ppchange and ppchange[transcript_id] != 'NA':
                logfc = ppchange[transcript_id]
                row = [transcript_id, chr_, pos, strand, logfc]
                strand_info[transcript_id] = strand

                left_boundary = pos_int - region_size - half_bin_size
                right_boundary = pos_int + region_size + half_bin_size
                if left_boundary < 0:
                    continue
                row_tmp = [chr_, left_boundary, right_boundary, transcript_id, strand]
                print('\t'.join(map(str, row_tmp)), file=o)
                # tss_res.append(row)
                tss_with_padding_n += 1
    if tss_with_padding_n == 0:
        logger.error(f'no gene in tss file match with the gene list, please check if the transcript names in the GTF file match with the TSS file')
        return 1
    
    def makewindows(fn_lb, fn_tss_padding, bin_size, chr_map):
        # rewrite the bedtools makewindows, because the original one does not contain the strand information
        fno = f'{pwout}/intermediate/{fn_lb}.split_regions.bed'
        with open(fn_tss_padding) as f, open(fno, 'w') as o:
            for i in f:
                chr_, s, e, ts, strand = i[:-1].split('\t')
                s, e = int(s), int(e)
                bin_sn = 0
                chr_orig = chr_map.get(chr_, chr_)
                for pos in range(s, e, bin_size):
                    bin_sn += 1
                    print(f'{chr_orig}\t{pos}\t{pos + bin_size}\t{ts}_{bin_sn}\t{strand}', file=o)
        # logger.info('sort split region bed 2')
        fntmp = f'{fno}.tmp'
        os.system(f'bedtools sort -i {fno} > {fntmp} && mv {fntmp} {fno}')

        return fno
    
    # build the new tss for the first file, and use it  for the rest
    fn_lb = fls_ctrl[0][0]
    fn_chr_map = f'{pw_bed}/{fn_lb}.chr_map.pkl'
    with open(fn_chr_map, 'rb') as f:
        chr_map = pickle.load(f)
    split_bed_per_file = makewindows(fn_lb, fn_tss_padding, bin_size, chr_map)
    
    # load the normalization factors (nf.txt)
    fn_nf = f'{pwout}/intermediate/nf.txt'
    factors_dict = {}
    with open(fn_nf) as f:
        f.readline() # header
        for i in f:
            sample, factor = i.strip().split('\t')
            factors_dict[sample] = float(factor)
    
    # logger.debug(factors_dict)
    # main part, get the overlap of the input bed files with the split bin bed
    # coverage_by_strand_flag = '-s'  # by_strand, if used, need to update the reference tss file
    coverage_by_strand_flag = '' # current setting
    for condition, fls in {'control': fls_ctrl, 'case': fls_case}.items():
        fn_count_sum = f'{pwout}/intermediate/{condition}.count'
        count = {}   # the sum of norm count of each bin across all samples in this condition
        n_sam_condition = len(fls)
        for fn_lb, fn_bed  in fls:
            norm_factor = factors_dict[fn_lb]
            
            fn_coverage_tmp = f'{pwout}/intermediate/{fn_lb}.coverage_count.tmp'
            # do we need to calculate the coverage by strand??? 

            # if not os.path.exists(fn_coverage_tmp):
            logger.info(f'getting coverage for {fn_lb}')
            s = time.time()
            cmd = f'bedtools coverage -a {split_bed_per_file} -b {fn_bed} -counts {coverage_by_strand_flag} -sorted > {fn_coverage_tmp}'
            # logger.debug(cmd)
            retcode = run_shell(cmd)
            if retcode:
                logger.error(f'failed to run bedtools coverage for {fn_lb}')
                return 1
            # process the coverage file
            # chr1	323932	323952	NR_028322.1_8   10
            dur = time.time() - s
            time_bedtools_coverage += dur
            logger.debug(f'done, time used = {dur:.2f}s')

            with open(fn_coverage_tmp) as f:
                for i in f:
                    _, transcript_id_chunk, ict = i.strip().rsplit('\t', 2)
                    if ict == '0':
                        continue
                    transcript_id, bin_sn = transcript_id_chunk.rsplit('_', 1)
                    # bin_sn = int(bin_sn)
                    ict = int(ict)
                    ict_norm = ict * norm_factor
                    count.setdefault(transcript_id, {}).setdefault(bin_sn, 0)
                    count[transcript_id][bin_sn] += ict_norm


        with open(fn_count_sum, 'w') as o:
            half_bin_number = (bin_number - 1) // 2
            # header = ['Transcript'] + [f'up_{i}' for i in range(half_bin_number, 0, -1)] + ['tss'] + [f'down_{i}' for i in range(1, half_bin_number + 1)]
            header = [f'up_{i}' for i in range(half_bin_number, 0, -1)] + ['tss'] + [f'down_{i}' for i in range(1, half_bin_number + 1)]
            print('\t'.join(header), file=o)
            bin_num_plus = [str(_) for _ in range(1, bin_number + 1)]
            bin_num_minus = [str(_) for _ in range(bin_number, 0, -1)]
            for transcript_id, v1 in count.items():
                strand = strand_info[transcript_id]
                bin_count_list = []
                bin_number_order = bin_num_plus if strand == '+' else bin_num_minus
                
                for bin_sn in bin_number_order:
                    bin_count_list.append(str(round(v1.get(bin_sn, 0) / n_sam_condition, 2))) # average the count across all samples
                tmp = "\t".join(bin_count_list)
                print(f'{transcript_id}\t{tmp}', file=o)
    
    logger.debug(f'total time used for bedtools coverage: {time_bedtools_coverage:.2f}s')

    # heatmap.R  --args file=\"$data_bed\" outdir=\"$inter_dir\" pname=\"$outname\" window=$bin_size region=$region_size
    # each row is [transcript_id, chr_, pos, strand, logfc]
    # read the count table
    df_case = pd.read_csv(f'{pwout}/intermediate/case.count', sep='\t')
    df_ctrl = pd.read_csv(f'{pwout}/intermediate/control.count', sep='\t')
    case_log = np.log2(df_case + 1)
    ctrl_log = np.log2(df_ctrl + 1)
    df_delta = case_log - ctrl_log
    abs_values = np.abs(df_delta.values.flatten())
    cutoff = np.quantile(abs_values, 0.75) # verified, match with R code
    cutoff1 = round(cutoff, 1)
    df_delta = df_delta.clip(lower=-cutoff, upper=cutoff)
    # logger.info(df_delta.head())
    # plot

    # Set up the layout and figure size
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [10, 1]}, figsize=(6, 18))
    plt.subplots_adjust(hspace=0)

    lab = round(region_size / 1000, 1) # label for the heatmap
    # Define color palette
    colors = ["blue", "yellow"]

    # Create a colormap with 5 discrete colors
    my_palette = LinearSegmentedColormap.from_list("custom_palette", colors, N=5)
    
    # Bottom color bar (equivalent to the first image call in R)
    axs[1].imshow(np.arange(1, 6).reshape(-1, 1), cmap=my_palette, aspect='auto')
    axs[1].axis('off')
    axs[1].xaxis.set_ticks_position('none')
    axs[1].set_xticks([0, 2, 4])
    axs[1].set_xticklabels([f"-{cutoff1}", "0", f"+{cutoff1}"], fontsize=12)

    # Main heatmap (equivalent to the second image call in R)
    im = axs[0].imshow(df_delta, cmap=my_palette, vmin=-cutoff, vmax=cutoff, aspect='auto')
    axs[0].axis('off')
    axs[0].set_xticks([0, df_delta.shape[1] // 2, df_delta.shape[1] - 1])
    axs[0].set_xticklabels([f"-{lab}K", "0", f"+{lab}K"], fontsize=12)

    # Save the figure as PDF
    plt.savefig(fn_out_pdf, bbox_inches='tight', dpi=300)
    plt.close()
    
    # os.system(f'rm {pwout}/intermediate/*.tmp')


def groseq_fisherexact_pausing_for_gene(c1, l1, c2, l2):
    """
    Perform Fisher's exact test to calculate the p-value for the pausing index of each gene. x is a row of the pausing count dataset
    per run, is about 1ms, so it is quite slow., 1-2ms per run
    user new fisher package, the speed is 100x faster (6us per call)
    must make sure the input are all non-zero
    the pvalue is slightly different from R fisher.test(x, alternative='two.sided') esp. for small pvalues
    """
    if not all([c1, l1, c2, l2]):
        return np.nan
    # null hypothesis: the pause read count is uniformly distributed in the gene body and the promoter region
    expectC1 = round(l1 * (c1 + c2) / (l1 + l2))
    expectC2 = c1 + c2 - expectC1
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
        return np.nan
    
    # default = two-sided
    return fisher.pvalue(c1, c2, c3, c4).two_tail


def get_FDR_per_sample(fn_lb, fn_peak, fn_fdr, pause_index_str):
    """
    get the pvalue and FDR based on the peak file, retain only 1 transcript per gene, (keep the transcript with the highest pp count)
    """
    logger.debug('loading raw peak table')
    # df_peak = pd.read_csv(fn_peak, sep='\t')
    df_peak_gene = pd.read_csv(fn_peak, sep='\t')
    
    # below code is for collapse genes
    # but actually, do not collapse them. this will lead to lots of problems for combining FDR files from different samples. e.g. building the pindex.txt
    # Transcript	Gene	ppc	ppm	ppd	pps	gbc	gbm	gbd	pauseIndex
    # NM_000014.5	A2M	40	50	0.8	12-:0	43	47522	0.0009	884.13023
    # idx = df_peak.groupby('Gene')['ppc'].idxmax()
    # df_peak_gene = df_peak.loc[idx]
    logger.debug('get pvalue by fisher exact test')
    df_peak_gene['pvalue'] = df_peak_gene.apply(lambda x: groseq_fisherexact_pausing_for_gene(x['ppc'], x['ppm'], x['gbc'], x['gbm']), axis=1)
    
    pvals = df_peak_gene['pvalue'].dropna()
    logger.debug('get FDR')
    fdr = multipletests(pvals, method='fdr_bh')[1]
    df_peak_gene.loc[pvals.index, 'FDR'] = fdr
    logger.debug(f'dump to file {fn_fdr}')
    df_peak_gene.to_csv(fn_fdr, sep='\t', index=False, na_rep='NA')
    
    logger.debug('building pause_index_str')
    col_ts = df_peak_gene.columns[0]
    df_tmp= df_peak_gene[[col_ts, 'pauseIndex', 'pvalue', 'FDR']].fillna('NA').astype(str)
    for row in df_tmp.itertuples(index=False, name=None):
        ts, *values = row
        pause_index_str.setdefault(ts, {})[fn_lb] = values
    return pause_index_str


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
        log2fc = np.log2(data.iloc[:, n_gene_cols + 1] / size_factors[1]) / (data.iloc[:, n_gene_cols] / size_factors[0])
        res_df = data.iloc[:, idx_gene_cols]
        res_df['log2FoldChange'] = log2fc
        return res_df, size_factors

    if size_factors_in is None:  # gb change
        dds.deseq2()
    else:
        # logger.debug(f'custom size factor specified: {size_factors_in}')
        counts = dds.X
        (dds.logmeans, dds.filtered_genes) = deseq2_norm_fit(counts)
        dds.layers["normed_counts"] = counts / size_factors[:, None]
        
        # added below for 0.4.10 version of pydeseq2
        dds.varm["_normed_means"] = dds.layers["normed_counts"].mean(0)


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

    sam_list = [re.sub('^ppc_', '', _) for _ in data.columns[idx_ppc_combined]]
    
    if skip_filtering:
        data_pass = data
        data_drop = pd.DataFrame()
        return col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop
    
    tmp1 = sum([1 if cols_raw[i].startswith('gbc_') else 0 for i in idx_gbc_combined])
    tmp2 = sum([1 if cols_raw[i].startswith('gbd_') else 0 for i in idx_gbd_combined])
    if tmp1 != n_sam or tmp2 != n_sam:
        logger.error('the number of samples in gbc or gbd is not correct')
        sys.exit(1)
        return None
    
    
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
    
    data_pass.iloc[:, idx_gene_cols].to_csv(f'{out_dir}/intermediate/active_gene.txt', sep='\t', index=False, header=False, na_rep='NA')
    
    data_pass_pp = data_pass.iloc[:, idx_gene_cols + idx_ppc_combined] # datapp in R
    data_drop_pp = data_drop.iloc[:, idx_gene_cols + idx_ppc_combined] # data_pp in R
    
    data_pass_gb = data_pass.iloc[:, idx_gene_cols + idx_gbc_combined] # datagb in R
    data_drop_gb = data_drop.iloc[:, idx_gene_cols + idx_gbc_combined] # data_gb in R
    
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
        logger.info('running DESeq2')
        with HiddenPrints():
            res_df, size_factors = run_deseq2(n_gene_cols, data_pass_gb, metadata, ref_level, size_factors_in=size_factors)
            res_df_full = pd.concat([res_df, data_drop_gb.iloc[:, idx_gene_cols]])
            
            res_df_full.to_csv(f'{out_dir}/known_gene/gb_change.txt', sep='\t', index=False, na_rep='NA')
            
            # pp change
            try:
                res_df, _ = run_deseq2(n_gene_cols, data_pass_pp, metadata, ref_level=ref_level, size_factors_in=size_factors)
            except:
                fn_deseq2_input = f'{outdir}/intermediate/deseq2.pkl'
                
                logger.error(f'fail to run deseq2, dumping the arguments to {fn_deseq2_input}')
                with open(fn_deseq2_input, 'wb') as o:
                    pickle.dump([n_gene_cols, data_pass_pp, metadata, ref_level, size_factors], o)
                
                sys.exit(1)
                
            res_df_full = pd.concat([res_df, data_drop_pp.iloc[:, idx_gene_cols]])
            res_df_full.to_csv(f'{out_dir}/known_gene/pp_change.txt', sep='\t', index=False, na_rep='NA')
    elif rep1 == 1:
        logger.debug('single contrl sample only')
        size_factors = 1 # size factor is 1 for the only sample
    elif rep2 == 0:
        # control only, multiple ctrol samples
        _, size_factors = run_deseq2(n_gene_cols, data_pass_gb, metadata, ref_level, size_factors_in=size_factors, size_factor_only=True)

    norm_factors = 1 / size_factors
    # get the normalized data in other function, because we need the chr start end strand information, and the data here is already processed
    # save normalization factors
    nf = pd.DataFrame({'sample': sam_list, 'nfactor': norm_factors})
    nf.to_csv(f'{out_dir}/intermediate/nf.txt', sep='\t', index=False, na_rep='NA')
    return size_factors, sam_list
    

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
    data.to_csv(f'{out_dir}/known_gene/alternative_isoform{norm_flag}.txt', index=False, sep='\t', na_rep='NA')

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
    fn_norm = f'{out_dir}/known_gene/normalized_pp_gb.txt'

    
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
    size_factors, sam_list = change_pp_gb_with_case(n_gene_cols, rep1, rep2, data, out_dir, window_size, factor_flag, factor1, factor2)

    # get the normalized data
    n_extra_cols = 4 # chr, start, end, strand
    n_prev_cols = n_gene_cols + n_extra_cols
    col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop = filter_pp_gb(data_raw, n_prev_cols, rep1, rep2, skip_filtering=True)
    
    norm_factors = 1 / size_factors
    ppc_norm = data_raw.iloc[:, idx_ppc_combined] * norm_factors
    ppd_norm = ppc_norm / window_size
    ppd_norm.columns = [f'ppd_{_}' for _ in sam_list]
    
    gbc_norm = data_raw.iloc[:, idx_gbc_combined] * norm_factors
    gbd_norm = data_raw.iloc[:, idx_gbd_combined] * norm_factors
    data_normalized = pd.concat([data_raw.iloc[:, :n_prev_cols], ppc_norm, ppd_norm, gbc_norm, gbd_norm], axis=1)
    data_normalized.to_csv(fn_norm, sep='\t', index=False, na_rep='NA')


def cmhtest(row, rep1, rep2):
    """
    performs the Mantel-Haenszel test on a certain gene
    the row is like ppc_sam1, ppc_sam2, ppc_sam3, ppc_sam4, gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2, ...
    """
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
            if pp1 + pp2 + gb1 + gb2 > 0:
                arr.append([[pp2, gb2], [pp1, gb1]])
    arr = np.array(arr).T
    cmt = contingency_tables.StratifiedTable(tables=arr)
    odds_ratio = cmt.oddsratio_pooled
    test_res = cmt.test_null_odds(correction=True)
    pvalue = test_res.pvalue
    # statistic = test_res.statistic
    # print(f'odds_ratio = {odds_ratio}, pvalue = {pvalue}, statistic = {statistic}') 
    # return log2(odd_ratio), pvalue
    return np.log2(odds_ratio) if odds_ratio != 0 else 'NA', pvalue


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

def add_FDR_col(df, col_pvalue, col_FDR='FDR'):
    pvals = df[col_pvalue].replace('NA', np.nan).dropna()
    fdr = multipletests(pvals, method='fdr_bh')[1]
    df.loc[pvals.index, 'FDR'] = fdr
    return df


def change_pindex(fno_prefix, n_gene_cols, fn, out_dir, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0):
    """
    fno_prefix, default is empty for known genes, longeRNA- for longeRNA
    adapted from Rscript change_pindex.R
    input = intermediate/ count_pp_gb.txt or longeRNA-count_pp_gb.txt  Transcript, Gene, chr, start, end, strand, (ppc_sam1 - ppc_sum2 ...), (gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ...)
    
    The transcripts were not collapsed to gene level
    """
    # for longeRNA, the count_pp_gb header is like
    # long-eRNA	ppc_ ...., so there is no Transcript, gene, chr, start, end, and strand
    # compared with previous results, matches
    if fno_prefix not in {'', 'longeRNA-'}:
        logger.error(f'invalid file prefix for pindex_change, should be empty or "longeRNA-", input="{fno_prefix}"')
        return 1
    if not rep2:
        logger.warning(f'Only control samples were specified, skip performing change_pindex...')
        return 
    is_erna = True if fno_prefix == 'longeRNA-' else False
    data = pd.read_csv(fn, sep='\t')
    n_sam = rep1 + rep2
    if is_erna:
        # the data does not contains the chr, start, end, strand columns, and there is no Gene column
        # no filtering needed
        ana_type = 'eRNA'
        n_cols_prev = n_gene_cols = 1
    else:
        ana_type = 'Known genes'
        n_cols_prev = n_gene_cols + 4
        cols_keep = list(range(n_gene_cols)) + list(range(n_cols_prev, len(data.columns))) # drop the chr, start, end, strand columns
        data = data.iloc[:, cols_keep]

    col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop = filter_pp_gb(data, n_gene_cols, rep1, rep2, skip_filtering=is_erna)
    
    if n_sam == 2:
        # use fisher exact test
        idx_ppc1 = col_idx['ppc']['cond1'][0]
        idx_ppc2 = col_idx['ppc']['cond2'][0]
        idx_gbc1 = col_idx['gbc']['cond1'][0]
        idx_gbc2 = col_idx['gbc']['cond2'][0]
        data_out = data_pass.iloc[:, idx_gene_cols].copy()
        data_out['pvalue'] = data_pass.apply(lambda x: get_pvalue_2_sample(x, idx_ppc1, idx_ppc2, idx_gbc1, idx_gbc1), axis=1)
        data_out = add_FDR_col(data_out, 'pvalue')
        # logger.info(data_out.head())

        
        # odds_ratio = (ppc2/ppc1) / (gbc2/gbc1) = ppc2 * gbc1 / (ppc1 * gbc2)
        def get_odds_ratio(row):
            ppc1 = row[idx_ppc1]
            ppc2 = row[idx_ppc2]
            gbc1 = row[idx_gbc1]
            gbc2 = row[idx_gbc2]
            tmp = ppc1 * gbc2
            if all([ppc1, ppc2, gbc1, gbc2, tmp]):
                return ppc2 * gbc1 / (tmp)
            return np.nan

        try:
            # log2fc = data_pass.apply(lambda x: np.log2(x[idx_ppc2] * x[idx_gbc1] / (x[idx_ppc1] * x[idx_gbc2])), axis=1)
            log2fc = data_pass.apply(lambda x: np.log2(get_odds_ratio(x)), axis=1)
        except:
            e = traceback.format_exc()
            logger.warning(e)
            x = list(data_pass.iloc[0])
            logger.warning(x)
            logger.info(get_odds_ratio(x))
            logger.warning(f'{x[idx_ppc2]=}, {x[idx_gbc1]=}, {x[idx_ppc1]=}, {x[idx_gbc2]=}')
            sys.exit(1)
        
        data_out['log2fc'] = log2fc

    # run the cmhtest, the data should have 3 dimensions
    else:
        # logger.info(data_pass.columns)
        # logger.info(data_pass.head().values)
        # cmhtest(row, rep1, rep2)  # return = log2(odds_ratio), pvalue
        data_out = data_pass.iloc[:, idx_gene_cols].copy()
        data_out[['log2fc', 'pvalue']] = data_pass.iloc[:, n_gene_cols:].apply(lambda x: cmhtest(list(x), rep1, rep2), axis=1, result_type='expand')
        data_out = add_FDR_col(data_out, 'pvalue')
        data_out = data_out.sort_values('FDR')
        
    data_out_full = pd.concat([data_out, data_drop.iloc[:, idx_gene_cols]]).fillna('NA')
    data_out_full.to_csv(f'{out_dir}/known_gene/{fno_prefix}pindex_change.txt', sep='\t', index=False, na_rep='NA')
    # logger.info(f'pindex_change done : {ana_type}')
    # logger.info(data_out.head())


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
    
def add_value_to_gtf(gene_info, pro_up, pro_down, gb_down_distance):
    strand = gene_info['strand']
    gene_raw_s, gene_raw_e = gene_info['start'], gene_info['end']
    if strand == '+':
        pp_start = gene_raw_s - pro_up
        pp_end = gene_raw_s + pro_down - 1
        gb_start = gene_raw_s + gb_down_distance
        gb_end = gene_raw_e
        strand_idx = 0
    else:
        pp_start = gene_raw_e - (pro_down - 1)
        pp_end = gene_raw_e + pro_up
        gb_start = gene_raw_s
        gb_end = gene_raw_e - gb_down_distance
        strand_idx = 1
    
    # modify here
    # skip the get gene_seq step
    gene_seq = None
    gb_seq_N = 0

    # gene_seq = get_seqence_from_fa(fa_idx, fh_fa, chr_, gene_raw_s, gene_raw_e)# already upper case
    # if gene_seq.count('N') == 0:
    #     gene_seq = None  # no need to check the mappable sites
    #     gb_seq_N = 0
    # else:
    #     gb_seq_N = gene_seq[gb_start - gene_raw_s: gb_end - gene_raw_s + 1].count('N')
    
    gb_len_mappable = gb_end - gb_start + 1 - gb_seq_N

    new_info = dict(zip(
        ['pp_start', 'pp_end', 'gb_start', 'gb_end', 'strand_idx', 'gb_len_mappable', 'gene_seq'], 
        [ pp_start,   pp_end,   gb_start ,  gb_end ,  strand_idx ,  gb_len_mappable,   gene_seq ]
        ))
    gene_info.update(new_info)
    return gene_info

def build_tss(gtf_info, fn_tss):
    # create the TSS file
    # chr1	11874	11874	NR_046018.2	+
    logger.info('creating TSS file')
    with open(fn_tss, 'w') as f:
        for k, v in gtf_info.items():
            itss = v['start'] if v['strand'] == '+' else v['end']
            f.write(f'{v["chr"]}\t{itss}\t{itss}\t{k}\t{v["strand"]}\n')

def process_gtf(fn_gtf, pwout):
    """
    process the gtf file, get the gene start and end position, and other info
    gtf position is 1-idx, full closed
    """

    err = {'no_transcript_id': 0, 'no_gene_name': 0, 'invalid_line_format': 0, 'invalid_strand': 0}
    fn_gtf = os.path.abspath(fn_gtf)
    fn_gtf_lb = os.path.basename(fn_gtf).replace('.gz', '').replace('.gtf', '')

    fn_tss = f'{pwout}/intermediate/{fn_gtf_lb}.tss.txt'

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
        return None, fn_tss, err
    
    fn_gtf_pkl = f'{gtf_out_dir}/{fn_gtf_lb}.gtf_info.pkl'
    if os.path.exists(fn_gtf_pkl):
        logger.debug('loading gtf from pickle')
        
        with open(fn_gtf_pkl, 'rb') as f:
            gtf_info = pickle.load(f)
        if not os.path.exists(fn_tss):
            build_tss(gtf_info, fn_tss)
        return gtf_info, fn_tss, err
    
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
    res_raw = {}  # k1 = gene name, k2 = transcript_id, v = {'chr': '', 'strand': '', 'gene_id': '', 'gene_name': '', 'start': 0, 'end': 0}
    
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
            
            # {'chr': '', 'strand': '', 'gene_name': '', 'start': 0, 'end': 0}
            res_raw.setdefault(gene_name, {}).setdefault(transcript_id, {'chr': chr_, 'strand': strand, 'gene_name': gene_name, 'start': start, 'end': end})
            
            ires = res_raw[gene_name][transcript_id]
            
            # deal with the case when the same transcript-ID with different chr or strand
            if chr_ != ires['chr'] or strand != ires['strand']:
                transcript_id_new = f'{transcript_id}_{chr_}_{strand}'
                res_raw.setdefault(gene_name, {}).setdefault(transcript_id_new, {'chr': chr_, 'strand': strand, 'gene_name': gene_name, 'start': start, 'end': end})
                ires = res_raw[gene_name][transcript_id_new]

            if start < ires['start']:
                ires['start'] = start
            if end > ires['end']:
                ires['end'] = end

    # merge transcripts with same start and end position
    res = {}
    n_merged = 0
    for gn, v1 in res_raw.items():
        tmp = {}  # key = unique_id
        for transcript_id, v2 in v1.items():
            unique_id = f'{v2["chr"]}_{v2["start"]}_{v2["end"]}_{v2["strand"]}'
            tmp.setdefault(unique_id, []).append(transcript_id) # these transcripts have the same start and end position
        for transcript_list in tmp.values():
            if len(transcript_list) > 1:
                n_merged += 1
            transcript_id_new = ';'.join(transcript_list)
            res[transcript_id_new] = v1[transcript_list[0]]

    
    logger.warning(f'merged {n_merged} transcripts with same start and end position')
    
    with open(fn_gtf_pkl, 'wb') as o:
        pickle.dump(res, o)

    build_tss(res, fn_tss)
    err_total = sum(err.values())
    if err_total:
        logger.info(f'error in parsing gtf file: {err}')
    return res, fn_tss, err



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


def run_shell(cmd, echo=False):
    """
    run shell command and check the return code and stdout stderr
    """
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    retcode = p.returncode
    if retcode and echo:
        logger.error(f"Error running command: {cmd}")
        logger.error(f"stdout: {stdout.decode()}")
        logger.error(f"stderr: {stderr.decode()}")
    return retcode

def get_lb(fn):
    return fn.rsplit('/', 1)[-1].rsplit('.', 1)[0]



def pre_count_for_bed(fn_lb, fn_out_bed, pw_bed, bin_size, reuse=True):
    """
    fh_bed can be the real file handle or a subprocess.Popen.stdout object, so do not use seek or tell
    process the bed file, and get the count of the read end regarding strand.
    will export 2 files, 
    1. dict, key = chr, k2 = strand, k3 = chunk_number. value = count of reads with read end in this chunk.  the chunk_number is got from position // bin_size  (default chunksize = 10bp)
    2. bed file, record the read end chr, read_end position, strand and the count of reads with this read end position
    
    bin_size = 200 # compared bin size of 10, 20, 50, 100 and 200 and 500.  200 is the best, maybe because when getting the gb count, the size fit the distribution
    return chr_map is used to make the chr pattern in the gtf file and input bed file consistant, esp. when running bedtools

    """
    fn_count_bin = f'{pw_bed}/{fn_lb}.count.bin_of_{bin_size}.pkl'
    fn_count_per_base = f'{pw_bed}/{fn_lb}.count.per_base.pkl'
    fn_chr_map = f'{pw_bed}/{fn_lb}.chr_map.pkl'
    fn_n_lines = f'{pw_bed}/{fn_lb}.line_count.txt'
    
    if reuse and os.path.exists(fn_count_per_base) and os.path.getsize(fn_count_per_base) > 10:
        logger.info(f'Loading pre-counting data...')
        with open(fn_count_per_base, 'rb') as f:
            count_per_base = pickle.load(f)
        with open(fn_count_bin, 'rb') as f:
            count_bin = pickle.load(f)
        logger.info('Loaded.')
        return count_per_base, count_bin

    count_bin = {'bin_size': bin_size}
    count_per_base = {}
    logger.info(f'Pre-counting for {fn_lb}')
    n_lines = 0
    with gzip.open(fn_out_bed, 'rt') if fn_out_bed.endswith('.gz') else open(fn_out_bed, 'r') as fh_bed:
        for i in fh_bed:
            # chr1	10511	10569	A00758:60:HNYK7DSXX:4:2461:27181:23171	32	-
            chr_, s, e, _, _, strand = i[:-1].split('\t')
            strand_idx, read_end = (0, int(e)) if strand == '+' else (1, int(s) + 1)

            chunk = read_end // bin_size
            count_bin.setdefault(chr_, {}).setdefault(chunk, [0, 0])
            count_bin[chr_][chunk][strand_idx] += 1

            count_per_base.setdefault(chr_, {}).setdefault(read_end, [0, 0])  # for each read end, the count of + and - strand
            count_per_base[chr_][read_end][strand_idx] += 1
            n_lines += 1

        # refine the chr_
        chr_map = {refine_chr(k): k for k in count_bin}
        count_per_base = {refine_chr(k): v for k, v in count_per_base.items()}
        count_bin = {refine_chr(k): v for k, v in count_bin.items()}
    logger.info('done.')
    # dump the pickle
    with open(fn_count_per_base, 'wb') as f:
        pickle.dump(count_per_base, f)
    with open(fn_count_bin, 'wb') as f:
        pickle.dump(count_bin, f)
        
    with open(fn_chr_map, 'wb') as f:
        pickle.dump(chr_map, f)
    logger.debug('dump to pickle done.')
    
    with open(fn_n_lines, 'w') as f:
        f.write(f'{n_lines}\n')
    
    
    return count_per_base, count_bin


def process_input(pwout_raw, fls):
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
        gz_suffix = '.gz' if fn.endswith('.gz') else ''
        fn_out_bed = f'{pwout_raw}/bed/{fn_lb}.sorted.bed'
        fn_out_bed_gz = f'{pwout_raw}/bed/{fn_lb}.sorted.bed{gz_suffix}'

        if os.path.exists(fn_out_bed) and os.path.getsize(fn_out_bed) > 10:
            res.append([fn_lb, fn_out_bed])
            continue
        elif gz_suffix and os.path.exists(fn_out_bed_gz) and os.path.getsize(fn_out_bed_gz) > 10:
            res.append([fn_lb, fn_out_bed_gz])
            continue
        
        fn_for_check = re.sub(r'\.gz$', '', fn)
        if fn_for_check.endswith('.bam'):
            logger.info(f'Converting {fn}')
            # still need to do the sorting, because for the following steps using bedtools coverage, the sorted bed will be more memory efficient
            # bedtools sort is faster than linux sort
            cmd = f"""bedtools bamtobed -i {fn} |bedtools sort -i - > {fn_out_bed}"""
            os.system(cmd)
            ires = [fn_lb, fn_out_bed]
        elif fn_for_check.endswith('.bed'):
            # check if sorted
            if 'sorted' in fn:
                is_sorted = 1
            else:
                is_sorted = check_is_sorted(fn)
                
            if is_sorted:
                fn_abs = os.path.realpath(fn)
                os.system(f'ln -sf {fn_abs} {fn_out_bed_gz}')
                ires = [fn_lb, fn_out_bed_gz]
            else:
                logger.warning(f'input bed is not sorted, now sorting...')
                # bedtools can handle gzip format
                cmd = f'bedtools sort -i {fn} > {fn_out_bed}'
                os.system(cmd)
                ires = [fn_lb, fn_out_bed]
        else:
            logger.error(f"Input file '{fn}' should be in bed or bam format")
            err = 1
            continue
        res.append(ires)
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

