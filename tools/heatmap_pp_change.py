#! /data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12

import sys, os, re, pickle, time
import pandas as pd
import numpy as np
import traceback

def getlogger(fn_log=None, logger_name=None, nocolor=False):
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
            'black': '\u001b[30;1m',
            'red': '\u001b[31;1m',
            'r': '\u001b[31;1m',
            'bold_red': '\u001b[31;1m',
            'rb': '\u001b[31;1m',
            'green': '\u001b[32;1m',
            'g': '\u001b[32;1m',
            'gb': '\u001b[32;1m',
            'yellow': '\u001b[33;1m',
            'blue': '\u001b[34;1m',
            'b': '\u001b[34;1m',
            'purple': '\u001b[35;1m',
            'p': '\u001b[35;1m',
            'grey': '\u001b[38;1m',
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
        console.setLevel('INFO')
        console.name = 'console'
        logger.addHandler(console)

    if fn_log and 'file' not in handler_names:
        fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
        fh_file.setLevel('DEBUG')
        fh_file.setFormatter(CustomFormatter())
        fh_file.name = 'file'
        logger.addHandler(fh_file)
    return logger
logger = getlogger()

def refine_chr(chr_):
    chr_ = chr_.strip().lower()
    for _ in ['chromosome', 'chrom', 'chr']:
        chr_ = chr_.replace(_, '')
    return {'23': 'x', '24': 'y', '25': 'mt', 'm': 'mt'}.get(chr_) or chr_


def draw_heatmap_pp_change(n_gene_cols, pwout, pw_bed, fls_ctrl, fls_case, fn_tss, by_strand=False, region_size=5000, bin_size=200, outname='heatmap_pp_change'):
    """
    draw heatmap for pp change
    fls_ctrl, fls_case, foramt is like [fn_lb, fn_bed]
    region_size, upstream and downstream distance relative to TSS for plotting PRO-seq signal (bp, default: 5000),should can be divided by bin size
    bin_size (bp, default: 200),should can be divided by 2
    """
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap
    fls_case = fls_case or []
    
    ppchange = {}
    fn_pp_change = f'{pwout}/known_gene/pp_change.txt'
    fn_out_pdf = f'heatmap.pdf'
    # fn_data_bed = f'{pwout}/intermediate/data_bed.tmp'
    fn_tss_padding = f'infile_bed.tmp' # expand the fn_tss with upstream and downstream region by region_size
    time_bedtools_coverage = 0

    if bin_size % 2 or region_size % bin_size:
        logger.error(f'bin_size should be even and region_size should be divisible by bin_size')
        return 1

    bin_number = int(region_size * 2 / bin_size + 1)  # *2 because of the upstream and downstream

    def to_number(s):
        if not s or s.upper() == 'NA':
            return 'NA'
        if 'inf' in s:
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
        fno = f'{fn_lb}.split_regions.bed'
        with open(fn_tss_padding) as f, open(fno, 'w') as o:
            for i in f:
                chr_, s, e, ts, strand = i[:-1].split('\t')
                s, e = int(s), int(e)
                bin_sn = 0
                chr_orig = chr_map.get(chr_, chr_)
                for pos in range(s, e, bin_size):
                    bin_sn += 1
                    print(f'{chr_orig}\t{pos}\t{pos + bin_size}\t{ts}_{bin_sn}\t-1\t{strand}', file=o)
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
    
    if by_strand:
        coverage_by_strand_flag = '-s'  # by_strand, if used, need to update the reference tss file
        logger.warning(f'currently getting coverage by strand')
    else:
        coverage_by_strand_flag = '' # current setting
    for condition, fls in {'control': fls_ctrl, 'case': fls_case}.items():
        fn_count_sum = f'{condition}.count'
        count = {}   # the sum of norm count of each bin across all samples in this condition
        n_sam_condition = len(fls)
        for fn_lb, fn_bed  in fls:
            norm_factor = factors_dict[fn_lb]
            
            fn_coverage_tmp = f'{fn_lb}.coverage_count.tmp'
            # do we need to calculate the coverage by strand??? 

            # if not os.path.exists(fn_coverage_tmp):
            if not os.path.exists(fn_coverage_tmp):
                logger.info(f'getting coverage for {fn_lb}')
                s = time.time()
                cmd = f'bedtools coverage -a {split_bed_per_file} -b {fn_bed} -counts {coverage_by_strand_flag} -sorted > {fn_coverage_tmp}'
                # logger.debug(cmd)
                retcode = os.system(cmd)
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
                    try:
                        # chr13	51846175	51846375	NM_011817.2_24	+	9
                        _, transcript_id_chunk, _, strand, ict = i.strip().rsplit('\t', 4)
                        if ict == '0':
                            continue
                        transcript_id, bin_sn = transcript_id_chunk.rsplit('_', 1)
                        # bin_sn = int(bin_sn)
                        ict = int(ict)
                        ict_norm = ict * norm_factor
                        count.setdefault(transcript_id, {}).setdefault(bin_sn, 0)
                        count[transcript_id][bin_sn] += ict_norm
                    except:
                        e = traceback.format_exc()
                        logger.error(e)
                        logger.error(f'error when processing {fn_coverage_tmp}, line = \n{i}\n')
                        sys.exit(1)

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
    df_case = pd.read_csv(f'case.count', sep='\t')
    df_ctrl = pd.read_csv(f'control.count', sep='\t')
    case_log = np.log2(df_case + 1)
    ctrl_log = np.log2(df_ctrl + 1)
    df_delta = case_log - ctrl_log
    abs_values = np.abs(df_delta.values.flatten())
    cutoff = np.quantile(abs_values, 0.75) # verified, match with R code
    cutoff1 = round(cutoff, 1)
    df_delta = df_delta.clip(lower=-cutoff, upper=cutoff)
    # logger.info(df_delta.head())
    # plot
    n_color_block = 5
    with PdfPages(fn_out_pdf) as pdf:
    # Set up the layout and figure size
        plt.figure(figsize=(6, 18))
        grid = plt.GridSpec(2, 1, height_ratios=[14, 1])

        lab = round(region_size / 1000, 1) # label for the heatmap
        # Define color palette
        colors = ["blue", "yellow"]

        # Create a colormap with 5 discrete colors
        my_palette = LinearSegmentedColormap.from_list("custom_palette", colors, N=n_color_block)
        n_cols = df_delta.shape[1]
        
        ax1 = plt.subplot(grid[0, 0])
        ax1.imshow(df_delta, cmap=my_palette, vmin=-cutoff, vmax=cutoff, aspect='auto')
        ax1.set_xticks([0, n_cols // 2, n_cols - 1])
        ax1.set_xticklabels([f"-{lab}K", "0", f"+{lab}K"], fontsize=12)
        ax1.set_yticks([])
        ax1.set_xlabel("")
        ax1.set_ylabel("")
        
        
        ax2 = plt.subplot(grid[1, 0])
        ax2.imshow(np.arange(1, n_color_block + 1).reshape(1, -1), cmap=my_palette, aspect='auto')
        ax2.set_xticks([0, 2, 4])
        ax2.set_xticklabels([f'-{cutoff1}', '0', f'+{cutoff1}'])
        ax2.set_yticks([])
        ax2.set_xlabel("")
        ax2.set_ylabel("")
        
        
        pdf.savefig(bbox_inches='tight')
        plt.close()

def get_pwbed(pwout):
    # /nobackup/h_vangard_1/wangj52/eNRSA/test-Matt-IFN-R1/INF_Cre_30min_vs_INF_EV_30min
    tmp = f'{pwout}/bed'
    if os.path.exists(tmp):
        return tmp
    pwroot = os.path.dirname(pwout)
    pw_bed = f'{pwroot}/bed'
    if not os.path.exists(pw_bed):
        logger.error(f'pw_bed not found: {pw_bed}')
        sys.exit(1)
    return pw_bed

def get_tss(pwout):
    # RefSeq-mm10.tss.txt
    from glob import glob
    tmp = glob(f'{pwout}/intermediate/*.tss.txt')
    if len(tmp) == 0:
        logger.error(f'TSS file not found under {pwout}/intermediate')
        sys.exit(1)
    return tmp[0]

def get_fls(pwout, pw_bed, n_ctrl):
    fn_nf = f'{pwout}/intermediate/nf.txt'
    if not os.path.exists(fn_nf):
        logger.error(f'nf.txt not exist : {fn_nf}')
        sys.exit(1)

    with open(fn_nf) as f:
        f.readline()
        sam_list = [_.split('\t')[0] for _ in f]
    
    lbs_ctrl = sam_list[:n_ctrl]
    lbs_case = sam_list[n_ctrl:]
    
    if len(lbs_case) == 0:
        logger.error(f'no case samples found (total sam number = n_ctrl = {n_ctrl}), not applicable for heatmap')
        sys.exit(1)
    
    bed_map = {}
    from glob import glob
    bed_list = glob(f'{pw_bed}/*.bed')
    for fn in bed_list:
        lb = os.path.basename(fn).replace('.bed', '')
        lb1 = lb.replace('.sorted', '')
        if lb not in bed_map or 'sorted' not in bed_map[lb]:
            bed_map[lb] = fn
        if lb1 != lb and (lb1 not in bed_map or 'sorted' not in bed_map[lb1]):
            bed_map[lb1] = fn
    fls_ctrl = []
    fls_case = []
    err = []
    
    for lbs, out_l in [[lbs_ctrl, fls_ctrl], [lbs_case, fls_case]]:
        for lb in lbs:
            if lb not in bed_map:
                err.append(lb)
                continue
            out_l.append([lb, bed_map[lb]])
    
    if err:
        logger.error(f'the following bed files are not found in {pw_bed}:')
        print('\t' + '\n\t'.join(sorted(err)))
        sys.exit(1)
    return fls_ctrl, fls_case

if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('pwout', help="""project output folder, the one with known_gene and intermediate subfolders""")
    ps.add_argument('n_ctrl', help="""number of controls""", type=int)
    ps.add_argument('-strand', '-s', help="""coverage by strand""",action='store_true')
    args = ps.parse_args()
    
    n_gene_cols = 2
    by_strand = args.strand
    n_ctrl = args.n_ctrl
    pwout = args.pwout
    pw_bed = get_pwbed(pwout)
    fn_tss = get_tss(pwout)
    fls_ctrl, fls_case = get_fls(pwout, pw_bed, n_ctrl)
    draw_heatmap_pp_change(n_gene_cols, pwout, pw_bed, fls_ctrl, fls_case, fn_tss, by_strand)
