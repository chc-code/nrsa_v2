#! /data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12
"""
process PROseq data
"""
import time
s = time.time()
import sys, os, re
import pickle
import json
from types import SimpleNamespace


def getargs():
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('-design_table', '-design',  help="""Optional, desgin table in tsv format. 2 sections, first section is the sample information, 2 columns. col1=bed/bam full file path, col2 = group name. Second section is the comparison information, 2 columns, col1 should start with @@ used to identify this line as comparison definition. col1=group name used as case, col2 = group name used as control. e.g. @@case\\tcontrol. If only need to process a single group, use @@null\\tgroup_name""", nargs='?')
    ps.add_argument('-in1', help="""required, read alignment files in bed (6 columns) or bam format for condition1, separated by space""", nargs='+')
    ps.add_argument('-in2', help="""read alignment files in bed (6 columns) or bam format for condition2, separated by space""", nargs='*')
    ps.add_argument('-pwout', '-output', help="""required, output/work directory""", required=True)
    ps.add_argument('-organism', '-org',  help="""required, define the genome. Valid is hg19, hg38, mm10, dm3, dm6, ce10, or danRer10. default: hg19""", required=True)

    ps.add_argument('-gtf', help="""user specified GTF file, if not specified, will use the default GTF file for the organism""")
    ps.add_argument('-fa', help="""Full path for the fasta file for the genome. If not specified, will search under the fa folder under the package directory. e.g. organism is hg19, then the default fasta file should be fa/hg19.fa""")
    ps.add_argument('-f1', help="""normalization factors for samples of condition1, separated by space. If no normalization factor is specified, we will use DEseq2 default method to normalize""", nargs='*', type=float)
    ps.add_argument('-f2', help="""normalization factors for samples of condition2, same as -f1""", nargs='*', type=float)
    ps.add_argument('-up', '-pp_up', help="""define the upstream of TSS as promoter (bp, default: 500)""", type=int, default=500)
    ps.add_argument('-down', '-pp_down', help="""define the downstream of TSS as promoter (bp, default: 500)""", type=int, default=500)
    ps.add_argument('-gb', '-gb_start', help="""define the start of gene body density calculation (bp, default: 1000)""", type=int, default=1000)
    ps.add_argument('-min_gene_len', help="""define the minimum length of a gene to perform the analysis (bp, default: 1000). Genes with length less than this will be ignored""", type=int, default=1000)
    ps.add_argument('-window', '-window_size', help="""define the window size (bp, default: 50)""", type=int, default=50)
    ps.add_argument('-step', '-step_size', help="""define the step size (bp, default: 5)""", type=int, default=5)
    ps.add_argument('-bench', '-test', help="""bench mode, only process the first 10 genes for testing purpose""", action='store_true')
    ps.add_argument('-force', help="""ignore the previous pre-count result for the bed files""", action='store_true')
    # ps.add_argument('-overwrite', help="""if the mapped reads count file already exists, overwrite it, default is reuse the old ones""", action='store_true')
    ps.add_argument('-verbose', '-v', help="""verbose mode, print more information""", action='store_true')
    # ps.add_argument('-testfunc', help="""test the new functions,debug mode""", action='store_true')
    args = ps.parse_args()
    
    return args

args = getargs()


# print('after get args')
# print(args)

from utils import check_dependency, run_shell, process_input, get_seqence_from_fa, build_idx_for_fa, get_ref, process_gtf, get_peak,  change_pp_gb, change_pindex, draw_box_plot, draw_heatmap_pindex, draw_heatmap_pp_change, get_alternative_isoform_across_conditions, get_FDR_per_sample, pre_count_for_bed, add_value_to_gtf

sys.dont_write_bytecode = True

time_cost = {}
now = time.time()

s = now

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


class Analysis:
    def __init__(self, args, is_long_eRNA=False):
        
        self.bin_size = 200 # the bin size for building the pre-count dict of each bed file. compared bin size of 10, 20, 50, 100 and 200 and 500.  200 is the best, maybe because when getting the gb count, the size fit the distribution
        
        # folders
        self.bin_dir = os.path.dirname(os.path.realpath(__file__))
        self.out_dir = args.pwout or os.getcwd()
        self.inter_dir = os.path.join(self.out_dir, 'intermediate')
        self.known_gene_dir = os.path.join(self.out_dir, 'known_gene')
        self.longerna = is_long_eRNA # toggle if this is long eRNA
        # self.overwrite_intermediate = args.overwrite
        
        for lb, d in zip(['Output', 'Intermediate', 'Known Genes', 'bed file folder'], [self.out_dir, self.inter_dir, self.known_gene_dir, f'{self.inter_dir}/bed']):
            if not os.path.exists(d):
                logger.debug(f'Making {lb} directory: {d}')
                os.makedirs(d)

        self.status = 0

        self.organism = args.organism
        if self.organism not in {'hg19', 'hg38', 'mm10', 'dm3', 'dm6', 'ce10', 'danRer10'}:
            logger.error(f"Invalid organism provided: {self.organism}")
            self.status = 1

        # input files
        in1 = process_input(self.out_dir, args.in1) # return = fn_lb, fn_bed
        if in1 is None:
            logger.error("Invalid input files provided for condition1")
            self.status = 1
        in2 = process_input(self.out_dir, args.in2)
        self.control_bed = in1
        self.case_bed = in2

        # ref files
        fa_in = args.fa
        if fa_in and os.path.exists(fa_in):
            if not fa_in.endswith('.fa'):
                logger.error(f"Invalid fasta file provided, file extension should be .fa in plain text. input = {fa_in}")
                self.status = 1
        ref_fls = get_ref(self.organism, fa_in=fa_in, gtf=args.gtf)
        if ref_fls is None:
            logger.error("Error encountered while retrieving reference files")
            self.status = 1
            return
        self.fa = ref_fls['fa']
        self.gtf = ref_fls['gtf']
        if not os.path.exists(self.gtf):
            logger.error(f"GTF file not available: {self.gtf}")
            self.status = 1
        
        # verfify the normalization factors
        factors = [
            [args.f1, self.control_bed, 'condition1'],
            [args.f2, self.case_bed, 'condition2']
        ]
        
        if self.case_bed:
            if bool(args.f1) != bool(args.f2):
                logger.error("Please provide normalization factors for both or none of the conditions")
                self.status = 1
        for factor, bed, lb in factors:
            n_bed = len(bed) if bed else 0
            if factor:
                if len(factor) != n_bed:
                    logger.error(f"Number of normalization factors provided for {lb} does not match the number of input files")
                    self.status = 1
                for n, f in enumerate(factor):
                    if f == 0:
                        logger.error(f"{lb} - factore {n + 1} / {len(factor)}: normalize factor could not be zero!")
                        self.status = 1
        
        if self.status:
            return
        
        self.ref = ref_fls
        self.input_fls = self.control_bed + (self.case_bed or []) # fn_lb, fn_bed
        
        self.out_fls = {
            'bed_peaks': {},
            'fdr': {},
            'count_pp_gb': os.path.join(self.inter_dir, 'count_pp_gb.txt')
        }
        
        
        if self.longerna:
            self.n_gene_cols = 1
            self.gene_cols = ['long-eRNA']
            self.longerna_flag_str = '-long-eRNA'
        else:
            self.n_gene_cols = 2
            self.gene_cols = ['Transcript', 'Gene']
            self.longerna_flag_str = ''
        
        self.skip_get_mapped_reads = 0
        self.skip_count_pp_gb = 0

        fn_lb_uniq = set()
        fn_lb_dup = {}

        for fn_lb, fn_bed in self.input_fls:
            # tssCount, tssLength, tssRatio, tssSummit_pos, genebodyCount, genebodyLength, genebodyRatio, tss_vs_genebody_ratio
            if fn_lb not in fn_lb_uniq:
                fn_lb_uniq.add(fn_lb)
            else:
                fn_lb_dup.setdefault(fn_lb, []).append(fn_bed)
                continue
            header_bed_peaks = self.gene_cols + ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
            fn_bed_peaks = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_raw.txt')
            file_exist_flag = 0
            # if os.path.exists(fn_bed_peaks) and not self.overwrite_intermediate:
            #     # logger.warning(f"File already exists: {fn_bed_peaks}, wil reuse it")
            #     file_exist_flag = 1
            #     fh_bed_peaks = open(fn_bed_peaks, 'r')
            # else:
            
            fh_bed_peaks = open(fn_bed_peaks, 'w')
            fh_bed_peaks.write('\t'.join(header_bed_peaks) + '\n')
            self.out_fls['bed_peaks'][fn_lb] = {'fn': fn_bed_peaks, 'fh': fh_bed_peaks, 'header': header_bed_peaks, 'exist': file_exist_flag}
            
            header_fdr = header_bed_peaks + ['pvalue', 'FDR']
            fn_fdr = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_FDR.txt')
            self.out_fls['fdr'][fn_lb] = {'fn': fn_fdr, 'header': header_fdr}
        
        
        if fn_lb_dup:
            tmp = json.dumps(fn_lb_dup, indent=4)
            logger.error(f'Duplicate file labels found in input file list, please make the file name unique:\n{tmp}')
            sys.exit(1)
        
        # config
        self.config = {
            'mapped_sites_only': 0, # 1: only mapped sites are used to calculate the mappable sites, 0: all ATCG bases are used regardless of whether they are covered by reads. in the perl code, it is always set to 0
            'pro_up': args.up, # define the upstream of TSS as promoter (bp, default: 500)
            'pro_down': args.down, # define the downstream of TSS as promoter (bp, default: 500)
            'gb_start': args.gb, # define the start of gene body density calculation (bp, default: 1000)
            'min_gene_len': args.min_gene_len, # define the minimum length of a gene to perform the analysis (bp, default: 1000). Genes with length less than this will be ignored
            'window_size': args.window, # define the window size (bp, default: 50)
            'step': args.step, # define the step size (bp, default: 5)
            'normalize_factor1': args.f1, # normalization factors for samples of condition1
            'normalize_factor2': args.f2, # normalization factors for samples of condition2
        }
        

def bench(s, lb):
    time_cost.setdefault(lb, 0)
    now = time.time()
    time_cost[lb] += now - s
    return now


def test_get_peak(gene_info, fn_lb, fn_bed, pwout, bin_size):
    count_per_base, count_bin = pre_count_for_bed(fn_lb, fn_bed, pwout, bin_size)
    chr_, strand, gene_raw_s, pp_start, pp_end, gb_start, gb_end, strand_idx, gb_len_mappable, gene_seq = [gene_info[_] for _ in ['chr', 'strand', 'start', 'pp_start', 'pp_end', 'gb_start', 'gb_end', 'strand_idx', 'gb_len_mappable', 'gene_seq']]
    gbc, gbd, pp_res = get_peak(count_per_base, count_bin, chr_, strand, gene_raw_s, strand_idx, pp_start, pp_end, gb_start, gb_end, gb_len_mappable, gene_seq,  50, 5, bin_size, {})
    
    logger.debug(f'gbc = {gbc}, gbd={gbd}')
    tmp = json.dumps(pp_res, indent=3)
    logger.debug(tmp)


def process_bed_files(analysis, fls, gtf_info, fa_idx, fh_fa, reuse_pre_count=False):
    
    invalid_chr_transcript = 0
    bin_size = analysis.bin_size
    # ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
    fail_to_retrieve_seq = 0
    pwout = analysis.out_dir
    window_size, step_size = analysis.config['window_size'], analysis.config['step']
    pw_bed = args.pw_bed
    dummy_seq = 'ATCG'
    n = 0
    prev = [time.time(), 0] # time, get_mapped_reads_count count
    section_size = 1000 # print log every 1000 loops

    # add more values to gtf_info
    logger.debug(f'adding more values to gtf_info')
    s_before_loop = time.time()
    s = s_before_loop
    ts_list = list(gtf_info)

    pp_str = {ts: [] for ts in ts_list}
    gb_str = {ts: [] for ts in ts_list}
    
    # seq_pool = {} # key = transcript_id
    # count_pool = {}
    # prev_peak_pool = {}
    # genes_with_N = set()
    
    for fn_lb, fn_bed in fls:
        count_per_base, count_bin = pre_count_for_bed(fn_lb, fn_bed, pw_bed, bin_size, reuse=reuse_pre_count)
        prev_peak = {}
        
        fh_bed_peaks = analysis.out_fls['bed_peaks'][fn_lb]['fh']
        for transcript_id in ts_list:
            gene_info = gtf_info[transcript_id]
            chr_, strand, gene_raw_s, pp_start, pp_end, gb_start, gb_end, strand_idx, gb_len_mappable, gene_seq = [gene_info[_] for _ in ['chr', 'strand', 'start', 'pp_start', 'pp_end', 'gb_start', 'gb_end', 'strand_idx', 'gb_len_mappable', 'gene_seq']]

            n += 1
            if n % section_size == 0:
                now = time.time()
                prev_time, prev_count = prev
                time_gap = now - prev_time
                speed = time_gap * 1000  / (n - prev_count)
                logger.info(f'{n/1000}k - get_mapped_reads_count count, time_gap={time_gap:.2}s, speed={speed:.3f}ms')
                prev = [now, n]

            # pp_res: {'ppc': prev[0], 'ppd': ppd, 'mappable_sites': mappable_sites, 'summit_pos': summit_pos_str, 'summit_count': summit_count}
            s = time.time()
            gbc, gbd, pp_res = get_peak(count_per_base, count_bin, chr_, strand, gene_raw_s, strand_idx, pp_start, pp_end, gb_start, gb_end, gb_len_mappable, gene_seq,  window_size, step_size, bin_size, prev_peak)
            s = bench(s, 'get_peak')

            pp_str[transcript_id].append(str(pp_res['ppc']))
            gb_str[transcript_id] += [str(gbc), str(gbd)]
            
            s = bench(s, 'save pp_str')
            # write to file
            # row = ['Transcript', 'Gene', 'tssCount', 'tssLength', 'tssRatio', 'tssSummit_pos', 'genebodyCount', 'genebodyLength', 'genebodyRatio', 'tss_vs_genebody_ratio']
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gene_info['gene_name'])
            row += [pp_res['ppc'], pp_res['mappable_sites'], pp_res['ppd'], pp_res['summit_pos'], gbc, gb_len_mappable, gbd]
            pro_vs_pb = round(pp_res['ppd'] / gbd, 4) if gbd else 'NA'
            row.append(pro_vs_pb)
            s = bench(s, 'prepare row')
            
            print('\t'.join(map(str, row)), file=fh_bed_peaks)
            s = bench(s, 'write to file')

    if fail_to_retrieve_seq:
        logger.warning(f'failed to retrieve sequence for {fail_to_retrieve_seq} genes')

    if invalid_chr_transcript:
        logger.warning(f"Number of genes with invalid chromosome = {invalid_chr_transcript}")

    # tmp = json.dumps(time_cost, indent=3)
    # logger.info(tmp)

    return pp_str, gb_str

def collect_previous_count(analysis, fls):
    """
    get the previous done _raw, _FDR file
    """
    pp_str = {} # key = transcript_id, v = the peak count in promoter region for each file. 
    gb_str = {} # key = transcript_id, v = [the peak count in gene body region, ratio] for each file.

    for fn_lb, _ in fls:
        fh_bed_peaks = analysis.out_fls['bed_peaks'][fn_lb]['fh']
        # Transcript	Gene	ppc	ppm	ppd	pps	gbc	gbm	gbd	pauseIndex
        fh_bed_peaks.seek(0)
        header = fh_bed_peaks.readline()[:-1].split('\t') # header
        idx_ppc = header.index('ppc')
        for line in fh_bed_peaks:
            row = line.strip().split('\t')
            transcript_id = row[0]
            ppc, gbc, gbd = [row[_ + idx_ppc] for _ in [0, 4, 6]]
            pp_str.setdefault(transcript_id, [])
            pp_str[transcript_id].append(ppc)
            gb_str.setdefault(transcript_id, [])
            gb_str[transcript_id] += [gbc, gbd]

    return pp_str, gb_str


def main(args):
    """
    args: an object with attributes equiv to the argparse object, must specify 
    """

    # check if the attributes are present
    defined_attrs = vars(args)
    required_attrs = {'in1', 'pwout', 'organism'}
    
    optional_attrs = {
        'in2': None,
        'gtf': None,
        'fa': None,
        'f1': None,
        'f2': None,
        'up': 500,
        'down': 500,
        'gb': 1000,
        'min_gene_len': 1000,
        'window': 50,
        'step': 5,
        'bench': False,
        'ignore': False,
    }
    missing = required_attrs - set(defined_attrs)
    exist = required_attrs & set(defined_attrs)
    for attr in exist:
        if getattr(args, attr) is None:
            missing.add(attr)
    if missing:
        logger.error(f"Missing required attributes: {sorted(missing)}")
        return
    for attr, default in optional_attrs.items():
        if attr not in defined_attrs:
            setattr(args, attr, default)
        
    benchmode = args.bench
    # check dependencies
    dependency_status = check_dependency()
    if dependency_status:
        logger.error("Dependency check failed")
        sys.exit(1)

    analysis = Analysis(args)
    if analysis.status:
        logger.error("Exit now")
        sys.exit(1)
    
    logger.debug(vars(args))

    logger.debug(f'skip_count_pp_gb = {analysis.skip_count_pp_gb}, skip_get_mapped_reads = {analysis.skip_get_mapped_reads}')

    rep1 = len(analysis.control_bed)
    rep2 = len(analysis.case_bed) if analysis.case_bed else 0
    window_size = analysis.config['window_size']
    factor1 = analysis.config['normalize_factor1']
    factor2 = analysis.config['normalize_factor2']
    factor_flag = 0 if factor1 is None else 1

    pro_up, pro_down, gb_down_distance, min_gene_len = [analysis.config[_] for _ in ['pro_up', 'pro_down', 'gb_start', 'min_gene_len']]


    # process the GTF file
    logger.info(f"Processing GTF file: {analysis.ref['gtf']}")
    gtf_info, fn_tss, err = process_gtf(analysis.ref['gtf'], pwout=analysis.out_dir)

    # if args.testfunc:
    #     logger.warning('Debugging mode')
    #     transcript_id1 = 'NM_000051.3'
    #     transcript_id2 = 'NM_000051'
    #     gene_info = gtf_info.get(transcript_id1) or gtf_info[transcript_id2]
    #     gene_info = add_value_to_gtf(gene_info, pro_up, pro_down, gb_down_distance) 
    #     logger.debug(gene_info)
    #     fn_lb, fn_bed = analysis.input_fls[0]
        
    #     logger.debug(fn_lb)
    #     test_get_peak(gene_info, fn_lb, fn_bed, analysis.out_dir, bin_size=200)
    #     logger.debug('done')
    #     sys.exit(0)



    err_total = sum(err.values())
    if err_total:
        logger.warning(f'Error encountered while processing GTF file')
        print(json.dumps(err, indent=4))
    len1 = len(gtf_info)
    if len1 == 0:
        logger.error("No gene information found in the GTF file")
        sys.exit(1)


    # filter out the genes with length less than the minimum gene length
    gtf_info_new = {}
    if benchmode:
        ct = 0
        bench_max_gene = 3000
        logger.info(f"Bench mode is on, only process the first {bench_max_gene} genes")

    # prepare fasta file
    fn_fa = analysis.ref['fa']
    fa_idx = build_idx_for_fa(fn_fa)
    fh_fa = open(fn_fa, 'r')

    short_genes = 0

    tmp = sorted(gtf_info)
    # for k, v in gtf_info.items():
    for k in tmp:
        v = gtf_info[k]
        if v['chr'] not in fa_idx:
            continue
        if k[0] != 'N' and k[:3] != 'ENS':
            # unkown_transcript += 1
            # unkown_transcript_list.append(k)
            continue
        if v['end'] - v['start'] + 1 >= min_gene_len:
            gtf_info_new[k] = add_value_to_gtf(v, pro_up, pro_down, gb_down_distance) 
            if benchmode:
                ct += 1
                if ct == bench_max_gene:
                    break
        else:
            short_genes += 1
            
    if short_genes:
        logger.warning(f"Number of genes with length less than the minimum gene length ({min_gene_len}): n= {short_genes}, current total genes: {len(gtf_info_new)}")
    else:
        logger.info(f"Total number of genes after filtering: {len(gtf_info_new)}")
    gtf_info = gtf_info_new
    fls = analysis.input_fls # element is [fn_lb, fn_bed]
    sam_order = [_[0] for _ in fls]
    n_gene_cols = analysis.n_gene_cols
    
    reuse_pre_count = not vars(args).get('ignore', False) # modify here
    if not analysis.skip_get_mapped_reads:
        logger.info(f'Getting pp_gb count')
        pp_str, gb_str = process_bed_files(analysis, fls, gtf_info, fa_idx, fh_fa, reuse_pre_count=reuse_pre_count)

        # close file handle
        for fn_lb, fn_bed in fls:
            analysis.out_fls['bed_peaks'][fn_lb]['fh'].close()
        fh_fa.close()
        # calculate FDR
        logger.info('Calculating FDR')
        pause_index_str = {}
        for fn_lb, fn_bed in fls:
            fn_peak = analysis.out_fls['bed_peaks'][fn_lb]['fn']
            fn_fdr = analysis.out_fls['fdr'][fn_lb]['fn']
            pause_index_str = get_FDR_per_sample(fn_lb, fn_peak, fn_fdr,  pause_index_str)

    else:
        logger.debug('collecting pause index from FDR files')
        pause_index_str = {}  # key = transcript_id
        for fn_lb, _ in fls:
            fn_peak = analysis.out_fls['bed_peaks'][fn_lb]['fn']
            fn_fdr = analysis.out_fls['fdr'][fn_lb]['fn']
            if not os.path.exists(fn_fdr):
                pause_index_str = get_FDR_per_sample(fn_lb, fn_peak, fn_fdr,  pause_index_str)
            with open(fn_fdr) as f:
                header = f.readline().strip().split('\t')
                idx_pindex, idx_pval, idx_fdr = header.index('pauseIndex'), header.index('pvalue'), header.index('FDR')
                for i in f:
                    line = i[:-1].split('\t')
                    transcript_id, pindex, pval, fdr = [line[_] for _ in [0, idx_pindex, idx_pval, idx_fdr]]
                    pause_index_str.setdefault(transcript_id, {})[fn_lb] = [pindex, pval, fdr]
    
    pause_index_str_new = {}
    for transcript_id, v1 in pause_index_str.items():
        v2 = []
        for isam in sam_order:
            v3 = v1.get(isam, ['NA', 'NA', 'NA'])
            v2 += v3
        pause_index_str_new[transcript_id] = v2
    pause_index_str = pause_index_str_new
    
 
    # count_pp_gb
    fn_count_pp_gb = analysis.out_fls['count_pp_gb']
    # analysis.skip_count_pp_gb = 0
    if not analysis.skip_count_pp_gb:
        if analysis.skip_get_mapped_reads:
            logger.debug('Re-use previous mapped reads count')
            pp_str, gb_str = collect_previous_count(analysis, fls)
        
        logger.info('Building count_pp_gb')
        prefix = 'longeRNA-' if analysis.longerna else ''
        header = analysis.gene_cols + ['chr', 'start', 'end', 'strand']
        header_pp = []
        header_gb = []
        for fn_lb, fn_bed in fls:
            header_gb += [f'gbc_{fn_lb}', f'gbd_{fn_lb}']
            header_pp.append(f'ppc_{fn_lb}')
        header += header_pp + header_gb
        
        with open(fn_count_pp_gb, 'w') as o:
            print('\t'.join(header), file=o)
            for transcript_id in pp_str:
                gene_info = gtf_info[transcript_id]
                row = [transcript_id]
                if not analysis.longerna:
                    row.append(gene_info['gene_name'])
                row +=  [gene_info['chr'], str(gene_info['start']), str(gene_info['end']), gene_info['strand']]
                row += pp_str[transcript_id] + gb_str[transcript_id]
                print('\t'.join(row), file=o)

    # if time_cost:
    #     tmp = json.dumps(time_cost, indent=3)
    #     print(tmp)
    
    fno_prefix = 'longeRNA-' if analysis.longerna else ''

    # modify here
    logger.info('Change_pp_gb')
    # logger.warning('modify here')
    
    change_pp_gb(n_gene_cols, fn_count_pp_gb, analysis.out_dir, rep1, rep2, window_size, factor1=factor1, factor2=factor2, factor_flag=factor_flag)
    
    # dump pindex.txt
    # pause_index[fn_lb][transcript_id] = [pro_vs_pb, pvalue, fdr]
    
    logger.debug('Dump pindex.txt')
    header_extra = []
    for fn_lb, fn_bed in fls:
        header_extra += [f'{fn_lb}-pindex', f'{fn_lb}-pvalue', f'{fn_lb}-FDR']
    header = analysis.gene_cols + header_extra
    fno =  os.path.join(analysis.known_gene_dir, fno_prefix + 'pindex.txt') 
    with open(fno, 'w') as o:
        print('\t'.join(header), file=o)
        for transcript_id, data in pause_index_str.items():
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gtf_info[transcript_id]['gene_name'])
            row += data
            print('\t'.join(map(str, row)), file=o)

    # boxplot
    logger.info(f'plotting boxplot')
    draw_box_plot(n_gene_cols, analysis.out_dir, 'boxplot', rep1, rep2)

    if rep2 > 0:
        # change_pindex
        logger.debug('Running change_pindex')
        change_pindex(fno_prefix, n_gene_cols, fn_count_pp_gb, analysis.out_dir, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0)
        
        
        # simple_heatmap
        logger.info(f'plotting pindex heatmap')
        draw_heatmap_pindex(n_gene_cols, analysis.out_dir)
        

        # heatmap
        logger.info(f'plotting heatmap for pp_change')
        # "perl heatmap.pl -w $out_dir -i $list -in1 $cond1_str -in2 $cond2_str -m $genome -n $tname";

        draw_heatmap_pp_change(n_gene_cols, analysis.out_dir,  fls_ctrl=analysis.control_bed, fls_case=analysis.case_bed, fn_tss=fn_tss, region_size=5000, bin_size=200, outname='heatmap')

        # modify here
        # logger.info('debug exit')
        # sys.exit(0)

        # # get alternative isoforms
        # logger.info('Getting alternative TSS isoforms')
        # # get_alternative_isoform_across_conditions(fn_norm_count, out_dir, rep1, rep2)
        # fn_pp_gb_count_norm = os.path.join(analysis.known_gene_dir, prefix + 'normalized_pp_gb.txt')
        # if not os.path.exists(fn_pp_gb_count_norm):
        #     logger.error(f"Normalized pp_gb file not found: {fn_pp_gb_count_norm}")
        #     sys.exit(1)
        # get_alternative_isoform_across_conditions(fn_pp_gb_count_norm, analysis.out_dir, rep1, rep2)


def get_new_args(args, update_dict):
    args_dict = vars(args)
    args_dict.update(update_dict)
    return SimpleNamespace(**args_dict)

    
if __name__ == "__main__":
    verbose = args.verbose
    if verbose:
        # set the "console" file handler to debug level
        # get the file handles
        handlers = logger.handlers
        for h in handlers:
            if h.name == 'console':
                h.setLevel('DEBUG')
                break
        logger.debug(f'Verbose mode is on')
    pwout_raw = os.path.realpath(args.pwout)
    args.pwout = pwout_raw
    args.pw_bed = f'{pwout_raw}/bed'
    if not os.path.exists(args.pw_bed):
        os.makedirs(args.pw_bed, exist_ok=True)
    args.organism = args.organism

    if args.in1 is not None:
        # logger.info('Processing input files')
        main(args)
    elif args.design_table is not None:
        group_info = {}  # key = group name, v = file list
        comparison = []
        groups_in_comparison = set()
        group_info['null'] = []
        with open(args.design_table) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if line.startswith('@@'):
                    tmp = line[2:].strip().split('\t')
                    if len(tmp) != 2:
                        logger.error(f"Invalid comparison definition: {line}, should have 2 columns, col1 = case group, col2 = control group")
                        sys.exit(1)
                    comparison.append(tmp)
                    groups_in_comparison |= set(tmp)
                else:
                    tmp = line.strip().split('\t')
                    if len(tmp) != 2:
                        logger.error(f"Invalid line: {line}, should have 2 columns, col1 = file path, col2 = group name")
                        sys.exit(1)
                    group_info.setdefault(tmp[1], []).append(tmp[0])
            
            group_not_defined = groups_in_comparison - set(group_info)
            if group_not_defined:
                logger.error(f"Groups in comparison not defined in the sample section: {sorted(group_not_defined)}")
                sys.exit(1)
            
            # if no comparison defined, process each group separately
            if len(comparison) == 0:
                logger.error(f'No comparison defined in {args.design_table}, please use @@ to define the comparison .If need to process a group separately, use @@null\\tgroup_name')
                sys.exit(1)
            else:
                for case, control in comparison:
                    update_dict = {}
                    update_dict['in1'] = group_info[control]
                    update_dict['in2'] = group_info[case] or None
                    comp_str = control if case == 'null' else f'{case}_vs_{control}'
                    out_dir = os.path.join(pwout_raw, comp_str)
                    update_dict['pwout'] = out_dir
                    
                    args_new = get_new_args(args, update_dict)
                    # logger.info(vars(args_new))
                    logger.info(f'Processing {comp_str}')
                    main(args_new)
    else:
        logger.error("No input files provided, either use -in1 / in2 or provide a design table")
        sys.exit(1)

