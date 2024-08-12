#! /usr/bin/env python3
"""
process PROseq data
"""
import time
s = time.time()
import sys, os, re
import pickle
import json
from types import SimpleNamespace
from utils import check_dependency, run_shell, process_input, get_seqence_from_fa, build_idx_for_fa, get_ref, process_gtf, get_peak, get_summit, groseq_fisherexact_pausing_for_gene, change_pp_gb, change_pindex, draw_box_plot, draw_heatmap_pindex, draw_heatmap_pp_change, calc_FDR, get_alternative_isoform_across_conditions


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

def getargs():
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('design_table',  help="""Optional, desgin table in tsv format. 2 sections, first section is the sample information, 2 columns. col1=bed/bam full file path, col2 = group name. Second section is the comparison information, 2 columns, col1 should start with @@ used to identify this line as comparison definition. col1=group name used as case, col2 = group name used as control. e.g. @@case\\tcontrol. If only need to process a single group, use @@null\\tgroup_name""", nargs='?')
    ps.add_argument('-in1', help="""required, read alignment files in bed (6 columns) or bam format for condition1, separated by space""", nargs='+')
    ps.add_argument('-in2', help="""read alignment files in bed (6 columns) or bam format for condition2, separated by space""", nargs='*')
    ps.add_argument('-o', help="""required, output/work directory""", required=True)
    ps.add_argument('-organism', '-m', '-org',  help="""required, define the genome. Valid is hg19, hg38, mm10, dm3, dm6, ce10, or danRer10. default: hg19""", required=True)

    ps.add_argument('-gtf', help="""user specified GTF file, if not specified, will use the default GTF file for the organism""")
    ps.add_argument('-fa', help="""Full path for the fasta file for the genome. If not specified, will search under the fa folder under the package directory. e.g. organism is hg19, then the default fasta file should be fa/hg19.fa""")
    ps.add_argument('-f1', help="""normalization factors for samples of condition1, separated by space. If no normalization factor is specified, we will use DEseq2 default method to normalize""", nargs='*', type=float)
    ps.add_argument('-f2', help="""normalization factors for samples of condition2, same as -f1""", nargs='*', type=float)
    ps.add_argument('-up', '-u', help="""define the upstream of TSS as promoter (bp, default: 500)""", type=int, default=500)
    ps.add_argument('-down', '-d', help="""define the downstream of TSS as promoter (bp, default: 500)""", type=int, default=500)
    ps.add_argument('-gb', '-b', help="""define the start of gene body density calculation (bp, default: 1000)""", type=int, default=1000)
    ps.add_argument('-min_gene_len', '-l', help="""define the minimum length of a gene to perform the analysis (bp, default: 1000). Genes with length less than this will be ignored""", type=int, default=1000)
    ps.add_argument('-window', '-w', help="""define the window size (bp, default: 50)""", type=int, default=50)
    ps.add_argument('-step', '-s', help="""define the step size (bp, default: 5)""", type=int, default=5)
    ps.add_argument('-bench', '-test', help="""bench mode, only process the first 10 genes for testing purpose""", action='store_true')
    ps.add_argument('-overwrite', help="""if the mapped reads count file already exists, overwrite it, default is reuse the old ones""", action='store_true')
    ps.add_argument('-verbose', '-v', help="""verbose mode, print more information""", action='store_true')
    args = ps.parse_args()
    
    return args


def seek(pos, padding=1000):
    
    query_region_size = 2000
    print('regions before:')
    f.seek(pos - padding)
    print(f.read(padding))
    
    print('\n\nqueried region')
    f.seek(pos)
    print(f.read(query_region_size))
    
    print('\n\nregions after:')
    f.seek(pos + query_region_size)
    print(f.read(padding))


class Analysis:
    def __init__(self, args, is_long_eRNA=False):
        
        # folders
        self.bin_dir = os.path.dirname(os.path.realpath(__file__))
        self.out_dir = args.o or os.getcwd()
        self.inter_dir = os.path.join(self.out_dir, 'intermediate')
        self.known_gene_dir = os.path.join(self.out_dir, 'known_gene')
        self.longerna = is_long_eRNA # toggle if this is long eRNA
        self.overwrite_intermediate = args.overwrite
        
        for lb, d in zip(['Output', 'Intermediate', 'Known Genes', 'bed file folder'], [self.out_dir, self.inter_dir, self.known_gene_dir, f'{self.inter_dir}/bed']):
            if not os.path.exists(d):
                logger.info(f'Making {lb} directory: {d}')
                os.makedirs(d)

        self.status = 0

        self.organism = args.organism
        if self.organism not in {'hg19', 'hg38', 'mm10', 'dm3', 'dm6', 'ce10', 'danRer10'}:
            logger.error(f"Invalid organism provided: {self.organism}")
            self.status = 1

        bed_idx_step = 50
        self.bed_idx_step = bed_idx_step
        
        # input files
        in1 = process_input(self.out_dir, args.in1, bed_idx_step=bed_idx_step) # return = fn_lb, fn_bed_count, idx_d, res_chunk
        if in1 is None:
            logger.error("Invalid input files provided for condition1")
            self.status = 1
        in2 = process_input(self.out_dir, args.in2, bed_idx_step=bed_idx_step)
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
        self.input_fls = self.control_bed + (self.case_bed or [])
        
        self.out_fls = {
            'bed_peaks': {},
            'fdr': {},
            'count_file': {},
        }
        
        
        if self.longerna:
            self.n_gene_cols = 1
            self.gene_cols = ['long-eRNA']
            self.longerna_flag_str = '-long-eRNA'
        else:
            self.n_gene_cols = 2
            self.gene_cols = ['Transcript', 'Gene']
            self.longerna_flag_str = ''
        
        self.skip_get_mapped_reads = 1
        # fn_lb, fn_bed_count, idx_d, res_chunk
        for fn_lb, fn_bed_count, idx_d, res_chunk in self.input_fls:
            # tssCount, tssLength, tssRatio, tssSummit_pos, genebodyCount, genebodyLength, genebodyRatio, tss_vs_genebody_ratio
            self.out_fls['count_file'][fn_lb] = {'fn': fn_bed_count, 'idx_d': idx_d, 'res_chunk': res_chunk}
            fh_count_file = open(fn_bed_count, 'r')
            self.out_dir['count_file'][fn_lb]['fh'] = fh_count_file
            
            header_bed_peaks = self.gene_cols + ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
            fn_bed_peaks = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_raw.txt')
            file_exist_flag = 0
            if os.path.exists(fn_bed_peaks) and not self.overwrite_intermediate:
                # logger.warning(f"File already exists: {fn_bed_peaks}, wil reuse it")
                file_exist_flag = 1
                fh_bed_peaks = open(fn_bed_peaks, 'r')
            else:
                fh_bed_peaks = open(fn_bed_peaks, 'w')
                fh_bed_peaks.write('\t'.join(header_bed_peaks) + '\n')
                self.skip_get_mapped_reads = 0
            self.out_fls['bed_peaks'][fn_lb] = {'fn': fn_bed_peaks, 'fh': fh_bed_peaks, 'header': header_bed_peaks, 'exist': file_exist_flag}
            
            header_fdr = header_bed_peaks + ['pvalue', 'FDR']
            fn_fdr = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_FDR.txt')
            fh_fdr = open(fn_fdr, 'w')
            fh_fdr.write('\t'.join(header_fdr) + '\n')
            self.out_fls['fdr'][fn_lb] = {'fn': fn_fdr, 'fh': fh_fdr, 'header': header_fdr}
        
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


def process_bed_files(analysis, fls, gtf_info, fa_idx, fh_fa, ):

    already_parsed_windows = {fn_lb: {} for fn_lb, fn, _2, _3 in fls} # fn_lb, fn_bed_count, idx_d, res_chunk
    n_transcript = 0
    invalid_chr_transcript = 0
    pp_str = {} # key = gene_id, v = the peak count in promoter region for each file. 
    gb_str = {} # key = gene_id, v = [the peak count in gene body region, ratio] for each file.
    
    # ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
    fn_peak_float_cols = {n_gene_cols + 2, n_gene_cols + 6} # used for output, to round the float values to 5 decimal places
    pause_pvalue = {}  # key = fn_bed, used for retrieving the pause raw count file, used as base to add 2 columns for p_value and FDR. v = list of pvalue in each row, used to calculate FDR
    gbc_sum = {} # key = fn_bed, v = sum of gene body count for each gene
    fn_peak_rows = {} # k1 = fn_bed, k2 = row_idx, v = row_str
    fail_to_retrieve_seq = 0


    for transcript_id, gene_info in gtf_info.items():
        # gene_info: {'chr': '', 'strand': '', 'gene_id': '', 'gene_name': '', 'start': 0, 'end': 0}
        chr_ = gene_info['chr']
        strand = gene_info['strand']
        n_transcript += 1
        if n_transcript % 1000 == 0:
            # modify here            
            logger.info(f'Processing gene {n_transcript} / {len(gtf_info)}')

        # skip the gene if the chromosome is not in the fasta file
        if chr_ not in fa_idx:
            invalid_chr_transcript += 1
            continue
        gene_raw_s, gene_raw_e = gene_info['start'], gene_info['end']
        up_sign = -1 if strand == '+' else 1
        down_sign = 1 if strand == '+' else -1
        gene_strand_s = gene_raw_s if strand == '+' else gene_raw_e
        gene_strand_e = gene_raw_e if strand == '+' else gene_raw_s
        s = bench(s, 'general_gene_info')
        
        gene_seq = get_seqence_from_fa(fa_idx, fh_fa, chr_, gene_raw_s, gene_raw_e)
        s = bench(s, 'get_gene_seq')
        if gene_seq == 1:
            fail_to_retrieve_seq += 1
            continue
        pp_str[transcript_id] = []
        gb_str[transcript_id] = []
        s = bench(s, 'before bed loop')
        for fn_lb, fn_bed_count, idx_d, res_chunk in fls:
            s = time.time()
            
            count_file_info = analysis.out_fls['count_file'][fn_lb] # fn, fh, idx_d, res_chunk
            # promoter peak
            promoter_s = gene_strand_s + up_sign * analysis.config['pro_up']
            promoter_e = gene_strand_s + down_sign * (analysis.config['pro_down'] - 1)
            regioin_s, region_e = sorted([promoter_s, promoter_e])

            pro_peak = get_peak(count_file_info, chr_, strand, regioin_s, regioin_e, gene_raw_s, gene_seq, analysis.config['window_size'], analysis.config['step'], analysis.config['mapped_sites_only'], time_cost, already_parsed_windows=already_parsed_windows[fn_lb])

            s = bench(s, 'get_pro_peak')
            # pro_peak = [max_read_count, mappable_sites, ratio, window_start]
            # but we need the summit position in this window
            if pro_peak[0] > 0:
                summit_pos = get_summit(count_file_info, strand, chr_, pro_peak[3], window_size=analysis.config['window_size'])
                pro_peak[3] = summit_pos
            s = bench(s, 'get_pro_summit')

            # gene body peak
            gene_body_s = gene_strand_s + down_sign * analysis.config['gb_start']
            gene_body_e = gene_strand_e
            region_s, region_e = sorted([gene_body_s, gene_body_e])

            gene_body_len = gene_raw_e - gene_raw_s + 1 - analysis.config['gb_start']
            gene_body_peak = get_peak(count_file_info, chr_, strand, region_s, region_e, gene_raw_s, gene_seq, gene_body_len, analysis.config['step'], analysis.config['mapped_sites_only'], time_cost, already_parsed_windows=already_parsed_windows[fn_lb])
            # gene_body_peak = [100, 1000, 0.1, 11111] # modify here, remove
            # s = bench(s, 'get_gb_peak')
            
            pp_str[transcript_id].append(str(pro_peak[0]))
            gb_str[transcript_id] += [str(gene_body_peak[0]), str(gene_body_peak[2])]
            
            # write to file
            fh_bed_peaks = analysis.out_fls['bed_peaks'][fn_lb]['fh']
            
            # row = ['Transcript', 'Gene', 'tssCount', 'tssLength', 'tssRatio', 'tssSummit_pos', 'genebodyCount', 'genebodyLength', 'genebodyRatio', 'tss_vs_genebody_ratio']
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gene_info['gene_name'])
            row += pro_peak + gene_body_peak[:3]
            pro_vs_pb = round(pro_peak[2] / gene_body_peak[2], 5) if gene_body_peak[2] else 'NA'
            row.append(pro_vs_pb)
            try:
                row_str = '\t'.join([str(round(_, 5)) if idx in fn_peak_float_cols and _ != 'NA' else str(_) for idx, _ in enumerate(row)])
            except:
                logger.error(f'invalid row, row = {list(enumerate(row))}, {[row[_] for _ in sorted(fn_peak_float_cols)]}')
                logger.error(f"transcript_id = {transcript_id}, gene_info = {gene_info}, pro_peak = {pro_peak}, gene_body_peak = {gene_body_peak}")
                
                sys.exit(1)
            print(row_str, file=fh_bed_peaks)
            
            # calculate p-value and FDR for certain rows
            ppc = pro_peak[0]
            gbd = gene_body_peak[2]
            gbc = gene_body_peak[0]
            gbc_sum.setdefault(fn_lb, 0)
            gbc_sum[fn_lb] += gbc
            
            # ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
            if ppc > 0:
                ppc, pp_len, gbc, gb_len = [row[n_gene_cols + _] for _ in [0, 1, 4, 5]] # 2,3,6,7 for normal genes
                pvalue = groseq_fisherexact_pausing_for_gene(ppc, pp_len, gbc, gb_len)
            else:
                pvalue = 'NA'
            # s = bench(s, 'fisher_exact')
            
            pause_pvalue.setdefault(fn_lb, []).append([pvalue, ppc, gbd])
            fn_peak_rows.setdefault(fn_lb, {'row_idx': 0})
            row_idx = fn_peak_rows[fn_lb]['row_idx']
            fn_peak_rows[fn_lb][row_idx] = [row_str, transcript_id, pro_vs_pb]
            fn_peak_rows[fn_lb]['row_idx'] += 1

    
    if fail_to_retrieve_seq:
        logger.warning(f"Failed to retrieve sequence for {fail_to_retrieve_seq} genes")
    
    if invalid_chr_transcript:
        logger.warning(f"Number of genes with invalid chromosome = {invalid_chr_transcript}")

    return pp_str, gb_str, pause_pvalue, gbc_sum, fn_peak_rows

def collect_previous_count(analysis, fls):
    """
    get the previous done _raw, _FDR file
    """
    pp_str = {} # key = transcript_id, v = the peak count in promoter region for each file. 
    gb_str = {} # key = transcript_id, v = [the peak count in gene body region, ratio] for each file.

    pause_pvalue = {}  # key = fn_bed, used for retrieving the pause raw count file, used as base to add 2 columns for p_value and FDR. v = list of pvalue in each row, used to calculate FDR

    gbc_sum = {} # key = fn_bed, v = sum of gene body count for each gene

    fn_peak_rows = {} # k1 = fn_bed, k2 = row_idx, v = [row_str, transcript_id, pro_vs_pb]

    fisher_cache = {} # key = fisher parameters concatenated, v = pvalue
    n_loop = 0
    fisher_run_count = 0
    prev = [time.time(), 0]  # the timestamp and the fisher count of prev logging
    for fn_lb, fn_bed_count, idx_d, res_chunk in fls:
        fh_bed_peaks = analysis.out_fls['bed_peaks'][fn_lb]['fh']
        # Transcript	Gene	ppc	ppm	ppd	pps	gbc	gbm	gbd	pauseIndex
        fh_bed_peaks.seek(0)
        header = fh_bed_peaks.readline().strip().split('\t') # header
        idx_ppc = header.index('ppc')
        gbc_sum.setdefault(fn_lb, 0)
        pause_pvalue.setdefault(fn_lb, [])
        fn_peak_rows.setdefault(fn_lb, {})
        for row_idx, line in enumerate(fh_bed_peaks, 1):
            n_loop += 1
            row_str = line.strip()
            row = line.strip().split('\t')
            transcript_id = row[0]
            ppc, pp_len, ppd, pps, gbc, gb_len, gbd, pro_vs_pb = row[idx_ppc:idx_ppc + 8]
            pp_str.setdefault(transcript_id, [])
            pp_str[transcript_id].append(ppc)
            gb_str.setdefault(transcript_id, [])
            gb_str[transcript_id] += [gbc, gbd]

            fisher_key = (ppc, pp_len, gbc, gb_len)
            ppc = int(ppc)
            gbc = int(gbc)
            pp_len = int(pp_len)
            gb_len = int(gb_len)
            gbd = float(gbd)
            gbc_sum[fn_lb] += gbc
            
            if all([ppc, pp_len, gbc, gb_len]):
                
                if fisher_key in fisher_cache:
                    pvalue = fisher_cache[fisher_key]
                else:
                    pvalue = groseq_fisherexact_pausing_for_gene(ppc, pp_len, gbc, gb_len)
                    # pvalue = 'NA'
                    fisher_cache[fisher_key] = pvalue
                    fisher_run_count += 1
            else:
                pvalue = 'NA'
            pause_pvalue[fn_lb].append([pvalue, ppc, gbd])
            fn_peak_rows[fn_lb][row_idx] = [row_str, transcript_id, pro_vs_pb]
            
            if n_loop % 10000 == 0:
                now = time.time()
                prev_time, prev_fisher_count = prev
                time_gap = now - prev_time
                speed_fisher = time_gap * 1000 / (fisher_run_count - prev_fisher_count)
                logger.info(f'{n_loop/10000}w - Fisher count: {fisher_run_count/10000}w, time_gap={time_gap:.2}s, speed_fisher={speed_fisher:.3f}ms')
                prev = [now, fisher_run_count]

    return pp_str, gb_str, pause_pvalue, gbc_sum, fn_peak_rows


def main(args=None):
    """
    args: an object with attributes equiv to the argparse object, must specify 
    """
    if args is None:
        args = getargs()
    else:
        # check if the attributes are present
        defined_attrs = vars(args)
        required_attrs = {'in1', 'o', 'organism'}
        
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
    # logger.info('building obj done')
    if analysis.status:
        logger.error("Exit now")
        sys.exit(1)

    # process the GTF file
    logger.info(f"Processing GTF file: {analysis.ref['gtf']}")
    gtf_info, err = process_gtf(analysis.ref['gtf'])
    # gtf_info: key = transcript_id, v = {'chr': '', 'strand': '', 'gene_id': '', 'gene_name': '', 'start': 0, 'end': 0}
    err_total = sum(err.values())
    if err_total:
        logger.warning(f'Error encountered while processing GTF file')
        print(json.dumps(err, indent=4))
    len1 = len(gtf_info)
    if len1 == 0:
        logger.error("No gene information found in the GTF file")
        sys.exit(1)
    logger.info(f"Total number of genes in the GTF file: {len1}")

    # filter out the genes with length less than the minimum gene length
    gtf_info_new = {}
    if benchmode:
        ct = 0
        bench_max_gene = 1000
        logger.info(f"Bench mode is on, only process the first {bench_max_gene} genes")

    short_genes = 0
    unkown_transcript = 0
    unkown_transcript_list = []
    # modify here, bench mode
    tmp = sorted(gtf_info)
    # for k, v in gtf_info.items():
    for k in tmp:
        v = gtf_info[k]
        if k[0] != 'N' and k[:3] != 'ENS':
            # unkown_transcript += 1
            # unkown_transcript_list.append(k)
            continue
        
        if v['end'] - v['start'] + 1 >= analysis.config['min_gene_len']:
            gtf_info_new[k] = v
            if benchmode:
                ct += 1
                if ct == bench_max_gene:
                    break
        else:
            short_genes += 1
            
    if short_genes:
        logger.warning(f"Number of genes with length less than the minimum gene length ({analysis.config['min_gene_len']}): {short_genes}, current total genes: {len(gtf_info_new)}")
    else:
        logger.info(f"Total number of genes after filtering: {len(gtf_info_new)}")
    gtf_info = gtf_info_new
    # logger.info(sorted(gtf_info)[:10])
    
    # with open('unkown_transcript.txt', 'w') as o:
    #     print('\n'.join(unkown_transcript_list), file=o)
    
    # sys.exit(1)
    
    # prepare fasta file
    fn_fa = analysis.ref['fa']
    n_gene_cols = analysis.n_gene_cols
    fh_fa = open(fn_fa, 'r')
    lb = fn_fa.rsplit('.', 1)[0]
    fn_idx = f'{lb}.idx.pkl'
    if os.path.exists(fn_idx) and os.path.getsize(fn_idx) > 100:
        with open(fn_idx, 'rb') as o:
            fa_idx = pickle.load(o)
    else:
        logger.info(f'building index for fasta file: {fn_fa}')
        fa_idx = build_idx_for_fa(fn_fa)
 
    fls = analysis.input_fls # element is [fn_bed, idx_bed, fn_lb]

    if not analysis.skip_get_mapped_reads:
        logger.error(f'fail to load prev mapped reads count')
        sys.exit(1)
        pp_str, gb_str, pause_pvalue, gbc_sum, fn_peak_rows = process_bed_files(gtf_info, fa_idx)
    else:
        logger.debug('Re-use previous mapped reads count')
        pp_str, gb_str, pause_pvalue, gbc_sum, fn_peak_rows = collect_previous_count(analysis, fls)


    logger.info('debug')
    sys.exit(1)
    logger.info(f'processing bed files done')
    # close file handle
    for fn_lb, fn_bed_count, idx_d, res_chunk in fls:
        analysis.out_fls['bed_peaks'][fn_lb]['fh'].close()
        analysis.out_fls['count_file'][fn_lb]['fh'].close()
    fh_fa.close()
    
    # s = bench(s, 'close file handle')
    
    # calculate FDR
    logger.info('Calculating FDR')
    pause_index = {} # key = transcript_id, v = [pause_index, p_value, fdr], pause_index = ratio_pp / ratio_gb

    for fn_lb, pvalues in pause_pvalue.items():
        pause_index[fn_lb] = {}
        file_gbc_sum = gbc_sum[fn_lb]
        valid_rows = [[], []] # 0 = pvalue, 1 = idx
        invalid_rows = set()
        
        # only filter for regular genes, long eRNA will keep all rows
        if analysis.longerna:
            valid_rows[0] = [_[0] for _ in pvalues]
            valid_rows[1] = list(range(len(pvalues)))
        else:
            for idx, (pvalue, ppc, gbd) in enumerate(pvalues):
                if ppc > 0 and gbd * 10_000_000/file_gbc_sum > 0.004 and pvalue is not None:
                    valid_rows[0].append(pvalue)
                    valid_rows[1].append(idx)
                else:
                    invalid_rows.add(idx)
        if len(valid_rows[0]) < 10:
            logger.warning(f'Number of rows with pvalue  less than 10, fn_lb = {fn_lb}')
        fdr_all = calc_FDR(valid_rows[0])
        pvalue_fdr_str = sorted(zip(valid_rows[1], valid_rows[0], fdr_all), key=lambda x: x[2]) # sort by FDR

        fh_fdr = analysis.out_fls['fdr'][fn_lb]['fh']
        header = analysis.out_fls['fdr'][fn_lb]['header']
        print('\t'.join(header), file=fh_fdr)
        for idx, pvalue, fdr in pvalue_fdr_str:
            row_str, transcript_id, pro_vs_pb = fn_peak_rows[fn_lb][idx]
            print(row_str + f'\t{pvalue}\t{fdr}', file=fh_fdr)
            pause_index[fn_lb][transcript_id] = [pro_vs_pb, pvalue, fdr]
        for idx in invalid_rows:
            row_str, transcript_id, pro_vs_pb = fn_peak_rows[fn_lb][idx]
            print(row_str + '\tNA\tNA', file=fh_fdr)
            pause_index[fn_lb][transcript_id] = [pro_vs_pb, 'NA', 'NA']
        fh_fdr.close()
    # s = bench(s, 'calculate_FDR')
    
    logger.error('exit')
    sys.exit(1)
    
    
    # count_pp_gb
    logger.info('Building count_pp_gb')
    prefix = 'longeRNA-' if analysis.longerna else ''
    fn_count_pp_gb = os.path.join(analysis.inter_dir, prefix + 'count_pp_gb.txt')
    header = analysis.gene_cols + ['chr', 'start', 'end', 'strand']
    header_pp = []
    header_gb = []
    for fn_lb, fn_bed_count, idx_d, res_chunk in fls:
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
            row +=  [gene_info['chr'], gene_info['start'], gene_info['end'], gene_info['strand']]
            row += pp_str[transcript_id] + gb_str[transcript_id]
            print('\t'.join(row), file=o)
    # s = bench(s, 'count_pp_gb')
    
    # run change_pp_gb
    rep1 = len(analysis.control_bed)
    rep2 = len(analysis.case_bed) if analysis.case_bed else 0
    window_size = analysis.config['window_size']
    factor1 = analysis.config['normalize_factor1']
    factor2 = analysis.config['normalize_factor2']
    factor_flag = 0 if factor1 is None else 1
    change_pp_gb(n_gene_cols, fn_count_pp_gb, analysis.out_dir, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0)
    
    
    # get alternative isoforms
    if rep2 > 0:
        logger.info('Getting alternative isoforms')
        # get_alternative_isoform_across_conditions(fn_norm_count, out_dir, rep1, rep2)
        fn_pp_gb_count_norm = os.path.join(analysis.known_gene_dir, prefix + 'normalized_pp_gb.txt')
        if not os.path.exists(fn_pp_gb_count_norm):
            logger.error(f"Normalized pp_gb file not found: {fn_pp_gb_count_norm}")
            sys.exit(1)
        get_alternative_isoform_across_conditions(fn_pp_gb_count_norm, analysis.out_dir, rep1, rep2)
    
    # change_pindex
    fno_prefix = 'longeRNA-' if analysis.longerna else ''
    change_pindex(fno_prefix, n_gene_cols, fn_count_pp_gb, analysis.out_dir, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0)
    
    # dump pindex.txt
    # pause_index[fn_lb][transcript_id] = [pro_vs_pb, pvalue, fdr]
    header_extra = []
    for fn_lb, fn_bed_count, idx_d, res_chunk in fls:
        header_extra += [f'{fn_lb}-pindex', f'{fn_lb}-pvalue', f'{fn_lb}-FDR']
    header = analysis.gene_cols + header_extra
    fno =  os.path.join(analysis.known_gene_dir, prefix + 'pindex.txt') 
    with open(fno, 'w') as o:
        print('\t'.join(header), file=o)
        for transcript_id, gene_info in gtf_info.items():
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gene_info['gene_name'])
            for fn_lb, fn_bed_count, idx_d, res_chunk in fls:
                for k in pause_index[fn_lb][transcript_id]:
                    row += k
            print('\t'.join(map(str, row)), file=o)

    # boxplot
    logger.info(f'plotting boxplot')
    draw_box_plot(n_gene_cols, analysis.out_dir, 'boxplot', rep1, rep2)
    
    # simple_heatmap
    logger.info(f'plotting pindex heatmap')
    draw_heatmap_pindex(n_gene_cols, analysis.out_dir)

    # heatmap
    logger.info(f'plotting heatmap for pp_change')
    fn_glist = os.path.join(analysis.inter_dir, prefix + "genes_for_heatmap.txt")
    fn_pp_change = os.path.join(analysis.known_gene_dir, prefix + "pp_change.txt")
    with open(fn_pp_change) as f, open(fn_glist, 'w') as o:
        f.readline()
        # Transcript	Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        for i in f:
            line = i.split('\t', 4)
            transcript_id, logfc = line[0], line[3]
            if logfc != 'NA':
                print(transcript_id, file=o)
    
    # "perl heatmap.pl -w $out_dir -i $list -in1 $cond1_str -in2 $cond2_str -m $genome -n $tname";
    draw_heatmap_pp_change(n_gene_cols, analysis.out_dir, fn_glist, fls_ctrl=analysis.control_bed, fls_case=analysis.case_bed, ref_fls=analysis.ref, region_size=5000, bin_size=200, outname='heatmap')
    

    # modify here
    logger.info('debug exit')
    sys.exit(0)

def get_new_args(args, update_dict):
    args_dict = vars(args)
    args_dict.update(update_dict)
    return SimpleNamespace(**args_dict)

    
if __name__ == "__main__":
    args = getargs()
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
    
    args.o = os.path.realpath(args.o)
    args.organism = args.organism
    
    if args.in1 is not None:
        # logger.info('Processing input files')
        main()
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
                    out_dir = args.o or os.getcwd()
                    comp_str = control if case == 'null' else f'{case}_vs_{control}'
                    out_dir = os.path.join(out_dir, comp_str)
                    update_dict['o'] = out_dir
                    
                    args_new = get_new_args(args, update_dict)
                    # logger.info(vars(args_new))
                    logger.info(f'Processing {comp_str}')
                    main(args_new)
    else:
        logger.error("No input files provided, either use -in1 / in2 or provide a design table")
        sys.exit(1)

