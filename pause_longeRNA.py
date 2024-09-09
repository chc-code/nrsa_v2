#! /data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12
"""
process PROseq data
"""
import time
s = time.time()
import sys, os, re
import logging
import pickle
import json
import traceback
import gc
sys.dont_write_bytecode = True


from utils import check_dependency, run_shell, process_input, get_seqence_from_fa, build_idx_for_fa, get_ref, get_peak,  change_pp_gb, change_pindex, draw_box_plot, draw_heatmap_pindex, draw_heatmap_pp_change, get_alternative_isoform_across_conditions, get_FDR_per_sample, pre_count_for_bed, add_value_to_gtf, get_ref_erna

from utils import Analysis, process_bed_files


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

logger = getlogger('NRSA')


def main(args):
    """
    args: an object with attributes equiv to the argparse object, must specify 
    """

    # check if the attributes are present
    defined_attrs = vars(args)
    required_attrs = {'in1', 'pwout', 'organism', 'pw_bed', 'gtf'}
    # f1, f2 = normalization factors
    
    optional_attrs = {
        'in2': None,
        'fa': None,
        'up': 500,
        'down': 500,
        'gb': 1000,
        'min_gene_len': 1000,
        'window': 50,
        'step': 5,
        'ignore': False,
        'skip_get_mapped_reads': False,
    }
    pwout = args.pwout
    pw_bed = args.pw_bed
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
    # check dependencies
    
    skip_get_mapped_reads = args.skip_get_mapped_reads
    dependency_status = check_dependency()
    if dependency_status:
        logger.error("Dependency check failed")
        sys.exit(1)

    analysis = Analysis(args, is_long_eRNA=True, raw_input=False, skip_get_mapped_reads=skip_get_mapped_reads) # raw_input set as False, will not run the process_input function
    if analysis.status:
        logger.error("Exit now")
        sys.exit(1)
    
    logger.debug(vars(args))

    rep1 = len(analysis.control_bed)
    rep2 = len(analysis.case_bed) if analysis.case_bed else 0
    window_size = analysis.config['window_size']
    
    # get the normalization factors
    factor_flag = 1
    fn_nf = f'{pwout}/intermediate/nf.txt'
    nf_d = {}
    try:
        with open(fn_nf) as f:
            f.readline() # header
            for i in f:
                lb, v = i[:-1].split('\t')
                nf_d[lb] = float(v)
    except:
        e = traceback.format_exc()
        logger.error(f'fail to get normalization factor, please make sure pause_PROseq is executed')
        logger.error(e)
        sys.exit(1)
    
    factor1 = [nf_d[lb] for lb, _ in analysis.control_bed]
    factor2 = [nf_d[lb] for lb, _ in analysis.case_bed]

    pro_up, pro_down, gb_down_distance, min_gene_len = [analysis.config[_] for _ in ['pro_up', 'pro_down', 'gb_start', 'min_gene_len']]

    # process the GTF file
    fn_gtf = analysis.ref['gtf']
    logger.debug(f"Processing longeRNA gtf file: {fn_gtf}")
    
    gtf_info = {}
    with open(fn_gtf) as f:
        f.readline()  # header
        # chr     start   end     strand
        # chr1    27365605        27376193        -
        # chr1    99872187        99904465        -
        # chr1    103653695       103664771       -
        # {'chr': chr_, 'strand': strand, 'gene_name': gene_name, 'start': start, 'end': end}
        for i in f:
            chr_, s, e, strand = i[:-1].split('\t')
            ts = f'{chr_}:{s}-{e}:{strand}'
            ires = {'chr': chr_, 'strand': strand, 'start': int(s), 'end': int(e)}
            gtf_info[ts] = add_value_to_gtf(ires, pro_up, pro_down, gb_down_distance) 
    
    # prepare fasta file
    fn_fa = analysis.ref['fa']
    fa_idx = build_idx_for_fa(fn_fa)
    fh_fa = open(fn_fa, 'r')

    fls = analysis.input_fls # element is [fn_lb, fn_bed]
    sam_order = [_[0] for _ in fls]
    n_gene_cols = analysis.n_gene_cols

    def close_fh():
        # close file handle
        for fn_lb, fn_bed in fls:
            analysis.out_fls['bed_peaks'][fn_lb]['fh'].close()
        fh_fa.close()

    fn_count_pp_gb = analysis.out_fls['count_pp_gb']
    if skip_get_mapped_reads and os.path.exists(fn_count_pp_gb):
        logger.warning(f'demo mode, skip getting pp_gb count, reuse prev res')
        close_fh()
    else:
        reuse_pre_count = True # modify here
        logger.info(f'Getting pp_gb count')
        pp_str, gb_str = process_bed_files(analysis, fls, gtf_info, fa_idx, fh_fa, reuse_pre_count=reuse_pre_count)
        close_fh()
    
        # count_pp_gb
        logger.info('Building count_pp_gb')
        prefix = 'longeRNA-' if analysis.longerna else ''
        header = analysis.gene_cols
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
                row += pp_str[transcript_id] + gb_str[transcript_id]
                print('\t'.join(row), file=o)

        
    # calculate FDR
    close_fh()
    logger.info('Calculating FDR')
    pause_index_str = {}
    for fn_lb, fn_bed in fls:
        fn_peak = analysis.out_fls['bed_peaks'][fn_lb]['fn']
        fn_fdr = analysis.out_fls['fdr'][fn_lb]['fn']
        pause_index_str = get_FDR_per_sample(fn_lb, fn_peak, fn_fdr,  pause_index_str)

    pause_index_str_new = {}
    for transcript_id, v1 in pause_index_str.items():
        v2 = []
        for isam in sam_order:
            v3 = v1.get(isam, ['NA', 'NA', 'NA'])
            v2 += v3
        pause_index_str_new[transcript_id] = v2
    pause_index_str = pause_index_str_new
    

    fno_prefix = 'longeRNA-' if analysis.longerna else ''

    # modify here
    logger.info('Change_pp_gb')
    # logger.warning('modify here')
    change_pp_gb(n_gene_cols, fn_count_pp_gb, analysis.out_dir, rep1, rep2, window_size, factor1=factor1, factor2=factor2, factor_flag=factor_flag, islongerna=True)
    
    logger.debug('Dump pindex.txt')
    header_extra = []
    for fn_lb, fn_bed in fls:
        header_extra += [f'{fn_lb}-pindex', f'{fn_lb}-pvalue', f'{fn_lb}-FDR']
    header = analysis.gene_cols + header_extra
    fno =  f'{pwout}/eRNA/{fno_prefix}pindex.txt'
    fno_pindex_change = f'{pwout}/eRNA/{fno_prefix}pindex_change.txt'
    with open(fno, 'w') as o:
        print('\t'.join(header), file=o)
        for transcript_id, data in pause_index_str.items():
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gtf_info[transcript_id]['gene_name'])
            row += data
            print('\t'.join(map(str, row)), file=o)

    if rep2 > 0:
        # change_pindex
        logger.info('Running change_pindex')
        change_pindex(fno_prefix, n_gene_cols, fn_count_pp_gb, fno_pindex_change, rep1, rep2, window_size, factor1=factor1, factor2=factor2, factor_flag=factor_flag)
