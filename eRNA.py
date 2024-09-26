#! /data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12


import os, sys, time
import re
import traceback
import pickle
import logging
from types import SimpleNamespace

# global var
bin_dir = os.path.dirname(os.path.realpath(__file__))
pw_data = os.path.dirname(bin_dir)

# change here, just for testing
# pw_data = '/Users/files/work/bigfile/NRSA'

def red(s):
    return f'\u001b[31m{s}\u001b[0m'
def green(s):
    return f'\u001b[32m{s}\u001b[0m'



def getargs():
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('-in1', help="""required, read alignment files in bed (6 columns) or bam format for condition1, separated by space""", nargs='*')
    ps.add_argument('-in2', help="""read alignment files in bed (6 columns) or bam format for condition2, separated by space""", nargs='*')
    ps.add_argument('-design_table', '-design',  help="""Optional, desgin table in tsv format. 2 sections, first section is the sample information, 2 columns. col1=bed/bam full file path, col2 = group name. Second section is the comparison information, 2 columns, col1 should start with @@ used to identify this line as comparison definition. col1=group name used as case, col2 = group name used as control. e.g. @@case\\tcontrol. If only need to process a single group, use @@null\\tgroup_name""", nargs='?')

    ps.add_argument('-pwout', help="""work directory, should be the same of pause_PROseq.pl\'s output/work directory""", required=True)
    ps.add_argument('-organism', '-m', '-org',  help="""define the genome: hg19, hg38, mm10, dm3, dm6, ce10, or danRer10. default: hg19""", choices=['hg19', 'hg38', 'mm10', 'dm3', 'dm6', 'ce10', 'danRer10'], required=True)
    ps.add_argument('-gtf', help='Customized GTF file path, if not specified, will use the default one for the organism')

    # ps.add_argument('-overlap_frac', '-p', help="""percentage of overlap for calling annotated gene, float, default=0.2""", type=float, default=0.2)
    ps.add_argument('-cutoff', '-c', help="""distance cutoff for divergent transcripts for enhancer detection (bp, default: 400)""", type=int, default=400)
    ps.add_argument('-distance', '-d', help="""distance within which two eRNAs are merged (bp, default: 500)""", type=int, default=500)
    ps.add_argument('-le', '-long_ena', help="""length cutoff for long-eRNA identification (bp, default: 10000)""", type=int, default=10000)
    ps.add_argument('-filter', '-f', help="""whether to filter for enhancers: 0 (not filter) or 1 (filter), deflault: 1""", choices=[0, 1], default=1)
    ps.add_argument('-dtss', help="""if filter enhancers, the minimum distance from enhancer to TSS (Transcription Start Site) (bp, default: 2000)""", type=int, default=2000)
    ps.add_argument('-dtts', help="""if filter enhancers, the minimum distance from enhancer to TTS(Transcription Termination Site) (bp, default: 20000)""", type=int, default=20000)
    ps.add_argument('-wd', help="""distance within which associate active genes of enhancer (bp, default: 50000)""", type=int, default=50000)
    ps.add_argument('-pri', help="""if set, will prioritize enhancer. default is false""", action='store_true')
    ps.add_argument('-peak', help="""required if -pri set. ChIP-seq peak file for the regulator of interest in bed format""")
    ps.add_argument('-lk', help=""" valid when -pri is set for two conditions. Transcriptional changes to be used for enhancer prioritization, the value can be pp (changes in promoter-proximal region), gb (changes in gene body region), or pindex (changes in pausing index) (default: gb)""", choices=["pp", "gb", "pindex"] , default='gb')
    ps.add_argument('-dirction', help="""valid when -pri set  for two conditions. The expected simultaneous change of expression between enhancers and target genes. 1: the expression changes of enhancers and target genes are in the same direction; -1: the expression changes of enhancers and target genes are in the opposite direction; 0: the expression changes of enhancers and target genes are either in the same or in the opposite direction (default: 0)""", type=int, choices=[-1, 0, 1], default=0)
    ps.add_argument('-wt', '-weight', help="""valid when -pri set for two conditions. Weight to balance the impact of binding and function evidence. The higher the weight, the bigger impact does the binding evidence have (default: 0.5, i.e., binding and functional evidence have equal impact on enhancer prioritization)""", type=float, default=0.5)
    ps.add_argument('-fdr', '-cf', help="""valid when -pri set for two conditions. Cutoff to select significant transcriptional changes (default: FDR < 0.05). Use Foldchange instead if no FDR is generated (default: Foldchange > 1.5)""", type=float, default=0.05)
    ps.add_argument('-demo', help="""skip HOMER makeTagDirectory if already done""", action='store_true')
    ps.add_argument('-verbose', '-v', help="""verbose mode, show debug level logging info""", action='store_true')
    args = ps.parse_args()
    return args

args = getargs()

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


def updatelogger(logger, fn_log, terminal_level=None):
    handlers = list(logger.handlers)
    handler_name_d = {v.name: idx for idx, v in enumerate(handlers)}
    fh_console = handlers[handler_name_d['console']]
    formatter = fh_console.formatter
    valid_levels = {'DEBUG', 'INFO', 'WARNING', 'ERROR', 'FATAL'}
    if terminal_level:
        terminal_level = terminal_level.upper()
        if terminal_level not in valid_levels:
            logger.error(f'invalid logging level: {terminal_level} ')
            return logger
    if fn_log and 'file' not in handler_name_d:
        fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
        fh_file.setLevel('DEBUG')
        fh_file.setFormatter(formatter)
        fh_file.name = 'file'
        logger.addHandler(fh_file)
    if terminal_level is not None:
        fh_console.setLevel(terminal_level)
    return logger


logger = getlogger(logger_name='NRSA')


from utils import process_input, check_dependency, get_ref_erna, process_gtf, gtf_compare, get_other_region, get_enhancer, refine_chr, run_shell, sort_bed_like_file, change_enhancer, draw_signal, time_cost_util, parse_design_table, pause_longeRNA_main, prioritize_enhancer


def main(args):
    args_d = vars(args)
    demo = args.demo
    organism = args.organism
    # in1, in2, pwout, organism, overlap_frac, cutoff, le, filter(0,1), 
    # dtss, dtts, 
    # wd (distance within active gene), 
    # pri(0,1), peak (file), 
    # lk (pp, gb, pindex), dir (1, -1, 0), wt (weight), fdr (cutoff)
    pwout = args.pwout
    pwout = os.path.realpath(pwout)
    fn_gtf = args.gtf
    pwout_raw = args_d.get('pwout_raw', pwout)

    status = 0
    logger.debug('start running')
    ref_fls = get_ref_erna(args.organism)
    if ref_fls is None:
        logger.error("Error encountered while retrieving reference files")
        status = 1
        return 1


    folders = ['eRNA', 'intermediate']
    for ipw in folders:
        if not os.path.exists(os.path.join(pwout, ipw)):
            os.makedirs(os.path.join(pwout, ipw))
            
    dependency_status = check_dependency()
    if dependency_status:
        logger.error("Dependency check failed")
        sys.exit(1)
    pw_homer = f'{pwout}/intermediate/HOMER_tag'
    args.pw_homer = pw_homer

    in1 = process_input(pwout_raw, args.in1) # return = fn_lb, fn_bed
    if in1 is None:
        logger.error("Invalid input files provided for condition1")
        status = 1
    in2 = process_input(pwout_raw, args.in2) or []
    control_bed = in1
    case_bed = in2
    logger.debug(f'processed in1 = {in1}, processed in2 = {in2}')

    
    if status:
        logger.error("Error encountered while processing input files")
        return status
    
    # make homder tags
    fn_peak_txt, fn_peak_gtf = f'{pwout}/intermediate/transcript.txt', f'{pwout}/intermediate/transcript.gtf'
    if not (demo and os.path.exists(f'{pw_homer}/tagInfo.txt')):
        bed_list = ' '.join([_[1] for _ in in1 + in2])
        logger.info('makeTagDirectory...')
        s = time.time()
        cmd_homer = f'makeTagDirectory {pw_homer} -format bed -forceBED  {bed_list}'
        status = os.system(cmd_homer)
        if status:
            logger.error("Error encountered while creating tag directory")
            return status

        # find peaks
        logger.info(f'Find Peaks...')
        cmd_findpeaks = f'findPeaks {pw_homer} -style groseq -o {fn_peak_txt} -gtf {fn_peak_gtf}'
        status = os.system(cmd_findpeaks)
        if status:
            logger.error("Error encountered while running findPeaks")
            return status
    
        time_cost_util['homer'] = time.time() - s
    else:
        logger.warning(f'demo mode, Skip HOMER commands.')

    # cnofigs
    # overlap_frac, cutoff, le, filter(0,1), 
    # pri(0,1), peak (file), 
    # lk (pp, gb, pindex), dir (1, -1, 0), wt (weight), fdr (cutoff)
    # overlap_frac = args.overlap_frac
    
    window_active_genes = args.wd   #distance within which associate active genes of enhancer (bp, default: 50000)
    lcut = args.cutoff #cutoff of 5' distance;
    
    filter_tss = args.filter
    d_tss = args.dtss # if filter enhancers, the minimum distance from enhancer to TSS (Transcription Start Site) (bp, default: 2000)
    d_tts = args.dtts # if filter enhancers, the minimum distance from enhancer to TTS(Transcription Termination Site) (bp, default: 20000)
    thres_long_eRNA = args.le #   length cutoff for long-eRNA identification (bp, default: 10000), if eRNA > this value, will be identified as long eRNA
    cutoff = args.cutoff # distance cutoff for divergent transcripts for enhancer detection (bp, default: 400)
    thres_merge = args.distance # distance between two eRNAs for merging(bp, default: 500)
    
    config = {
        # 'overlap_frac': overlap_frac,
        'window_active_genes': window_active_genes,
        'lcut': lcut,
        'filter_tss': filter_tss,
        'd_tss': d_tss,
        'd_tts': d_tts,
        'thres_long_eRNA': thres_long_eRNA,
        'cutoff': cutoff,
        'thres_merge': thres_merge,
    }
    
    
    # compare the refseq gtf with peak_gtf
    logger.debug('process gtf')
    gtf_info, fn_tss, fn_tss_tts, err = process_gtf(ref_fls['gtf'], pwout)
    logger.debug('gtf info loaded')
    fn_enhancer_raw = f'{pwout}/eRNA/Enhancer.raw.txt'
    fn_enhancer = f'{pwout}/eRNA/Enhancer.txt'
    fn_enhancer_center = f'{pwout}/eRNA/Enhancer_center.txt'
    fn_pkl_other_region = f'{pwout}/intermediate/eRNA.other_region.pkl'
    fn_other_region = f'{pwout}/intermediate/eRNA.other_region.txt'
    fn_active_genes = f'{pwout}/intermediate/active_gene.txt'
    fn_active_tss = f'{pwout}/intermediate/active_tss.txt'
    fn_fantom, fn_association = [ref_fls[_] for _ in ['fantom', 'association']]
    fno_longerna = f'{pwout}/eRNA/long_eRNA.txt'
    fno_closest = f'{pwout}/intermediate/closest.txt'
    fn_4d = ref_fls['4d']
    fn_enhancer_short = f'{pwout}/intermediate/Enhancer_temp.bed'
    fn_count_enhancer = f'{pwout}/intermediate/count_enhancer.txt'

    
    if demo and os.path.exists(fn_pkl_other_region):
        logger.debug(f'loading previous gtf_compare results: {fn_pkl_other_region}')
        with open(fn_pkl_other_region, 'rb') as f:
            other_region = pickle.load(f)
    else:
        logger.info('gtf_compare')
        other_genes = gtf_compare(gtf_info, fn_peak_gtf) # the transcript without overlap with refseq gtf
        other_region = get_other_region(other_genes, fn_peak_txt) # is a list
        logger.debug('dump other region')
        with open(fn_pkl_other_region, 'wb') as o:
            pickle.dump(other_region, o)
    
    # other regions, [chr, s, e, strand], s and e are int
    # this list is exactly the same as perl code
    if not os.path.exists(fn_other_region):
        other_region = sorted(other_region, key=lambda x: (x[0], x[1]))
        tmp = ['\t'.join(map(str, i)) for i in other_region]
        with open(fn_other_region, 'w') as o:
            print('\n'.join(tmp), file=o)
    logger.debug(f'other regions, n = {len(other_region)}')

    # load active genes
    # 2 columns, col2 = ts, col2 = gn
    ts_to_gn = {}
    active_genes = set()
    with open(fn_active_genes) as f:
        for i in f:
            ts, gn = i[:-1].split('\t')
            active_genes.add(ts)
            ts_to_gn[ts] = gn

    active_tss = {}
    logger.debug(f'getting active tss from {fn_tss}')
    with open(fn_tss) as f, open(fn_active_tss, 'w') as o:
        # chr1	11874	11874	NR_046018.2	+
        for i in f:
            line = i.split('\t')
            chr_, tss, ts = line[0], line[1], line[3]
            chr_ = refine_chr(chr_)
            
            if ts in active_genes:
                p1, p2 = i.split('\t', 1)
                o.write(f'{chr_}\t{p2}')
                
                # o.write(i)
                gn = ts_to_gn.get(ts, ts)
                active_tss.setdefault(chr_, []).append([int(tss), gn])
    
    # sort
    status = sort_bed_like_file(fn_active_tss)
    if status:
        return 1
    
    # find the enhancer region


    if not (demo and os.path.exists(fn_enhancer)):
        logger.info(f'Finding enhancer region')
        enh_out, lerna_out = get_enhancer(other_region, fn_fantom, fn_association, fn_tss_tts, lcut=lcut, thres_long_eRNA=thres_long_eRNA, distance_for_merge=thres_merge, filter=filter_tss)
        logger.debug(f'enh_out = {len(enh_out)}, lerna count = {len(lerna_out)}')
        
        # enh_out = [chr, start, end, center_list, fantom_list, asso_list]
        if enh_out:
            with open(fn_enhancer_raw, 'w') as o:
                print('\n'.join(enh_out), file=o)
        else:
            logger.warning(f'No enhancer found')
            return 1

        status = sort_bed_like_file(fn_enhancer_raw)
        if status:
            return 1
        
        # save long_eRNA
        # "chr","start","end","strand"
        if lerna_out:
            with open(fno_longerna, 'w') as o:
                print('\t'.join(["chr","start","end","strand"]), file=o)
                for i in lerna_out:
                    print('\t'.join(map(str, i)), file=o)

        # find closest gene for enhancer
        logger.info(f'Finding closest TSS for enhancer')
        cmd = f'bedtools closest -a {fn_enhancer_raw} -b {fn_active_tss} -d > {fno_closest}'
        logger.debug(cmd)
        status = run_shell(cmd, echo=True)
        if status:
            logger.error(f'fail to run bedtools closest')
            return 1

        # process closest result
        # last column will be the distance, if overlap, will be 0
        logger.debug('processing bedtools closest results')
        res_closest = {}  # k = enhancer_info, v = {'closeclosest_gene_gene': [],  distance: distance}
        center_str = {} # key = each center_point, v = chr, center point, fantom
        center_enhancer = {}  # k = each center point, v = whole enhancer info
        with open(fno_closest) as f:
            for i in f:
                line = i[:-1].rsplit('\t', 6) # first 6 cols = enhancer info, then tss info [chr(col7), s, s , ts(col10), strand (col11)],  col12 = distance
                distance = int(line[-1])
                enhancer_info = line[0]
                ts = line[-3]
                
                if distance == -1:
                    continue
                res_closest.setdefault(enhancer_info, {'closest_gene': [], 'distance': distance})
                gn =  ts_to_gn.get(ts, ts)
                res_closest[enhancer_info]['closest_gene'].append(gn)
                res_closest[enhancer_info]['distance'] = distance
                chr_, _, _, center_list, fantom_list, asso_list = enhancer_info.split('\t')
                for icenter, ifantom in zip(center_list.split(','), fantom_list.split(',')):
                    k_center = f'{chr_}:{icenter}'
                    if k_center not in center_str:
                        center_str[k_center] = f'{chr_}\t{icenter}\t{ifantom}'
                        center_enhancer[k_center] = enhancer_info
        logger.debug(f'enhancer_closest gene, n = {len(res_closest)}, center_count = {len(center_enhancer)}')

        d_4d = {}
        if fn_4d is not None:
            # chr1	557489	560146	chr5	134284878	134293544	PCBD2;	MCF7	ChIA-PET	19890323
            with open(fn_4d) as f:
                # no header
                for i in f:
                    line = i[:-1].split('\t')
                    chr_ = refine_chr(line[0])
                    d_4d.setdefault(chr_, []).append([line[_] for _ in [1, 2, 6, 7, 8, 9]]) #  s, e, gn, cell, method, pmid 
        
        logger.info(f'Finding associated genes within 50kb and 4DGenome')
        tss50k = {} # window_active_genes = 50k
        res_4d = {}
        for enhancer_info in sorted(res_closest):
            # closest_info = res_closest[enhancer_info]
            chr_, start, end, center_list, fantom_list, asso_list = enhancer_info.split('\t')
            start, end = int(start), int(end)
            for tss_pos, gn in active_tss[chr_]:
                if (tss_pos > end and tss_pos - end < window_active_genes) or (tss_pos < start and start - tss_pos < window_active_genes):
                    tss50k.setdefault(enhancer_info, set()).add(gn)
                    
            if chr_ in d_4d:
                for info_4d in d_4d[chr_]:
                    # s, e, gn, cell, method, pmid 
                    s_4d, e_4d, gn_list, cell, method, pmid = info_4d
                    s_4d, e_4d = int(s_4d), int(e_4d)
                    if (start <= s_4d <= end) or (start <= e_4d <= end) or (s_4d <= start <= e_4d) or (s_4d <= end <= e_4d):  # there are overlap
                        ires = res_4d.setdefault(enhancer_info, {'gn': set(), 'cell': set(), 'method': set(), 'pmid': set()})
                        ires['gn'] |= set(gn_list.strip(';').split(';'))
                        ires['cell'].add(cell)
                        ires['method'].add(method)
                        ires['pmid'].add(pmid)
                        
        # save to Enhancer.txt again
        enhancer_id_map = {}
        enhancer_id = 1
        with open(fn_enhancer, 'w') as o:
            print("Enhancer_ID\tchr\tstart\tend\tcenter\tFANTOM5\tassociated_gene-FANTOM5\tassociated_gene-50kb\tassociated_gene-4DGenome\tCell_Tissue\tDetection_method\tPMID\tclosest_gene\tdistance", file=o)
            for enhancer_info in sorted(res_closest):
                col_tss50k = tss50k.get(enhancer_info, ['NA'])
                col_tss50k = ','.join(sorted(col_tss50k))
                
                info_4d = res_4d.get(enhancer_info, None)
                if info_4d:
                    cols_4d = '\t'.join([','.join(sorted(info_4d[_])) for _ in ['gn', 'cell', 'method', 'pmid']])
                else:
                    cols_4d = 'NA\tNA\tNA\tNA'
                
                closest_info = res_closest.get(enhancer_info)
                if closest_info:
                    col_closest = closest_info['closest_gene']
                    col_distance = str(closest_info['distance'])
                    if col_closest:
                        col_closest = ','.join(sorted(col_closest))
                    else:
                        col_closest = 'NA'
                else:
                    col_closest,col_distance = 'NA', '-1'
                print('\t'.join([str(enhancer_id), enhancer_info, col_tss50k, cols_4d, col_closest, col_distance]), file=o)
                enhancer_id_map[enhancer_info] = str(enhancer_id)
                enhancer_id += 1
                
        # save center_str
        # center_str[k_center] = f'{chr_}\t{icenter}\t{ifantom}'
        center_out = []
        with open(fn_enhancer_center, 'w') as o:
            print('\t'.join(["chr","position","FONTOM5","Enhancer_ID"]), file=o)
            for k_center, v in center_str.items():
                enhancer_info = center_enhancer[k_center]
                enhancer_id = enhancer_id_map[enhancer_info]
                print(f'{v}\t{enhancer_id}', file=o)
        
        # currently, the results excatly match with version1
        os.unlink(fn_enhancer_raw)
        # logger.info(f'uncomment above')
    else:
        logger.warning(f'debug mode, skip get enhancer steps')
        lerna_out = 1 if os.path.exists(fno_longerna) else 0

        
        
    # enhancer region short
    k_enhancer_map = {}
    
    # the chrom patter is different for this file and the raw bed file,
    # to run bedtools coverage, need to convert the refined chr back
    # 4104-AW-1_sorted_rRNArm-F4q10.chr_map.pkl
    
    logger.debug('loading the chr_map to previous raw bed files')
    fn_lb = in1[0][0] # use the first file as the fn_lb
    fn_chr_map = f'{pwout_raw}/bed/{fn_lb}.chr_map.pkl'
    with open(fn_chr_map, 'rb') as f:
        chr_map_to_raw_bed = pickle.load(f)
    
    chr_not_converted = set()
    with open(fn_enhancer) as f, open(fn_enhancer_short, 'w') as o:
        f.readline()
        for i in f:
            line = i.split('\t', 4)
            # Enhancer_ID	chr	start	end	center
            k_enhancer = '\t'.join(line[2:4])
            chr_ = line[1]
            if chr_ in chr_map_to_raw_bed:
                chr_raw = chr_map_to_raw_bed[chr_]
            else:
                chr_raw = chr_
                chr_not_converted.add(chr_)
            enhancer_id = line[0]
            k_enhancer = f'{chr_raw}\t{k_enhancer}'
            print(k_enhancer, file=o)
            k_enhancer_map[k_enhancer] = enhancer_id
    
    if len(chr_not_converted) > 0:
        logger.debug(f'{len(chr_not_converted)} lines , the chrom not converted back to raw:  {sorted(chr_not_converted)}, chr_map = {chr_map_to_raw_bed}')
    status = sort_bed_like_file(fn_enhancer_short)
    if status:
        return 1

    # logger.warning(f'enhancer count = {len(k_enhancer_map)}')

    logger.info('Detecting Enhancer change...')
    n_ctrl, n_case = len(in1), len(in2)
    factors_d = {}  # key = fn_lb, v = nf value
    fn_nf = f'{pwout}/intermediate/nf.txt'
    if not os.path.exists(fn_nf):
        logger.error(f'normalization factor file {fn_nf} not found, please run pause_PROSeq first')
        return 1
    with open(fn_nf) as f:
        f.readline()  # header
        tmp = [i[:-1].split('\t') for i in f if i.strip()]
        factors_d = dict([(fn_lb, float(v)) for fn_lb, v in tmp])
    
    fls_in = in1 + in2
    enhancer_count = {}  # key = enhancer info, chr, start, end, v = list of number,  length is same as sample number
    n_sam = n_ctrl + n_case
    sam_idx = -1
    logger.info(f'Getting enhancer read count...')
    for fn_lb, fn_bed in fls_in:
        sam_idx += 1
        fn_coverage_res = f'{pwout}/intermediate/enhancer_cov.{fn_lb}.bed'
        cmd = f'bedtools coverage -a {fn_enhancer_short} -b {fn_bed} -counts -sorted > {fn_coverage_res}'
        if demo and os.path.exists(fn_coverage_res):
            logger.debug(f'skip bedtools coverage due to demo mode')
        else:
            logger.debug(cmd)
            retcode = run_shell(cmd, echo=True)
            if retcode:
                return 1
        # process
        if os.path.getsize(fn_coverage_res) < 10:
            logger.warning(f'Empty bedtools coverage file, cmd = \n{cmd}')
            continue
        with open(fn_coverage_res) as f:
            for i in f:
                try:
                    k_enhancer, ict = i[:-1].rsplit('\t', 1)
                except:
                    logger.error(f'invalid coverage line: {i[:-1]}')
                    sys.exit(1)
                if k_enhancer not in enhancer_count:
                    enhancer_count[k_enhancer] = ['0' for _ in range(n_sam)]
                enhancer_count[k_enhancer][sam_idx] = ict

    if len(enhancer_count) == 0:
        logger.error(f'no Enhancer count found')
        return 1
    # logger.warning(len(enhancer_count))
    # dump the count
    with open(fn_count_enhancer, 'w') as o:
        print(f'Enhancer_ID\tchr\tstart\tend\t' + '\t'.join([_[0] for _ in fls_in]), file=o)
        for k_enhancer, v in enhancer_count.items():
            enhancer_id = k_enhancer_map[k_enhancer]
            print(f'{enhancer_id}\t{k_enhancer}\t' + '\t'.join(v), file=o)
    os.system(f'sort -k 1,1n {fn_count_enhancer} > {fn_count_enhancer}.tmp && mv {fn_count_enhancer}.tmp {fn_count_enhancer}')
    
    # change_enhancer
    logger.info(f'change_enhancer')
    sam_ctrl = [_[0] for _ in in1]
    sam_case = [_[0] for _ in in2]
    change_enhancer(pwout, fn_count_enhancer, factors_d, n_ctrl, n_case, sam_ctrl, sam_case, flag=1)
    
    # draw signal
    logger.info('Drawing signal around enhancer center...')
    signal_type = 'p'
    draw_signal(pwout, fn_enhancer_center, in1, in2, chr_map_to_raw_bed, signal_type=signal_type, demo=demo)
    
    # pause longeRNA
    
    if lerna_out:
        logger.info(f'Detecting long eRNA change...')
        # perl pause_longeRNA.pl -o $out_dir -a $out3 -in1 $cond1_str -in2 $cond2_str -m $genome";
        args = SimpleNamespace()
        args.in1 = in1
        args.in2 = in2
        args.pwout = pwout
        args.organism = organism
        args.pw_bed = f'{pwout_raw}/bed'
        args.pwout_raw = pwout_raw
        args.gtf = fno_longerna
        args.f1 = None
        args.f2 = None  # the f1 and f2 are get from nf.txt
        args.skip_get_mapped_reads = 1 if demo else 0
        logger.info('running pause longeRNA')
        pause_longeRNA_main(args)

    
    # ps.add_argument('-pri', help="""if set, will prioritize enhancer. default is false""", action='store_true')
    # ps.add_argument('-peak', help="""required if -pri set. ChIP-seq peak file for the regulator of interest in bed format""")
    # ps.add_argument('-lk', help=""" valid when -pri is set for two conditions. Transcriptional changes to be used for enhancer prioritization, the value can be pp (changes in promoter-proximal region), gb (changes in gene body region), or pindex (changes in pausing index) (default: gb)""", choices=["pp", "gb", "pindex"] , default='gb')
    # ps.add_argument('-dirction', help="""valid when -pri set  for two conditions. The expected simultaneous change of expression between enhancers and target genes. 1: the expression changes of enhancers and target genes are in the same direction; -1: the expression changes of enhancers and target genes are in the opposite direction; 0: the expression changes of enhancers and target genes are either in the same or in the opposite direction (default: 0)""", type=int, choices=[-1, 0, 1], default=0)
    # ps.add_argument('-wt', '-weight', help="""valid when -pri set for two conditions. Weight to balance the impact of binding and function evidence. The higher the weight, the bigger impact does the binding evidence have (default: 0.5, i.e., binding and functional evidence have equal impact on enhancer prioritization)""", type=float, default=0.5)
    # ps.add_argument('-fdr', '-cf', help="""valid when -pri set for two conditions. Cutoff to select significant transcriptional changes (default: FDR < 0.05). Use Foldchange instead if no FDR is generated (default: Foldchange > 1.5)""", type=float, default=0.05)

    
    if args.pri:
        fn_peak = args.peak
        if not os.path.exists(fn_peak):
            logger.error(f'ChIP-seq peak file not exist: {fn_peak}')
            return 1
        logger.info(f'g@Now running prioritizing...')
        prioritize_enhancer(pwout, fn_peak, n_ctrl, n_case, args.direction, args.wt, args.fdr, args.lk)

   
    return 0

if __name__ == "__main__":
    args = getargs()
    pwout = args.pwout
    os.makedirs(pwout, exist_ok=True)
    fn_log = f'{pwout}/eRNA.run.log'
    fn_log_base = os.path.basename(fn_log)
    terminal_level = 'DEBUG' if args.verbose else None
    logger = updatelogger(logger, fn_log, terminal_level=terminal_level)

    if os.path.exists(fn_log_base):
        os.unlink(fn_log_base)
    if os.path.exists(fn_log):
        os.symlink(fn_log, fn_log_base)

    logger.debug(f'working in {os.getcwd()}')
    logger.debug(f'inpu args = {vars(args)}')


    pwout_raw = os.path.realpath(args.pwout)
    logger.debug(f'pw_out_raw = {pwout_raw}')
    args.pwout = pwout_raw
    args.pwout_raw = pwout_raw
    args.pw_bed = f'{pwout_raw}/bed'
    if not os.path.exists(args.pw_bed):
        os.makedirs(args.pw_bed, exist_ok=True)

    arg_list = parse_design_table(args)
    retcode = 0
    for comp_str, iargs in arg_list:
        if comp_str:
            logger.info(f'g@now running {comp_str}')
        logger.debug(vars(iargs))
        try:
            retcode = main(iargs) or retcode
        except:
            e = traceback.format_exc()
            logger.error(f'Error found during running, please check the log file for more information: {fn_log}')
            logger.debug(e)
            sys.exit(1)

    if not retcode:
        logger.debug(f'g@script finished without error')
    else:
        logger.debug(f'error encountered during running')

    sys.exit(retcode)
