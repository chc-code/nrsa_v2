#! /data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12


import os, sys
import re

# global var
bin_dir = os.path.dirname(os.path.realpath(__file__))
pw_data = os.path.dirname(bin_dir)

# change here, just for testing
# pw_data = '/Users/files/work/bigfile/NRSA'

def red(s):
    return f'\u001b[31m{s}\u001b[0m'
def green(s):
    return f'\u001b[32m{s}\u001b[0m'


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

logger = getlogger('eRNA.run.log', 'NRSA')

def getarg():
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('-in1', help="""required, read alignment files in bed (6 columns) or bam format for condition1, separated by space""", nargs='+', required=True)
    ps.add_argument('-in2', help="""read alignment files in bed (6 columns) or bam format for condition2, separated by space""", nargs='*')
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
    args = ps.parse_args()
    return args


def main():
    from utils import process_input, check_dependency, get_ref_erna, process_gtf, gtf_compare, get_other_region, get_enhancer, refine_chr
    from types import SimpleNamespace
    
    
    args = getarg()
    args_d = vars(args)
    # in1, in2, pwout, organism, overlap_frac, cutoff, le, filter(0,1), 
    # dtss, dtts, 
    # wd (distance within active gene), 
    # pri(0,1), peak (file), 
    # lk (pp, gb, pindex), dir (1, -1, 0), wt (weight), fdr (cutoff)
    analysis = SimpleNamespace()
    analysis.pwout = pwout = args.pwout
    fn_gtf = args.gtf
    analysis.pwout_raw = args_d.get('pwout_raw', pwout)

    status = 0
    logger.debug('start running')
    ref_fls = get_ref_erna(args.organism)
    if ref_fls is None:
        logger.error("Error encountered while retrieving reference files")
        status = 1
        return
    analysis.ref_fls = ref_fls


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

    in1 = process_input(analysis.pwout_raw, args.in1) # return = fn_lb, fn_bed
    if in1 is None:
        logger.error("Invalid input files provided for condition1")
        status = 1
    in2 = process_input(analysis.pwout_raw, args.in2) or []
    analysis.control_bed = in1
    analysis.case_bed = in2

    
    if status:
        logger.error("Error encountered while processing input files")
        return status
    
    # make homder tags
    bed_list = ' '.join([_[1] for _ in in1 + in2])
    logger.info('makeTagDirectory...')
    cmd_homer = f'makeTagDirectory {pw_homer} -format bed -forceBED  {bed_list}'
    status = os.system(cmd_homer)
    if status:
        logger.error("Error encountered while creating tag directory")
        return status
    
    # find peaks
    logger.info(f'Find Peaks...')
    fn_peak_gtf, fn_peak_txt = f'{pwout}/intermediate/trascript.txt', f'{pwout}/intermediate/trascript.gtf'
    cmd_findpeaks = f'findPeaks {pw_homer} -style groseq -o {fn_peak_txt} -gtf {fn_peak_gtf}'
    status = os.system(cmd_findpeaks)
    if status:
        logger.error("Error encountered while running findPeaks")
        return status

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
    
    analysis.config = {
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
    gtf_info, fn_tss, fn_tss_tts, err = process_gtf(ref_fls['gtf'], pwout)
    other_genes = gtf_compare(gtf_info, fn_peak_gtf) # the transcript without overlap with refseq gtf
    other_region = get_other_region(other_genes, fn_peak_txt)

    # load active genes
    # 2 columns, col2 = ts, col2 = gn
    logger.debug('loading active TSS')
    fn_active_genes = f'{pwout}/intermediate/active_gene.txt'
    fn_active_tss = f'{pwout}/intermediate/active_tss.txt'
    with open(fn_active_genes) as f:
        active_genes = {_.split('\t')[0] for _ in f if _.strip()}

    active_tss = {}
    with open(fn_tss) as f, open(fn_active_tss, 'w') as o:
        # chr1	11874	11874	NR_046018.2	+
        for i in f:
            line = i.split('\t')
            chr_, tss, ts = line[0], line[1], line[3]
            if ts in active_genes:
                o.write(i)
                gn = gtf_info.get(ts, {}).get('gene_name', ts)
                active_tss.setdefault(chr_, []).append([tss, gn])
    
    
    # find the enhancer region
    logger.info(f'Finding enhancer region')
    fn_fantom, fn_association = [ref_fls[_] for _ in ['fantom', 'association']]
    enh_out, lerna_out = get_enhancer(other_region, fn_fantom, fn_association, fn_tss_tts, lcut=lcut, thres_long_eRNA=thres_long_eRNA, distance_for_merge=thres_merge, filter=filter_tss)
    
    # enh_out = [chr, start, end, center_list, fantom_list, asso_list]
    fn_enhancer = f'{pwout}/eRNA/Enhancer.txt'
    if enh_out:
        with open(fn_enhancer, 'w') as o:
            print('\n'.join(enh_out), file=o)
    else:
        logger.warning(f'No enhancer found')
        return 1


    # find closest gene for enhancer
    logger.info(f'Finding closest TSS for enhancer')
    fno_closest = f'{pwout}/intermediate/closest.txt'
    cmd = f'bedtools closest -a {fn_enhancer} -b {fn_active_tss} -d > $fno_closest'
    status = os.system(cmd)
    if status:
        logger.error(f'fail to run bedtools closest')
        return 1

    # process closest result
    # last column will be the distance, if overlap, will be 0
    logger.debug('processing bedtools closest results')
    res_closest = {}  # k = enhancer_info, v = {'close_gene': [], }
    center_str = {} # key = each center_point, v = chr, center point, fantom
    center_enhancer = {}  # k = each center point, v = whole enhancer info
    with open(fno_closest) as f:
        for i in f:
            line = i[:-1].rsplit('\t', 6) # first 6 cols = enhancer info, then tss info [chr(col7), s, s , ts(col10), strand (col11)],  col12 = distance
            distance = int(line[-1])
            enhancer_info = line[0]
            ts = line[-3]
            res_closest.setdefault(enhancer_info, {'close_gene': [], 'distance': distance})
            if distance == -1:
                continue
            gn = gtf_info.get(ts, {}).get('gene_name', ts)
            res_closest[enhancer_info]['closest_gene'].append(gn)
            res_closest[enhancer_info]['distance'] = distance
            chr_, _, _, center_list, fantom_list, asso_list = enhancer_info.split('\t')
            for icenter, ifantom in zip(center_lis.split(','), fantom_lis.split(',')):
                k_center = f'{chr_}:{icenter}'
                if k_center not in center_str:
                    center_str[k_center] = f'{chr_}\t{icenter}\t{ifantom}'
                    center_enhancer[k_center] = enhancer_info
            

    fn_4d = ref_fls['4d']
    d_4d = {}
    if fn_4d is not None:
        # chr1	557489	560146	chr5	134284878	134293544	PCBD2;	MCF7	ChIA-PET	19890323
        with open(fn_4d) as f:
            # no header
            for i in f:
                line = i[:-1].split('\t')
                chr_ = refine_chr(line[0])
                d_4d.setdefault(chr_, []).append([line[_] for _ in [1, 2, 6, 7, 8, 9]]) #  s, e, gn, cell, method, pmid 
    
    logger.info(f'Fiinding associated genes within 50kb and 4DGenome')
    tss50k = {} # window_active_genes = 50k
    res_4d = {}
    for enhancer_info in sorted(res_closest):
        # closest_info = res_closest[enhancer_info]
        chr_, start, end, center_list, fantom_list, asso_list = enhancer_info.split('\t')
        for tss_pos, gn in active_tss[chr_]:
            if (tss_pos > end and tss_pos - end < window_active_genes) or (tss_pos < strand and start - tss_pos < window_active_genes):
                tss50k.setdefault(enhancer_info, set()).add(gn)
                
        if chr_ in d_4d:
            for info_4d in d_4d[chr_]:
                # s, e, gn, cell, method, pmid 
                s_4d, e_4d, gn_list, cell, method, pmid = info_4d
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
                cols_4d = '\t'.join([info_4d[_] for _ in ['gn', 'cell', 'method', 'pmid']])
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
            print('\t'.join([enhancer_id, enhancer_info, col_tss_50k, cols_4d, col_closest, col_distance]))
            enhancer_id_map[enhancer_info] = str(enhancer_id)
            enhancer_id += 1
            
    # save center_str
    # center_str[k_center] = f'{chr_}\t{icenter}\t{ifantom}'
    fn_center = f'{pwout}/eRNA/Enhancer_cener.txt'
    center_out = []
    with open(fn_center, 'w') as o:
        print('\t'.join(["chr","position","FONTOM5","Enhancer_ID"]), file=o)
        for k_center, v in center_str.items():
            enhancer_info = center_enhancer[k_center]
            enhancer_id = enhancer_id_map[enhancer_info]
            print(f'{v}\t{enhancer_id}', file=o)
    
    
    # save long_eRNA
    # "chr","start","end","strand"
    fno_longerna = f'{pwout}/eRNA/long_eRNA.txt'
    if lerna_out:
        with open(fno_longerna, 'w') as o:
            print('\t'.join(["chr","start","end","strand"]), file=o)
            for i in lerna_out:
                print('\t'.join(map(str, i)), file=o)
    
    
    
    
    
    # temp_str = ""
    # temp_id = {}
    # line = 0
    # with open(out1, 'r') as in_file:
    #     for line_num, line in enumerate(in_file):
    #         if line_num == 0:
    #             continue
    #         temp = line.strip().split('\t')
    #         temp_str += f"{temp[1]}\t{temp[2]}\t{temp[3]}\n"
    #         temp_id[f"{temp[1]}-{temp[2]}-{temp[3]}"] = temp[0]

    # out4 = os.path.join(inter_dir, "Enhancer_temp.bed")
    # with open(out4, 'w') as out_file:
    #     out_file.write(temp_str)

    # print("Detecting Enhancer change...")
    # if case_bed:
    #     factor1 = []
    #     factor2 = []
    #     nffile = os.path.join(inter_dir, "nf.txt")
    #     factor = {}
    #     with open(nffile, 'r') as in_file:
    #         for line in in_file:
    #             temp = line.strip().split('\t')
    #             factor[temp[0]] = temp[1]

    #     count_str = {}
    #     rep_str = ""

    #     for rep in cond1:
    #         temp1 = rep.rindex("/")
    #         temp2 = rep.rindex(".")
    #         temp3 = rep[temp1 + 1:temp2]

    #         if temp3 in factor:
    #             factor1.append(factor[temp3])
    #         else:
    #             print("input file name doesn't match the file name for pausing analysis!")
    #             exit()

    #         rep_str += f"\t{temp3}"
    #         cov_temp = os.path.join(inter_dir, "cov_temp.bed")
    #         cov_cmd1 = f"bedtools coverage -a {out4} -b {rep} > {cov_temp}"
    #         subprocess.run(cov_cmd1, shell=True, check=True)
    #         with open(cov_temp, 'r') as in_file:
    #             for line in in_file:
    #                 temp = line.strip().split()
    #                 key = f"{temp[0]}-{temp[1]}-{temp[2]}"
    #                 count_str.setdefault(key, []).append(temp[3])

    #     for rep in cond2:
    #         temp1 = rep.rindex("/")
    #         temp2 = rep.rindex(".")
    #         temp3 = rep[temp1 + 1:temp2]

    #         if temp3 in factor:
    #             factor2.append(factor[temp3])
    #         else:
    #             print("input file name doesn't match the file name for pausing analysis!")
    #             exit()

    #         rep_str += f"\t{temp3}"
    #         cov_temp = os.path.join(inter_dir, "cov_temp.bed")
    #         cov_cmd1 = f"bedtools coverage -a {out4} -b {rep} > {cov_temp}"
    #         subprocess.run(cov_cmd1, shell=True, check=True)
    #         with open(cov_temp, 'r') as in_file:
    #             for line in in_file:
    #                 temp = line.strip().split()
    #                 key = f"{temp[0]}-{temp[1]}-{temp[2]}"
    #                 count_str.setdefault(key, []).append(temp[3])

    #     temp_count = 0
    #     temp_str2 = "Enhancer_ID\tchr\tstart\tend"
    #     for rep in cond1:
    #         temp_str2 += f"\t{factor1[temp_count]}"
    #         temp_count += 1
    #     for rep in cond2:
    #         temp_str2 += f"\t{factor2[temp_count]}"
    #         temp_count += 1
    #     temp_str2 += "\n"

    #     for key in sorted(count_str.keys()):
    #         temp_str2 += '\t'.join([temp_id[key], key.replace('-', '\t'), '\t'.join(count_str[key])])

    #     out5 = os.path.join(out2_dir, "Enhancer_FPKM.txt")
    #     with open(out5, 'w') as out_file:
    #         out_file.write(temp_str2)
    #     sort_file(out5)

    # outstr1 = ""
    # for item in key_str.keys():
    #     outstr1 += f"{item}\t{key_str[item]}\n"

    # out6 = os.path.join(out2_dir, "Enhancer_associated_peak.txt")
    # with open(out6, 'w') as out_file:
    #     out_file.write(outstr1)
    # sort_file(out6)

    # if prior == 1:
    #     out7 = os.path.join(out2_dir, "Enhancer_candidate.bed")
    #     with open(out7, 'w') as out_file:
    #         out_file.write(tempstr1)
    #     sort_file(out7)
    # import math

    # # Variables initialization
    # out_dir = "output_directory/"
    # enhcfile = out_dir + "eRNA/Enhancer_change.txt"
    # enhchange = {}
    # enhfdr = {}
    # flag23 = 0
    # temax = 0
    # emax = {}

    # # Reading the enhancer change file
    # with open(enhcfile, "r") as infile:
    #     for line in infile:
    #         line = line.strip()
    #         flag23 += 1
    #         if flag23 == 1:
    #             continue
    #         temp = line.split("\t")
    #         if (temp[5] not in ["NA", "Inf", "-Inf"]) and (abs(float(temp[5])) > temax):
    #             temax = abs(float(temp[5]))
    #         if temp[5] == "Inf":
    #             enhchange[temp[0]] = 1000
    #             emax[temp[0]] = 1
    #         elif temp[5] == "-Inf":
    #             enhchange[temp[0]] = -1000
    #             emax[temp[0]] = -1
    #         else:
    #             enhchange[temp[0]] = float(temp[5])
    #         enhfdr[temp[0]] = temp[-1]

    # # Update emax values based on temax
    # for key in emax:
    #     if emax[key] == 1:
    #         emax[key] = temax
    #     elif emax[key] == -1:
    #         emax[key] = -temax

    # # Variables for the second part of the script
    # flag24 = 0
    # outstr24 = ""

    # # Define other required variables
    # weight = 0.5  # Example value, update accordingly
    # fdr = 0.05  # Example value, update accordingly
    # direction = 1  # Example value, update accordingly
    # dis_to_p = {}  # Example dictionary, update accordingly
    # tchange = {}  # Example dictionary, update accordingly
    # tfdr = {}  # Example dictionary, update accordingly
    # out1 = "output1.txt"  # Example filename, update accordingly

    # # Reading the output file and processing
    # with open(out1, "r") as infile:
    #     for line in infile:
    #         line = line.strip()
    #         flag24 += 1
    #         if flag24 == 1:
    #             outstr24 += line + "\tFscore\tBscore\n"
    #             continue
    #         temp = line.split("\t")
    #         score = 0

    #         bscore = 0
    #         if temp[0] in dis_to_p:
    #             bscore = (weight * 2 / (1 + math.exp(0.0004054651 * dis_to_p[temp[0]])))

    #         fscore = 0
    #         if enhfdr[temp[0]] == "NA":
    #             score = bscore
    #         elif enhfdr[temp[0]] > fdr:
    #             score = bscore
    #         else:
    #             temp1 = temp[6].split(',;')
    #             f5 = {t: 1 for t in temp1}
    #             wd = temp[7].split(',')
    #             a4d = temp[8].split(',')
    #             cl = temp[12]

    #             f5gmax, wdgmax, a4dgmax, clgmax = 0, 0, 0, 0

    #             if direction == 1:
    #                 for key1 in f5:
    #                     if key1 == "NA":
    #                         continue
    #                     if key1 not in tchange:
    #                         continue
    #                     if tfdr[key1] == "NA":
    #                         continue
    #                     if tfdr[key1] > fdr:
    #                         continue
    #                     if enhchange[temp[0]] * tchange[key1] > 0 and abs(tchange[key1]) > abs(f5gmax):
    #                         f5gmax = tchange[key1]

    #                 for i in range(len(wd)):
    #                     if wd[i] == "NA":
    #                         continue
    #                     if tfdr[wd[i]] == "NA":
    #                         continue
    #                     if tfdr[wd[i]] > fdr:
    #                         continue
    #                     if enhchange[temp[0]] * tchange[wd[i]] > 0 and abs(tchange[wd[i]]) > abs(wdgmax):
    #                         wdgmax = tchange[wd[i]]

    #                 for i in range(len(a4d)):
    #                     if a4d[i] == "NA":
    #                         continue
    #                     if a4d[i] not in tchange:
    #                         continue
    #                     if tfdr[a4d[i]] == "NA":
    #                         continue
    #                     if tfdr[a4d[i]] > fdr:
    #                         continue
    #                     if enhchange[temp[0]] * tchange[a4d[i]] > 0 and abs(tchange[a4d[i]]) > abs(a4dgmax):
    #                         a4dgmax = tchange[a4d[i]]

    #                 if cl in tchange and tfdr[cl] != "NA" and tfdr[cl] <= fdr:
    #                     if enhchange[temp[0]] * tchange[cl] > 0:
    #                         clgmax = tchange[cl]

    #                 fscore = (enhchange[temp[0]] / temax) * ((f5gmax + wdgmax + a4dgmax + clgmax) / tgmax)

    #             elif direction == -1:
    #                 for key1 in f5:
    #                     if key1 == "NA":
    #                         continue
    #                     if key1 not in tchange:
    #                         continue
    #                     if tfdr[key1] == "NA":
    #                         continue
    #                     if tfdr[key1] > fdr:
    #                         continue
    #                     if enhchange[temp[0]] * tchange[key1] < 0 and abs(tchange[key1]) > abs(f5gmax):
    #                         f5gmax = tchange[key1]

    #                 for i in range(len(wd)):
    #                     if wd[i] == "NA":
    #                         continue
    #                     if tfdr[wd[i]] == "NA":
    #                         continue
    #                     if tfdr[wd[i]] > fdr:
    #                         continue
    #                     if enhchange[temp[0]] * tchange[wd[i]] < 0 and abs(tchange[wd[i]]) > abs(wdgmax):
    #                         wdgmax = tchange[wd[i]]

    #                 for i in range(len(a4d)):
    #                     if a4d[i] == "NA":
    #                         continue
    #                     if a4d[i] not in tchange:
    #                         continue
    #                     if tfdr[a4d[i]] == "NA":
    #                         continue
    #                     if tfdr[a4d[i]] > fdr:
    #                         continue
    #                     if enhchange[temp[0]] * tchange[a4d[i]] < 0 and abs(tchange[a4d[i]]) > abs(a4dgmax):
    #                         a4dgmax = tchange[a4d[i]]

    #                 if cl in tchange and tfdr[cl] != "NA" and tfdr[cl] <= fdr:
    #                     if enhchange[temp[0]] * tchange[cl] < 0:
    #                         clgmax = tchange[cl]

    #                 fscore = (-1) * (enhchange[temp[0]] / temax) * ((f5gmax + wdgmax + a4dgmax + clgmax) / tgmax)

    #             else:
    #                 for key1 in f5:
    #                     if key1 == "NA":
    #                         continue
    #                     if key1 not in tchange:
    #                         continue
    #                     if tfdr[key1] == "NA":
    #                         continue
    #                     if tfdr[key1] > fdr:
    #                         continue
    #                     if abs(tchange[key1]) > abs(f5gmax):
    #                         f5gmax = tchange[key1]

    #                 for i in range(len(wd)):
    #                     if wd[i] == "NA":
    #                         continue
    #                     if tfdr[wd[i]] == "NA":
    #                         continue
    #                     if tfdr[wd[i]] > fdr:
    #                         continue
    #                     if abs(tchange[wd[i]]) > abs(wdgmax):
    #                         wdgmax = tchange[wd[i]]

    #                 for i in range(len(a4d)):
    #                     if a4d[i] == "NA":
    #                         continue
    #                     if a4d[i] not in tchange:
    #                         continue
    #                     if tfdr[a4d[i]] == "NA":
    #                         continue
    #                     if tfdr[a4d[i]] > fdr:
    #                         continue
    #                     if abs(tchange[a4d[i]]) > abs(a4dgmax):
    #                         a4dgmax = tchange[a4d[i]]

    #                 if cl in tchange and tfdr[cl] != "NA" and tfdr[cl] <= fdr:
    #                     clgmax = tchange[cl]

    #                 fscore = (abs(enhchange[temp[0]]) / temax) * ((abs(f5gmax) + abs(wdgmax) + abs(a4dgmax) + abs(clgmax)) / tgmax)

    #         score = bscore
    #         outstr24 += f"{line}\t{fscore}\t{bscore/weight}\n"

    # # Writing the output
    # with open(out1, "w") as outfile:
    #     outfile.write(outstr24)


if __name__ == "__main__":
    args = getarg()

    main()
