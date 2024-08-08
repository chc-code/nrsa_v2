

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

logger = getlogger()

def getarg():
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('-in1', help="""required, read alignment files in bed (6 columns) or bam format for condition1, separated by space""", nargs='+')
    ps.add_argument('-in2', help="""read alignment files in bed (6 columns) or bam format for condition2, separated by space""", nargs='*')
    ps.add_argument('-o', help="""work directory, should be the same of pause_PROseq.pl\'s output/work directory""", nargs=1)
    ps.add_argument('-organism', '-m', '-org',  help="""define the genome: hg19, hg38, mm10, dm3, dm6, ce10, or danRer10. default: hg19""", choices=['hg19', 'hg38', 'mm10', 'dm3', 'dm6', 'ce10', 'danRer10'], default='hg19')


    ps.add_argument('-overlap_frac', '-p', help="""percentage of overlap for calling annotated gene, float, default=0.2""", type=float, default=0.2)
    ps.add_argument('-cutoff', '-c', help="""distance cutoff for divergent transcripts for enhancer detection (bp, default: 400)""", type=int, default=400)
    ps.add_argument('-distance', 'd', help="""distance within which two eRNAs are merged (bp, default: 500)""", type=int, default=500)
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

def check_args(args):
    code = 0
    overlap_frac = args.overlap_frac
    control_bed = args.in1
    case_bed = args.in2
    genome = args.genome
    direction = args.dir
    prior = args.pri
    peakfile = args.peak
    linkage = args.lk
    fdr = args.fdr
    weight = args.wt

    # reference files
    # my $gtf = $dir."\/hg19\/RefSeq-hg19-exon-formatted.gtf";		#reference genes in gtf format;
    # my $fantom=$dir."\/annotation\/human_permissive_enhancers_phase_1_and_2.bed";
    # my $association=$dir."\/annotation\/human_enhancer_tss_associations.bed";
    # my $tss=$dir."\/hg19\/RefSeq-hg19-tss.txt";
    # my $tss_tts=$dir."\/hg19\/RefSeq-hg19-tss-tts.txt";
    # my $fdgenome=$dir."\/annotation\/4DGenome-hg19.txt";


    # check intermediate folder
    pw_intermediate = os.path.join(pw_data, 'intermediate')
    if not os.path.exists(pw_intermediate):
        logger.error(f"Not find the results of pause_PROseq.pl, please check the work directory.")
        return 1, None


    # check HOMER_tag
    pw_homer = os.path.join(pw_intermediate, 'HOMER_tag')

    # 1) overlap_frac

    if overlap_frac <= 0 or overlap_frac > 1:
        logger.error("Overlap for calling annotated gene should be less than 0-1")
        code = 1
    # 4) bed files checking
    if control_bed:
        cond1 = control_bed
        for i in range(len(cond1)):
            if not os.path.exists(cond1[i]):
                logger.error(f"Input experiment file in bed/bam format (-in1) '{cond1[i]}' does not exist\n")

            temp1 = cond1[i].rindex(".")
            temp2 = cond1[i][temp1+1:]
            if temp2 not in ["bed", "bam"]:
                logger.error(f"Please input experiment file '{cond1[i]}' in bed or bam format!\n")
                code = 1
                
            cond1[i] = os.path.abspath(cond1[i])

    else:
        logger.error("Input experiment files in bed/bam format (-in1) has not been defined yet\n")
        code = 1

    if case_bed:
        cond2 = case_bed
        for i in range(len(cond2)):
            if not os.path.exists(cond2[i]):
                logger.error(f"Input experiment file in bed/bam format (-in2) '{cond2[i]}' does not exist\n")
                code = 1
            temp1 = cond2[i].rindex(".")
            temp2 = cond2[i][temp1+1:]
            if temp2 not in ["bed", "bam"]:
                logger.error(f"Please input experiment file '{cond2[i]}' in bed or bam format!\n")
                code = 1
            cond2[i] = os.path.abspath(cond2[i])

    # 5) normalize factors checking

    # 6) check genome
    valid_genomes = ["hg19", "mm10", "dm3", "ce10", "danRer10", "hg38", "dm6"]
    if genome not in valid_genomes:
        logger.error("Please set genome (-m) as hg19, hg38, mm10, dm3, dm6, ce10, or danRer10!\n")
        code = 1

    # get the path of the reference files
    ref_files = {}
    pw_ref = os.path.join(pw_data, 'ref')
    pw_fa = os.path.join(pw_data, 'fa')
    pw_annotation = os.path.join(pw_data, 'annotation')
    
    for igenome in valid_genomes:
        ref_files[igenome] = {}
        ref_files[igenome]["gtf"] = os.path.join(pw_ref, igenome, f"RefSeq-{igenome}-exon-formatted.gtf")
        ref_files[igenome]["tss"] = os.path.join(pw_ref, igenome, f"RefSeq-{igenome}-tss.txt")
        ref_files[igenome]["tss_tts"] = os.path.join(pw_ref, igenome, f"RefSeq-{igenome}-tss-tts.txt")
        
        ref_files[igenome]["genome_seq"] = os.path.join(pw_fa, igenome, f"{igenome}.fa")

        ref_files[igenome]["fdgenome"] = os.path.join(pw_annotation,  f"4DGenome-{igenome}.txt")
        ref_files[igenome]["fantom"] = os.path.join(pw_annotation, "human_permissive_enhancers_phase_1_and_2.bed")
        ref_files[igenome]["association"] = os.path.join(pw_annotation, "human_enhancer_tss_associations.bed")

        # validite file existence
        for key, value in ref_files[genome].items():
            if not os.path.exists(value):
                if pw_annotation not in value:
                    logger.warning(f"{igenome} - {key}: '{value}' does not exist\n")
                    ref_files[igenome][key] = None
                else:
                    logger.error(f"{igenome} - {key}: '{value}' does not exist\n")
                    code = 1

    # 7) peak file checking
    if prior == 1:
        if peakfile is None:
            logger.error("Peak file for sorting enhancers (-peak) has not been defined yet\n")
            code = 1
        else:
            if not os.path.exists(peakfile):
                logger.error(f"Peak file (-peak) '{peakfile}' does not exist\n")
                code = 1
            else:
                peakfile = os.path.abspath(peakfile)

    # 8) prioritize or not checking
    if prior not in [0, 1]:
        logger.error("(-pri) should be 0 (do not prioritize) or 1 (do prioritize)!\n")
        code = 1

    # 9) direction checking
    if direction not in [1, -1, 0]:
        logger.error("(-dir) should be 1 (consistent), -1 (inconsistent), or 0 (don't care)!\n")
        code = 1

    # 10) weight checking
    if prior == 1:
        if not (0 <= weight <= 1):
            logger.error("Weight (-wt) should be in [0, 1]!\n")
            code = 1

    # 11) cutoff of fdr checking
    if prior == 1:
        if not (0 <= fdr <= 1):
            logger.error("Cutoff of FDR (-cf) should be in [0, 1]!\n")
            code = 1

    # 12) linkage checking
    if prior == 1:
        if linkage not in ["pp", "gb", "pindex"]:
            logger.error("Type of linkage (-lk) should be pp, gb, or pindex!\n")
            code = 1

    return code, ref_files


def gtf_compare(gtf_in, gtf_ref, overlap_frac):
    pass

def region_compare(other_gene, out_txt):
    pass

def central(other_region, lcut, target_type):
    # lcut: cutoff of 5' distance;
    # target_type can be eRNA or long_eRNA
    
    fno_enhancer = f'{out_dir}/{target_type}/Enhancer.txt'
    enhancer_found = 0
    # build the enhancer region content here
    # sort the enhancer region
    
    
    if enhancer_found == 0:
        logger.error("No enhancer region was detected! Program exits!")
        return 1, None
    
    pass
findEnhancer = central

def main():
    overlap_frac = args.overlap_frac # percentage of overlap for calling annotated gene, default 0.2
    genome = args.genome
    lcut = args.cutoff #cutoff of 5' distance;
    d_e = args.distance #distance between two eRNAs for merging(bp, default: 500)
    l_erna = args.le #length cutoff for long-eRNA identification
    out_dir = args.w # work directory, should be the same of pause_PROseq.pl's output/work directory
    wd = args.wd #distance within which associate active genes of enhancer (bp, default: 50000)

    toggle_filter = args.filter  # whether to filter for enhancers: 0 (not filter) or 1 (filter)
    toggle_prior = args.pri
    
    filter_tss = args.dtss
    filter_tts = args.dtts
    wd = args.wd  #target active gene within distance
    peakfile = args.peak
    direction = args.dir
    weight = args.wt
    fdr = args.fdr
    linkage = args.lk; #pp, gb, or pindex

    dependency_exist = check_dependency()
    if not dependency_exist:
        return 1

    # check arguments
    check_return, ref_files = check_args(args)
    if check_return:
        sys.exit(1)

    # check the input files, if bam, convert to bed
    control_bed, case_bed = process_input(args)
    if control_bed is None or case_bed is None:
        return 1
    
    # create the tag files
    pw_intermediate = os.path.join(pw_data, 'intermediate')
    pw_tag = os.path.join(pw_intermediate, 'HOMER_tag')
    
    
    bed_files = control_bed + case_bed
    # use forceBED to ignore the 5th column
    cmd_make_tag = f'makeTagDirectory {pw_tag} -format bed -forceBED {" ".join(bed_files)}'
    
    logger.info(f'Creating tag directory using HOMER')
    run_status = os.system(cmd_make_tag)
    if run_status != 0:
        logger.error(f"Failed to create tag directory for HOMER, please check")
        return 1
    logger.info(f'done.')
    
    # findPeaks, this step is quite fast
    # will generate the psuedo-transcript according to the alignment
    out_txt = os.path.join(pw_intermediate, 'findpeaks.transcript.txt')
    out_gtf = os.path.join(pw_intermediate, 'findpeaks.transcript.gtf')
    logger.info(f'Finding peaks from input files')
    cmd_findpeaks = f'findPeaks {pw_tag} -style groseq -o {out_txt} -gtf {out_gtf}'
    run_status = os.system(cmd_findpeaks)
    if run_status != 0:
        logger.error(f"Failed to run findPeaks for HOMER, please check")
        return 1
    logger.info(f'done.')

    # compare the gtf
    other_gene = gtf_compare(out_gtf, ref_files[genome]["gtf"], overlap_frac)
    other_region = region_compare(other_gene, out_txt)
    
    # find the enhancer region
    logger.info(f'Finding enhancer region')
    fno_enhancer, outstr3 = central(other_region, lcut, 'eRNA')
    
    if fno_enhancer == 1:
        return 1
    
    out1 = os.path.join(out2_dir, "Enhancer.txt")
    with open(out1, "w") as out1_file:
        out1_file.write(outstr1)

    sort_file(out1)

    print("Finding closest genes...")
    gene_id = {}
    with open(gtf) as in9_file:
        for line in in9_file:
            gene = line.split()
            if gene[0].startswith("chrUn") or gene[0].startswith("Un") or gene[2] != "exon":
                continue
            gene[11] = gene[11].strip("\";")
            gene[13] = gene[13].strip("\";")
            gene_id[gene[11]] = gene[13]

    active = os.path.join(inter_dir, "active_gene.txt")
    active_gene = {}
    with open(active) as in3_file:
        for line in in3_file:
            temp = line.split()
            active_gene[temp[0]] = 1

    active_tss_str = ""
    with open(tss) as in4_file:
        for line in in4_file:
            temp = line.split("\t")
            if temp[3] in active_gene:
                active_tss_str += line

    active_tss = os.path.join(inter_dir, "active_tss.txt")
    with open(active_tss, "w") as out5_file:
        out5_file.write(active_tss_str)

    sort_file(active_tss)

    closest_out = os.path.join(inter_dir, "closest.txt")
    subprocess.run(f"bedtools closest -a {out1} -b {active_tss} -d > {closest_out}", shell=True, check=True)


    # Initialize dictionaries
    clos_gene = {}
    distance = {}

    center_str = {}
    center_enhancer = {}

    # Read and process the closest_out file
    with open(closest_out, 'r') as infile:
        for line in infile:
            temp = line.strip().split('\t')
            temp_str = "\t".join(temp[0:6])
            if int(temp[11]) == -1:
                if temp_str not in clos_gene:
                    clos_gene[temp_str] = []
                clos_gene[temp_str].append("NA")
                distance[temp_str] = temp[11]
                continue
            if temp[9] in gene_id and gene_id[temp[9]]:
                if temp_str not in clos_gene:
                    clos_gene[temp_str] = []
                clos_gene[temp_str].append(gene_id[temp[9]])
                distance[temp_str] = temp[11]
            else:
                if temp_str not in clos_gene:
                    clos_gene[temp_str] = []
                clos_gene[temp_str].append(temp[9])
                distance[temp_str] = temp[11]
            
            if ',' in temp[3]:
                temp1 = temp[3].split(',')
                temp2 = temp[4].split(',')
                for i in range(len(temp1)):
                    key = f"{temp[0]}-{temp1[i]}"
                    if key not in center_str:
                        center_str[key] = f"{temp[0]}\t{temp1[i]}\t{temp2[i]}"
                        center_enhancer[key] = temp_str
            else:
                key = f"{temp[0]}-{temp[3]}"
                if key not in center_str:
                    center_str[key] = f"{temp[0]}\t{temp[3]}\t{temp[4]}"
                    center_enhancer[key] = temp_str

    print("Finding associate genes within 50Kb...")
    tss_chr = []
    tss_pos = []
    tss_gene = []

    # Read active TSS file
    with open(active_tss, 'r') as infile:
        for line in infile:
            temp = line.strip().split('\t')
            if gene_id.get(temp[3]):
                tss_chr.append(temp[0])
                tss_pos.append(int(temp[1]))
                tss_gene.append(gene_id[temp[3]])

    # Initialize dictionaries for 50Kb proximity checks
    tss50k = {}
    etss50k = {}

    wd = 50000  # Window size of 50Kb

    # Process each enhancer
    for enhancer in sorted(distance.keys()):
        temp = enhancer.split('\t')
        chr = temp[0]
        start = int(temp[1])
        end = int(temp[2])
        
        flag1 = 0
        for i in range(len(tss_chr)):
            if flag1 == 0 and chr != tss_chr[i]:
                continue
            if flag1 == 1 and chr != tss_chr[i]:
                break
            if chr == tss_chr[i]:
                flag1 = 1
                if tss_pos[i] > end and (tss_pos[i] - end) < wd:
                    if enhancer not in etss50k or tss_gene[i] not in etss50k[enhancer]:
                        if enhancer in tss50k:
                            tss50k[enhancer] += f",{tss_gene[i]}"
                        else:
                            tss50k[enhancer] = tss_gene[i]
                        if enhancer not in etss50k:
                            etss50k[enhancer] = {}
                        etss50k[enhancer][tss_gene[i]] = 1
                if tss_pos[i] < start and (start - tss_pos[i]) < wd:
                    if enhancer not in etss50k or tss_gene[i] not in etss50k[enhancer]:
                        if enhancer in tss50k:
                            tss50k[enhancer] += f",{tss_gene[i]}"
                        else:
                            tss50k[enhancer] = tss_gene[i]
                        if enhancer not in etss50k:
                            etss50k[enhancer] = {}
                        etss50k[enhancer][tss_gene[i]] = 1


    chr_4d = []
    start_4d = []
    end_4d = []
    gene_4d = []
    cell_4d = []
    method_4d = []
    pmid_4d = []
    fd = {}
    e4d = {}
    cell_4d_map = {}
    ecell_4d = {}
    method_4d_map = {}
    emethod_4d = {}
    pmid_4d_map = {}
    epmid_4d = {}

    if genome not in ["ce10", "danRer10"]:
        print("Finding associate genes in 4DGenome...")

        with open(fdgenome, 'r') as file:
            for line in file:
                temp = line.strip().split('\t')
                chr_4d.append(temp[0])
                start_4d.append(temp[1])
                end_4d.append(temp[2])
                gene_4d.append(temp[6])
                cell_4d.append(temp[7])
                method_4d.append(temp[8])
                pmid_4d.append(temp[9])

        for enhancer in sorted(distance.keys()):
            temp = enhancer.split('\t')
            chr = temp[0]
            start = int(temp[1])
            end = int(temp[2])

            flag1 = 0
            for i in range(len(chr_4d)):
                if flag1 == 0 and chr != chr_4d[i]:
                    continue
                if flag1 == 1 and chr != chr_4d[i]:
                    break
                if chr == chr_4d[i]:
                    flag1 = 1
                    if (start_4d[i] >= start and start_4d[i] <= end) or (end_4d[i] >= start and end_4d[i] <= end) or (start >= start_4d[i] and start <= end_4d[i]) or (end >= start_4d[i] and end <= end_4d[i]):
                        temp = gene_4d[i].split(';')
                        for gene in temp:
                            if enhancer not in e4d or gene not in e4d[enhancer]:
                                fd[enhancer] = fd.get(enhancer, '') + ("," if enhancer in fd else '') + gene
                                e4d.setdefault(enhancer, {})[gene] = 1

                        if enhancer not in ecell_4d or cell_4d[i] not in ecell_4d[enhancer]:
                            cell_4d_map[enhancer] = cell_4d_map.get(enhancer, '') + ("," if enhancer in cell_4d_map else '') + cell_4d[i]
                            ecell_4d.setdefault(enhancer, {})[cell_4d[i]] = 1

                        if enhancer not in emethod_4d or method_4d[i] not in emethod_4d[enhancer]:
                            method_4d_map[enhancer] = method_4d_map.get(enhancer, '') + ("," if enhancer in method_4d_map else '') + method_4d[i]
                            emethod_4d.setdefault(enhancer, {})[method_4d[i]] = 1

                        if enhancer not in epmid_4d or pmid_4d[i] not in epmid_4d[enhancer]:
                            pmid_4d_map[enhancer] = pmid_4d_map.get(enhancer, '') + ("," if enhancer in pmid_4d_map else '') + pmid_4d[i]
                            epmid_4d.setdefault(enhancer, {})[pmid_4d[i]] = 1

    enhancer_id = 1
    enhancer_id_map = {}

    outstr5 = "Enhancer_ID\tchr\tstart\tend\tcenter\tFANTOM5\tassociated_gene-FANTOM5\tassociated_gene-50kb\tassociated_gene-4DGenome\tCell/Tissue\tDetection_method\tPMID\tclosest_gene\tdistance\n"
    for enhancer in sorted(distance.keys()):
        temp_50k = tss50k.get(enhancer, "NA")

        if genome not in ["ce10", "danRer10"]:
            if enhancer in fd:
                temp1_4d = fd[enhancer]
                temp2_4d = cell_4d_map[enhancer]
                temp3_4d = method_4d_map[enhancer]
                temp4_4d = pmid_4d_map[enhancer]
            else:
                temp1_4d = "NA"
                temp2_4d = "NA"
                temp3_4d = "NA"
                temp4_4d = "NA"
            outstr5 += f"{enhancer_id}\t{enhancer}\t{temp_50k}\t{temp1_4d}\t{temp2_4d}\t{temp3_4d}\t{temp4_4d}\t{','.join(clos_gene[enhancer])}\t{distance[enhancer]}\n"
        else:
            outstr5 += f"{enhancer_id}\t{enhancer}\t{temp_50k}\tNA\tNA\tNA\tNA\t{','.join(clos_gene[enhancer])}\t{distance[enhancer]}\n"

        enhancer_id_map[enhancer] = enhancer_id
        enhancer_id += 1

    with open(out1, 'w') as out_file:
        out_file.write(outstr5)

    outstr2 = ""
    tempstr1 = ""
    for center_key in sorted(center_str.keys()):
        outstr2 += f"{center_str[center_key]}\t{enhancer_id_map[center_enhancer[center_key]]}\n"
        if prior == 1:
            temp = center_key.split('-')
            tempstr1 += f"{temp[0]}\t{temp[1]}\t{temp[1]}\t{enhancer_id_map[center_enhancer[center_key]]}\tNA\t+\n"

    out2 = os.path.join(out2_dir, "Enhancer_center.txt")
    with open(out2, 'w') as out_file:
        out_file.write(outstr2)
    sort_file(out2)

    out3 = os.path.join(out2_dir, "long_eRNA.txt")
    if outstr3:
        with open(out3, 'w') as out_file:
            out_file.write(outstr3)
        sort_file(out3)

    temp_str = ""
    temp_id = {}
    line = 0
    with open(out1, 'r') as in_file:
        for line_num, line in enumerate(in_file):
            if line_num == 0:
                continue
            temp = line.strip().split('\t')
            temp_str += f"{temp[1]}\t{temp[2]}\t{temp[3]}\n"
            temp_id[f"{temp[1]}-{temp[2]}-{temp[3]}"] = temp[0]

    out4 = os.path.join(inter_dir, "Enhancer_temp.bed")
    with open(out4, 'w') as out_file:
        out_file.write(temp_str)

    print("Detecting Enhancer change...")
    if case_bed:
        factor1 = []
        factor2 = []
        nffile = os.path.join(inter_dir, "nf.txt")
        factor = {}
        with open(nffile, 'r') as in_file:
            for line in in_file:
                temp = line.strip().split('\t')
                factor[temp[0]] = temp[1]

        count_str = {}
        rep_str = ""

        for rep in cond1:
            temp1 = rep.rindex("/")
            temp2 = rep.rindex(".")
            temp3 = rep[temp1 + 1:temp2]

            if temp3 in factor:
                factor1.append(factor[temp3])
            else:
                print("input file name doesn't match the file name for pausing analysis!")
                exit()

            rep_str += f"\t{temp3}"
            cov_temp = os.path.join(inter_dir, "cov_temp.bed")
            cov_cmd1 = f"bedtools coverage -a {out4} -b {rep} > {cov_temp}"
            subprocess.run(cov_cmd1, shell=True, check=True)
            with open(cov_temp, 'r') as in_file:
                for line in in_file:
                    temp = line.strip().split()
                    key = f"{temp[0]}-{temp[1]}-{temp[2]}"
                    count_str.setdefault(key, []).append(temp[3])

        for rep in cond2:
            temp1 = rep.rindex("/")
            temp2 = rep.rindex(".")
            temp3 = rep[temp1 + 1:temp2]

            if temp3 in factor:
                factor2.append(factor[temp3])
            else:
                print("input file name doesn't match the file name for pausing analysis!")
                exit()

            rep_str += f"\t{temp3}"
            cov_temp = os.path.join(inter_dir, "cov_temp.bed")
            cov_cmd1 = f"bedtools coverage -a {out4} -b {rep} > {cov_temp}"
            subprocess.run(cov_cmd1, shell=True, check=True)
            with open(cov_temp, 'r') as in_file:
                for line in in_file:
                    temp = line.strip().split()
                    key = f"{temp[0]}-{temp[1]}-{temp[2]}"
                    count_str.setdefault(key, []).append(temp[3])

        temp_count = 0
        temp_str2 = "Enhancer_ID\tchr\tstart\tend"
        for rep in cond1:
            temp_str2 += f"\t{factor1[temp_count]}"
            temp_count += 1
        for rep in cond2:
            temp_str2 += f"\t{factor2[temp_count]}"
            temp_count += 1
        temp_str2 += "\n"

        for key in sorted(count_str.keys()):
            temp_str2 += '\t'.join([temp_id[key], key.replace('-', '\t'), '\t'.join(count_str[key])])

        out5 = os.path.join(out2_dir, "Enhancer_FPKM.txt")
        with open(out5, 'w') as out_file:
            out_file.write(temp_str2)
        sort_file(out5)

    outstr1 = ""
    for item in key_str.keys():
        outstr1 += f"{item}\t{key_str[item]}\n"

    out6 = os.path.join(out2_dir, "Enhancer_associated_peak.txt")
    with open(out6, 'w') as out_file:
        out_file.write(outstr1)
    sort_file(out6)

    if prior == 1:
        out7 = os.path.join(out2_dir, "Enhancer_candidate.bed")
        with open(out7, 'w') as out_file:
            out_file.write(tempstr1)
        sort_file(out7)
    import math

    # Variables initialization
    out_dir = "output_directory/"
    enhcfile = out_dir + "eRNA/Enhancer_change.txt"
    enhchange = {}
    enhfdr = {}
    flag23 = 0
    temax = 0
    emax = {}

    # Reading the enhancer change file
    with open(enhcfile, "r") as infile:
        for line in infile:
            line = line.strip()
            flag23 += 1
            if flag23 == 1:
                continue
            temp = line.split("\t")
            if (temp[5] not in ["NA", "Inf", "-Inf"]) and (abs(float(temp[5])) > temax):
                temax = abs(float(temp[5]))
            if temp[5] == "Inf":
                enhchange[temp[0]] = 1000
                emax[temp[0]] = 1
            elif temp[5] == "-Inf":
                enhchange[temp[0]] = -1000
                emax[temp[0]] = -1
            else:
                enhchange[temp[0]] = float(temp[5])
            enhfdr[temp[0]] = temp[-1]

    # Update emax values based on temax
    for key in emax:
        if emax[key] == 1:
            emax[key] = temax
        elif emax[key] == -1:
            emax[key] = -temax

    # Variables for the second part of the script
    flag24 = 0
    outstr24 = ""

    # Define other required variables
    weight = 0.5  # Example value, update accordingly
    fdr = 0.05  # Example value, update accordingly
    direction = 1  # Example value, update accordingly
    dis_to_p = {}  # Example dictionary, update accordingly
    tchange = {}  # Example dictionary, update accordingly
    tfdr = {}  # Example dictionary, update accordingly
    out1 = "output1.txt"  # Example filename, update accordingly

    # Reading the output file and processing
    with open(out1, "r") as infile:
        for line in infile:
            line = line.strip()
            flag24 += 1
            if flag24 == 1:
                outstr24 += line + "\tFscore\tBscore\n"
                continue
            temp = line.split("\t")
            score = 0

            bscore = 0
            if temp[0] in dis_to_p:
                bscore = (weight * 2 / (1 + math.exp(0.0004054651 * dis_to_p[temp[0]])))

            fscore = 0
            if enhfdr[temp[0]] == "NA":
                score = bscore
            elif enhfdr[temp[0]] > fdr:
                score = bscore
            else:
                temp1 = temp[6].split(',;')
                f5 = {t: 1 for t in temp1}
                wd = temp[7].split(',')
                a4d = temp[8].split(',')
                cl = temp[12]

                f5gmax, wdgmax, a4dgmax, clgmax = 0, 0, 0, 0

                if direction == 1:
                    for key1 in f5:
                        if key1 == "NA":
                            continue
                        if key1 not in tchange:
                            continue
                        if tfdr[key1] == "NA":
                            continue
                        if tfdr[key1] > fdr:
                            continue
                        if enhchange[temp[0]] * tchange[key1] > 0 and abs(tchange[key1]) > abs(f5gmax):
                            f5gmax = tchange[key1]

                    for i in range(len(wd)):
                        if wd[i] == "NA":
                            continue
                        if tfdr[wd[i]] == "NA":
                            continue
                        if tfdr[wd[i]] > fdr:
                            continue
                        if enhchange[temp[0]] * tchange[wd[i]] > 0 and abs(tchange[wd[i]]) > abs(wdgmax):
                            wdgmax = tchange[wd[i]]

                    for i in range(len(a4d)):
                        if a4d[i] == "NA":
                            continue
                        if a4d[i] not in tchange:
                            continue
                        if tfdr[a4d[i]] == "NA":
                            continue
                        if tfdr[a4d[i]] > fdr:
                            continue
                        if enhchange[temp[0]] * tchange[a4d[i]] > 0 and abs(tchange[a4d[i]]) > abs(a4dgmax):
                            a4dgmax = tchange[a4d[i]]

                    if cl in tchange and tfdr[cl] != "NA" and tfdr[cl] <= fdr:
                        if enhchange[temp[0]] * tchange[cl] > 0:
                            clgmax = tchange[cl]

                    fscore = (enhchange[temp[0]] / temax) * ((f5gmax + wdgmax + a4dgmax + clgmax) / tgmax)

                elif direction == -1:
                    for key1 in f5:
                        if key1 == "NA":
                            continue
                        if key1 not in tchange:
                            continue
                        if tfdr[key1] == "NA":
                            continue
                        if tfdr[key1] > fdr:
                            continue
                        if enhchange[temp[0]] * tchange[key1] < 0 and abs(tchange[key1]) > abs(f5gmax):
                            f5gmax = tchange[key1]

                    for i in range(len(wd)):
                        if wd[i] == "NA":
                            continue
                        if tfdr[wd[i]] == "NA":
                            continue
                        if tfdr[wd[i]] > fdr:
                            continue
                        if enhchange[temp[0]] * tchange[wd[i]] < 0 and abs(tchange[wd[i]]) > abs(wdgmax):
                            wdgmax = tchange[wd[i]]

                    for i in range(len(a4d)):
                        if a4d[i] == "NA":
                            continue
                        if a4d[i] not in tchange:
                            continue
                        if tfdr[a4d[i]] == "NA":
                            continue
                        if tfdr[a4d[i]] > fdr:
                            continue
                        if enhchange[temp[0]] * tchange[a4d[i]] < 0 and abs(tchange[a4d[i]]) > abs(a4dgmax):
                            a4dgmax = tchange[a4d[i]]

                    if cl in tchange and tfdr[cl] != "NA" and tfdr[cl] <= fdr:
                        if enhchange[temp[0]] * tchange[cl] < 0:
                            clgmax = tchange[cl]

                    fscore = (-1) * (enhchange[temp[0]] / temax) * ((f5gmax + wdgmax + a4dgmax + clgmax) / tgmax)

                else:
                    for key1 in f5:
                        if key1 == "NA":
                            continue
                        if key1 not in tchange:
                            continue
                        if tfdr[key1] == "NA":
                            continue
                        if tfdr[key1] > fdr:
                            continue
                        if abs(tchange[key1]) > abs(f5gmax):
                            f5gmax = tchange[key1]

                    for i in range(len(wd)):
                        if wd[i] == "NA":
                            continue
                        if tfdr[wd[i]] == "NA":
                            continue
                        if tfdr[wd[i]] > fdr:
                            continue
                        if abs(tchange[wd[i]]) > abs(wdgmax):
                            wdgmax = tchange[wd[i]]

                    for i in range(len(a4d)):
                        if a4d[i] == "NA":
                            continue
                        if a4d[i] not in tchange:
                            continue
                        if tfdr[a4d[i]] == "NA":
                            continue
                        if tfdr[a4d[i]] > fdr:
                            continue
                        if abs(tchange[a4d[i]]) > abs(a4dgmax):
                            a4dgmax = tchange[a4d[i]]

                    if cl in tchange and tfdr[cl] != "NA" and tfdr[cl] <= fdr:
                        clgmax = tchange[cl]

                    fscore = (abs(enhchange[temp[0]]) / temax) * ((abs(f5gmax) + abs(wdgmax) + abs(a4dgmax) + abs(clgmax)) / tgmax)

            score = bscore
            outstr24 += f"{line}\t{fscore}\t{bscore/weight}\n"

    # Writing the output
    with open(out1, "w") as outfile:
        outfile.write(outstr24)


if __name__ == "__main__":
    args = getarg()

    main()
