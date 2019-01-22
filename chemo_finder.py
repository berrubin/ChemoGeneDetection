import sys
import multiprocessing
from multiprocessing import Pool
from Bio import Seq
import subprocess
import os
from Bio import SeqIO
import numpy
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-p", "--num_threads", dest = "num_threads", type = int, default = 1, help = "Number of cores to use.")
parser.add_option("-c", "--genomes_file", dest = "genomes_file", type = str, help = "File with paths to genome fastas")
parser.add_option("-r", "--reference_file", dest = "reference_file", type = str, help = "Seed protein sequences")
parser.add_option("-f", "--family", dest = "family", default = "or", type = str, help = "Receptor family (e.g., or, gr)")
parser.add_option("-i", "--iteration", dest = "iteration", type = int, default = 0, help = "Iteration")
parser.add_option("-l", "--min_len", dest = "min_len", type = int, default = 350, help = "Minium length of gene")
parser.add_option("-b", "--base_dir", dest = "base_dir", type = str, help = "Directory in which to build output")
parser.add_option("-t", "--tblastn_dir", dest = "tblastn_dir", type = str, help = "tblastn output directory")
parser.add_option("-k", "--flank_size", dest = "flank_size", type = int, default = 1000, help = "length of flanking sequence beyond blast hit to include in protein alignment")
parser.add_option("-n", "--intron_size", dest = "intron_size", type = int, default = 12000, help = "Maximum intron size")
parser.add_option("-e", "--evalue", dest = "evalue", type = float, default = 100, help = "Maximum evalue of hits to include")
parser.add_option("-m", "--length_margin", dest = "length_margin", type = int, default = 40, help = "Acceptable difference in lengths between possible gene annotations. Otherwise, the longer sequence is always used even if it has a lower score.")
parser.add_option("-v", "--hit_overlap", dest = "hit_overlap", type = int, default = 30, help = "Maximum overlap in query sequence allowed for merging HSPs.")
parser.add_option("-s", "--prop_short", dest = "prop_short", type = float, default = None, help = "Proportion of the length of the seed protein that the inferred protein is allowed to be.")
parser.add_option("-x", "--exonerate_path", dest = "exonerate_path", type = str, default = "", help = "Absolute path to exonerate executable")

(options,args) = parser.parse_args()

CUR_ITER = options.iteration
MIN_LEN = options.min_len
PROP_SHORT = options.prop_short
if PROP_SHORT == None:
    USE_MIN_LEN = True
else:
    USE_MIN_LEN = False

FAMILY = options.family
PREDICTOR = "exonerate"
FLANK_SIZE = options.flank_size
INTRON_SIZE = options.intron_size
HIT_OVERLAP = options.hit_overlap
EVALUE = options.evalue
if EVALUE >= 1:
    EVALUE = int(EVALUE)
LENGTH_MARGIN = options.length_margin
EXONERATE_PATH = options.exonerate_path
GENES_DIR = "%s/%s_genes_iter_%s_%s_%s_%s_%s_%s" % (options.base_dir, FAMILY, CUR_ITER, EVALUE, FLANK_SIZE, INTRON_SIZE, HIT_OVERLAP, LENGTH_MARGIN)
REGIONAL_DIR = "%s/%s_regional_iter_%s_%s_%s_%s_%s_%s" % (options.base_dir, FAMILY, CUR_ITER, EVALUE, FLANK_SIZE, INTRON_SIZE, HIT_OVERLAP, LENGTH_MARGIN)
if EXONERATE_PATH == "":
    EXONERATE_PATH = "exonerate"
else:
    EXONERATE_PATH = "%s/exonerate" % EXONERATE_PATH

def main():
    if not os.path.exists(options.base_dir):
        os.mkdir(options.base_dir)
    work_list = []
    pool = multiprocessing.Pool(processes = options.num_threads)
    genome_dic = read_config(options.genomes_file)
    for species, genome in genome_dic.items():
        work_list.append([species, genome, FAMILY, CUR_ITER])
    pool.map_async(parse_blast, work_list).get(9999999)
    pool.close()
    pool.join()

def read_config(config_file):
    reader = open(config_file, 'rU')
    file_dic = {}
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.strip().split()
        file_dic[cur_line[0]] = cur_line[1]
    return file_dic

def parse_blast(param_list):
    species = param_list[0]
    genome_file = param_list[1]
    family = param_list[2]
    iteration = param_list[3]
    genome_dic = {}
    reader = SeqIO.parse(genome_file, format = 'fasta')
    
    for rec in reader:
        genome_dic[rec.id] = str(rec.seq)
    prot_dic = {}
    reader = SeqIO.parse(options.reference_file, format = 'fasta')
    for rec in reader:
        prot_dic[rec.id] = str(rec.seq)
    reader = open("%s/%s_%s_iter_%s_%s.txt" % (options.tblastn_dir, species, family, iteration, EVALUE), 'rU')
    file_dic = {}
    scaf_regions_dic = {}
    for line in reader:
        if "altto" in line:
            continue
        cur_line = line.split()
        if "PSE" in cur_line[0]:
            continue
        if float(cur_line[10]) > EVALUE:
            continue
        if cur_line[0] not in file_dic.keys():
            file_dic[cur_line[0]] = []
        file_dic[cur_line[0]].append(line.strip())
    for gene, lines in file_dic.items():
        hsp_scaf_dic = {}
        hsp_list = []
        for line in lines:
            cur_line = line.split()
            cur_scaf = cur_line[1]
            if cur_scaf not in hsp_scaf_dic.keys():
                hsp_scaf_dic[cur_scaf] = []
            hsp_scaf_dic[cur_scaf].append(HSP(int(cur_line[6]), int(cur_line[7]), int(cur_line[8]), int(cur_line[9]), cur_line[1], float(cur_line[13]), [float(cur_line[11])]))
        for scaf, hsp_list in hsp_scaf_dic.items():
            hsp_list = sort_hsps(hsp_list)
            temp_hsp_list = hsp_list
            cur_hsp_len = len(hsp_list)
            x = 0
            while x < len(hsp_list):
                temp_hsp_list = merge_HSPs(temp_hsp_list[x], temp_hsp_list)
                if len(temp_hsp_list) < cur_hsp_len:
                    x = 0
                else:
                    x += 1
                cur_hsp_len = len(temp_hsp_list)
            for hsp in temp_hsp_list:
                if iteration > 0:
                    if USE_MIN_LEN:
                        if hsp.qend - hsp.qstart < MIN_LEN:
                            continue
                    else:
                        if hsp.qend - hsp.qstart < PROP_SHORT*(len(prot_dic[gene])):
                            continue
                elif iteration == 0:
                    if USE_MIN_LEN:
                        if hsp.qend - hsp.qstart < MIN_LEN / 2.0:
                            continue
                    else:
                        if hsp.qend - hsp.qstart < PROP_SHORT*(len(prot_dic[gene]) / 2.0):
                            continue
                if hsp.scaf not in scaf_regions_dic.keys():
                    scaf_regions_dic[hsp.scaf] = []
                if USE_MIN_LEN:
                    if hsp.send - hsp.sstart < PROP_SHORT*(len(prot_dic[gene])*3.0):
                        continue
                else:
                    if hsp.send - hsp.sstart < 3*MIN_LEN:
                        continue
                scaf_regions_dic[hsp.scaf].append(GeneRegion(hsp.sstart, hsp.send, hsp.scaf, hsp.frame, numpy.sum(hsp.bits), gene))
    ors_list = []
    if not os.path.isdir(GENES_DIR):
        os.mkdir(GENES_DIR)
    if USE_MIN_LEN:
        length_param = MIN_LEN
    else:
        length_param = PROP_SHORT
    success_file = open("%s/%s_iter_%s_%s_%s_pep.faa" % (GENES_DIR, species, iteration, family, length_param), 'w')
    success_cds = open("%s/%s_iter_%s_%s_%s.fna" % (GENES_DIR, species, iteration, family, length_param), 'w')
    success_trans = open("%s/%s_iter_%s_%s_%s_trans.fna" % (GENES_DIR, species, iteration, family, length_param), 'w')
    success_gff = open("%s/%s_iter_%s_%s_%s.gff" % (GENES_DIR, species, iteration, family, length_param), 'w')
    if not os.path.isdir(REGIONAL_DIR):
        os.mkdir(REGIONAL_DIR)
    if not os.path.isdir("%s/%s" % (REGIONAL_DIR, species)):
        os.mkdir("%s/%s" % (REGIONAL_DIR, species))

    gene_index = 1
    for scaf, region_list in scaf_regions_dic.items():
        x = 0
        y = 0
        new_region_list = []
        for reg in region_list:
            reg_len, align_score, new_start, new_end = protein_length(reg, species, family, iteration, prot_dic, genome_dic)
            reg.prot_len = reg_len
            reg.align_score = align_score
            if USE_MIN_LEN:
                if reg_len > MIN_LEN:
                    reg.start = new_start
                    reg.end = new_end
                    new_region_list.append(reg)
            else:
                if reg_len > PROP_SHORT * (len(prot_dic[reg.query])):
                    reg.start = new_start
                    reg.end = new_end
                    new_region_list.append(reg)
        while x < len(new_region_list):
            y = 0
            while y < len(new_region_list):
                if x == y:
                    y += 1
                    continue
                region1 = new_region_list[x]
                region2 = new_region_list[y]
                if region1.overlap(region2):

                    if region1.prot_len > region2.prot_len - LENGTH_MARGIN and region1.prot_len < region2.prot_len + LENGTH_MARGIN: #if lengths are equal
                        region1_sublen = region1.end - float(region1.start)
                        region2_sublen = region2.end - float(region2.start)
                        if region1.align_score / region1_sublen > region2.align_score / region2_sublen:
                            del new_region_list[y]
                            x = 0
                            y = 0
                            continue
                        else:
                            del new_region_list[x]
                            x = 0
                            y = 0
                            continue

                    elif region1.prot_len > region2.prot_len:
                        del new_region_list[y]
                        x = 0
                        y = 0
                        continue
                    else:
                        del new_region_list[x]
                        x = 0
                        y = 0
                        continue
                y += 1
            x += 1
        print scaf
        for region in new_region_list:
            ref_file = open("%s/%s/%s_%s_%s_ref.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end), 'w')
            ref_file.write(">%s\n%s\n" % (region.query, prot_dic[region.query]))
            ref_file.close()

            outfile = open("%s/%s/%s_%s_%s.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end), 'w')

            cur_start = region.start - FLANK_SIZE
            if cur_start < 0:
                cur_start = 0
            cur_end = region.end + FLANK_SIZE
            if cur_end > len(genome_dic[region.scaf]):
                cur_end = len(genome_dic[region.scaf])
            outfile.write(">%s_%s_%s_%s\n%s\n" % (region.scaf, region.start, region.end, region.query[0:20], genome_dic[region.scaf][cur_start:cur_end]))
            outfile.close()
            if PREDICTOR == "exonerate":
                cmd = [EXONERATE_PATH, "--forcegtag", "--exhaustive", "y", "--model", "protein2genome", "--showtargetgff", "true", "--verbose", "0", "--showalignment", "no", "--showvulgar", "no", "--bestn", "1", "%s/%s/%s_%s_%s_ref.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end), "%s/%s/%s_%s_%s.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end)]
                FNULL = open(os.devnull, 'w')
                with open("%s/%s/%s_%s_%s.gff" % (REGIONAL_DIR, species, region.scaf, region.start, region.end), 'w') as outfile:
                    subprocess.call(cmd, stdout = outfile, stderr=FNULL)
                outfile.close()
                gene_id, prot_seq, cds_seq, align_score, new_start, new_end, trans_seq = reconstruct_gene("%s_%s_%s" % (region.scaf, region.start, region.end), species, iteration, family, region.scaf, region.start, region.end, gene_index)
                if "*" in prot_seq:
                    prot_split = prot_seq.split("*")
                    longest = ""
                    for prot in prot_split:
                        if len(prot) > len(longest):
                            longest = prot
                    prot_seq = longest
                passes_length = False
                if USE_MIN_LEN:
                    if len(prot_seq) >= int(MIN_LEN):
                        passes_length = True
                else:
                    if len(prot_seq) >= PROP_SHORT * (len(prot_dic[reg.query])):
                        passes_length = True
                if passes_length:
                    success_file.write(">%s_%s_%s\n%s\n"% (species, family, gene_index, prot_seq))
                    success_cds.write(">%s_%s_%s\n%s\n"% (species, family, gene_index, cds_seq))
                    success_trans.write(">%s_%s_%s\n%s\n"% (species, family, gene_index, trans_seq))
                    gene_index += 1
                    success_file.flush()
                    success_trans.flush()
                    success_cds.flush()
                    cur_gff = open("%s/%s/%s_%s_%s_calib.gff" % (REGIONAL_DIR, species, region.scaf, region.start, region.end), 'rU')
                    for gff_line in cur_gff:
                        success_gff.write(gff_line)
                    success_gff.flush()
            
    success_file.close()
#    run_muscle(family, iteration, species)
    success_gff.close()
    success_cds.close()
    success_trans.close()
    #need to run "module add hmmer"
    
#    if family == "ir":
#        cmd = ["/Genomics/kocherlab/berubin/local/src/interproscan-5.21-60.0/interproscan.sh", "-i", "%s_genes_iter_%s/%s_successes_%s_%s.fa" % (family, iteration, species, family, MIN_LEN), "--seqtype", "p", "-T", "%s_genes_iter_%s/%s_successes_%s_%s_temp" % (family, iteration, species, family, MIN_LEN), "-b", "%s_genes_iter_%s/%s_successes_%s_%s.iprscan" % (family, iteration, species, family, MIN_LEN)]
#        subprocess.call(cmd)
#        print "IPR over"
#        pass_seqs = check_ir_gos("%s_genes_iter_%s/%s_successes_%s_%s.fa" % (family, iteration, species, family, MIN_LEN), "%s_genes_iter_%s/%s_successes_%s_%s.iprscan.tsv" % (family, iteration, species, family, MIN_LEN))
#        print pass_seqs
#        ir_file = open("%s_genes_iter_%s/%s_successes_%s_%s_iprpass.fa" % (family, iteration, species, family, MIN_LEN), 'w')
#        for gene, seq in pass_seqs.items():
#            ir_file.write(">%s\n%s\n" % (gene, seq))
#        ir_file.close()
                       

def run_muscle(family, iteration, species):
    outfile = open("%s/%s_AM_successes_%s_%s.faa" % (GENES_DIR, species, family, MIN_LEN), 'w')
    reader = SeqIO.parse("%s/%s_iter_%s_%s_%s_pep.faa" % (GENES_DIR, species, iteration, family, MIN_LEN), format = 'fasta')
    for rec in reader:
        outfile.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
    reader = SeqIO.parse("am_%s.fa" % family, format = 'fasta')
    for rec in reader:
        outfile.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
    outfile.close()
    cmd = ["muscle", "-in", "%s/%s_AM_successes_%s_%s.faa" % (GENES_DIR, species, family, MIN_LEN), "-out", "%s/%s_AM_successes_%s_%s.afaa" % (GENES_DIR, species, family, MIN_LEN)]
    subprocess.call(cmd)

def protein_length(region, species, family, iteration, prot_dic, genome_dic):
    if region.prot_len != -9:
        return region.prot_len
    ref_file = open("%s/%s/%s_%s_%s_%s_ref.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query), 'w')
    ref_file.write(">%s\n%s\n" % (region.query, prot_dic[region.query]))
    ref_file.close()

    outfile = open("%s/%s/%s_%s_%s_%s.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query), 'w')

    cur_start = region.start - FLANK_SIZE
    if cur_start < 0:
        cur_start = 0
    cur_end = region.end + FLANK_SIZE
    if cur_end > len(genome_dic[region.scaf]):
        cur_end = len(genome_dic[region.scaf])
    outfile.write(">%s_%s_%s_%s\n%s\n" % (region.scaf, region.start, region.end, region.query[0:20], genome_dic[region.scaf][cur_start:cur_end]))
    outfile.close()
    if iteration == 0:
        cmd = [EXONERATE_PATH, "--exhaustive", "y", "--forcegtag", "--model", "protein2genome", "--showtargetgff", "true", "--verbose", "0", "--showalignment", "no", "--showvulgar", "no", "--bestn", "1", "%s/%s/%s_%s_%s_%s_ref.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query), "%s/%s/%s_%s_%s_%s.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query)]
    else:
        cmd = [EXONERATE_PATH, "--exhaustive", "y", "--forcegtag", "--model", "protein2genome", "--showtargetgff", "true", "--verbose", "0", "--showalignment", "no", "--showvulgar", "no", "--bestn", "1", "%s/%s/%s_%s_%s_%s_ref.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query), "%s/%s/%s_%s_%s_%s.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query)]
    FNULL = open(os.devnull, 'w')
    with open("%s/%s/%s_%s_%s_%s.gff" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query), 'w') as outfile:
        subprocess.call(cmd, stdout = outfile, stderr=FNULL)
    outfile.close()
    target_path = "%s_%s_%s_%s" % (region.scaf, region.start, region.end, region.query)
    gene_id, prot_seq, cds_seq, align_score, new_start, new_end, trans_seq = reconstruct_gene(target_path, species, iteration, family, region.scaf, region.start, region.end, -9)
    if len(cds_seq) % 3 != 0:
        prot_seq = ""
    if "*" in prot_seq:
        prot_split = prot_seq.split("*")
        longest = ""
        for prot in prot_split:
            if len(prot) > len(longest):
                longest = prot
        prot_seq = longest
    if len(prot_seq) > 0:
        if prot_seq[0] == "M":
            align_score = align_score + 20

    os.remove("%s/%s/%s_reconst.fna" % (GENES_DIR, species, target_path))
    os.remove("%s/%s/%s_reconst.faa" % (GENES_DIR, species, target_path))
    os.remove("%s/%s/%s_%s_%s_%s_ref.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query))
    os.remove("%s/%s/%s_%s_%s_%s.fa" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query))
    os.remove("%s/%s/%s_calib.gff" % (REGIONAL_DIR, species, target_path))
    os.remove("%s/%s/%s_%s_%s_%s.gff" % (REGIONAL_DIR, species, region.scaf, region.start, region.end, region.query))
    return len(prot_seq), align_score, new_start, new_end

def check_ir_gos(seq_file, ipr_file):
    seq_dic = {}
    reader = SeqIO.parse(seq_file, format = 'fasta')
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    reader = open(ipr_file, 'rU')
    pass_list = []
    for line in reader:
        cur_line = line.split()
        cur_id = cur_line[0]
        if "PF10613" in cur_line or "PF00060" in cur_line:
            if cur_id not in pass_list:
                pass_list.append(cur_id)
    pass_seq_dic = {}
    for gene in pass_list:
        pass_seq_dic[gene] = seq_dic[gene]
    return pass_seq_dic

def merge_HSPs(my_hsp, temp_hsp_list):
    i = 0
    included = False
    border_size = INTRON_SIZE
    query_overlap = HIT_OVERLAP

    while i < len(temp_hsp_list):
        cur_hsp = temp_hsp_list[i]
        if my_hsp.scaf != cur_hsp.scaf:
            i += 1
            continue
        if my_hsp.frame != cur_hsp.frame:
            i += 1
            continue
        if my_hsp.frame > 0:
            if my_hsp.sstart <= cur_hsp.sstart and my_hsp.send >= cur_hsp.sstart - border_size:
                if my_hsp.send >= cur_hsp.send:
                    new_hsp = my_hsp
                    my_hsp = new_hsp
                    del temp_hsp_list[i]
                    i = -1
                else:
                    if my_hsp.qend < cur_hsp.qstart + query_overlap:
                        new_hsp = HSP(my_hsp.qstart, cur_hsp.qend, my_hsp.sstart, cur_hsp.send, my_hsp.scaf, my_hsp.frame, my_hsp.bits + cur_hsp.bits)
                        my_hsp = new_hsp
                        del temp_hsp_list[i]
                        i = -1
            elif my_hsp.send >= cur_hsp.send and my_hsp.sstart <= cur_hsp.send + border_size:

                if my_hsp.sstart >= cur_hsp.sstart:
                    if cur_hsp.qend < my_hsp.qstart + query_overlap:
                        new_hsp = HSP(cur_hsp.qstart, my_hsp.qend, cur_hsp.sstart, my_hsp.send, my_hsp.scaf, my_hsp.frame, my_hsp.bits + cur_hsp.bits)
                        my_hsp = new_hsp
                        del temp_hsp_list[i]
                        i = -1

        elif my_hsp.frame < 0:
            if my_hsp.sstart <= cur_hsp.sstart and my_hsp.send >= cur_hsp.sstart - border_size:
                if my_hsp.send >= cur_hsp.send:
                    new_hsp = my_hsp
                    my_hsp = new_hsp
                    del temp_hsp_list[i]
                    i = -1
                else:
                    if my_hsp.qstart > cur_hsp.qend - query_overlap:
                        new_hsp = HSP(cur_hsp.qstart, my_hsp.qend, my_hsp.sstart, cur_hsp.send, my_hsp.scaf, my_hsp.frame, my_hsp.bits + cur_hsp.bits)
                        my_hsp = new_hsp
                        del temp_hsp_list[i]
                        i = -1
            elif my_hsp.send >= cur_hsp.send and my_hsp.sstart <= cur_hsp.send + border_size:
                if my_hsp.sstart >= cur_hsp.sstart:
                    if cur_hsp.qstart > my_hsp.qend - query_overlap:
                        new_hsp = HSP(my_hsp.qstart, cur_hsp.qend, cur_hsp.sstart, my_hsp.send, my_hsp.scaf, my_hsp.frame, my_hsp.bits + cur_hsp.bits)
                        my_hsp = new_hsp
                        del temp_hsp_list[i]
                        i = -1
        elif my_hsp.sstart >= cur_hsp.sstart and my_hsp.send <= cur_hsp.send:
            included = True
            break
        i += 1
    if not included:
        temp_hsp_list.append(my_hsp)
    return temp_hsp_list

def reconstruct_gene(target_path, species, iteration, family, scaf, orig_start, end, gene_index):
    reader = SeqIO.parse("%s/%s/%s.fa" % (REGIONAL_DIR, species, target_path), format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
        cur_seq = str(rec.seq)
    reader = open("%s/%s/%s.gff" % (REGIONAL_DIR, species, target_path), 'rU')
    seq_str = ""
    gene_counter = 1
    model_dic = {}
    score = 0
    new_start = 0
    new_end = 0
    gene_start = 0
    gene_end = 0
    trans_seq = ""
    trans_dic = {}
    new_gff = open("%s/%s/%s_calib.gff" % (REGIONAL_DIR, species, target_path), 'w')
    if orig_start < FLANK_SIZE:
        cur_flank_size = orig_start
    else:
        cur_flank_size = FLANK_SIZE
    start = orig_start
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.split()
        if cur_line[2] == "gene":
            cur_gene = gene_counter
            model_dic[cur_gene] = ""
            gene_counter += 1
            score = int(cur_line[5])
            new_start = int(cur_line[3])
            new_end = int(cur_line[4])
            gene_start = start + new_start - cur_flank_size
            gene_end = start + new_end - cur_flank_size
            cur_strand = cur_line[6]
            if cur_strand == "+":
                trans_dic[cur_gene] = cur_seq[new_start-1:new_end]
            else:
                trans_dic[cur_gene] = str(Seq.Seq(cur_seq[new_start-1:new_end]).reverse_complement())

            new_gff.write("%s\tprotein2genome\tgene\t%s\t%s\t.\t%s\t.\tID=%s_%s_%s;Name=%s_%s_%s;\n" % (scaf, start + new_start - cur_flank_size, start + new_end - cur_flank_size, cur_strand, species, family, gene_index, species, family, gene_index))
            new_gff.write("%s\tprotein2genome\tmRNA\t%s\t%s\t.\t%s\t.\tID=%s_%s_%s-RA;Parent=%s_%s_%s;Name=%s_%s_%s-RA;\n" % (scaf, start + new_start - cur_flank_size, start + new_end - cur_flank_size, cur_strand, species, family, gene_index, species, family, gene_index,species, family, gene_index))
        if cur_line[2] == "cds":
            cur_strand = cur_line[6]
            cur_start = int(cur_line[3])
            cur_end = int(cur_line[4])
            new_gff.write("%s\tprotein2genome\tCDS\t%s\t%s\t.\t%s\t.\tID=%s_%s_%s-RA:cds;Parent=%s_%s_%s-RA;\n" % (scaf, start + cur_start - cur_flank_size , start + cur_end - cur_flank_size, cur_strand, species, family, gene_index, species, family, gene_index))
            new_gff.write("%s\tprotein2genome\texon\t%s\t%s\t.\t%s\t.\tID=%s_%s_%s-RA;Parent=%s_%s_%s-RA;\n" % (scaf, start + cur_start - cur_flank_size, start + cur_end - cur_flank_size, cur_strand, species, family, gene_index, species, family, gene_index))
            if cur_strand == "+":
                model_dic[cur_gene] = model_dic[cur_gene] + cur_seq[cur_start-1:cur_end]
            else:
                model_dic[cur_gene] = model_dic[cur_gene] + str(Seq.Seq(cur_seq[cur_start-1:cur_end]).reverse_complement())
    longest = 0
    longest_index = -1
    for gene_index, gene_seq in model_dic.items():
        if len(gene_seq) > longest:
            seq_str = gene_seq
            longest = len(gene_seq)
            trans_seq = trans_dic[gene_index]
    if not os.path.isdir("%s/%s" % (GENES_DIR, species)):
        os.mkdir("%s/%s" % (GENES_DIR, species))
    new_gff.close()
    outfile = open("%s/%s/%s_reconst.fna" % (GENES_DIR, species, target_path), 'w')
    outfile.write(">%s\n%s\n" % (rec.id, seq_str))
    outfile.close()
    outfile = open("%s/%s/%s_reconst.faa" % (GENES_DIR, species, target_path), 'w')
    outfile.write(">%s\n%s\n" % (rec.id, str(Seq.Seq(seq_str).translate())))
    outfile.close()
    return rec.id, str(Seq.Seq(seq_str).translate()), seq_str, score, gene_start, gene_end, trans_seq

def sort_hsps(hsp_list):
    new_hsp_list = []
    hsp_count = len(hsp_list)
    while True:
        smallest_hsp = hsp_list[0]
        for hsp in hsp_list:
            if hsp.less_than(smallest_hsp):
                smallest_hsp = hsp
        new_hsp_list.append(smallest_hsp)
        hsp_list.remove(smallest_hsp)
        if len(new_hsp_list) == hsp_count:
            break
    return new_hsp_list

class GeneRegion:
    def __init__(self, start, end, scaf, frame, bits, query):
        self.start = start
        self.end = end
        self.scaf = scaf
        self.frame = frame
        self.bits = bits
        self.query = query
        self.prot_len = -9
        self.align_score = -9
    
    def overlap(self, other_region):
        if self.scaf != other_region.scaf:
            return False
        if self.start <= other_region.start and self.end >= other_region.start:
            return True
        if self.start >= other_region.start and self.start <= other_region.end:
            return True
        if self.start >= other_region.start and self.end <= other_region.end:
            return True
        return False

    def change_start(self, new_start):
        self.start = new_start

    def equals(self, other_region):
        if self.start == other_region.start and self.end == other_region.end and self.scaf == other_region.scaf:
            return True

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (self.start, self.end, self.scaf, self.prot_len, self.end - self.start, self.query)
            

class HSP:
    def __init__(self, qstart, qend, sstart, send, scaf, frame, bits):
        if frame < 0:
            self.frame = -1
        else:
            self.frame = 1
        if sstart > send:
            self.sstart = send
            self.send = sstart
        else:
            self.sstart = sstart
            self.send = send
        self.qstart = qstart
        self.qend = qend
        self.scaf = scaf
        self.bits = bits
        
    def less_than(self, other_hsp):
        if self.sstart < other_hsp.sstart:
            return True
        return False

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (self.sstart, self.send, self.qstart, self.qend, self.scaf, self.bits)


if __name__ == '__main__':
    main()
