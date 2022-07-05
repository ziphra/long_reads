#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, pprint, gzip, vcf, argparse, subprocess, tabix, time, csv
from pyfaidx import Fasta

""" Author: Julien BURATTI
	julien.buratti@aphp.fr
	Hopital Pitie-Salpetriere, AP-HP """

pp = pprint.PrettyPrinter(indent=4)

########## PATHS ##########
# Gene based
omim_f = "/media/euphrasie/Elements/Annotations/OMIM/apr2022/genemap2.txt"
sysid_f = "/media/euphrasie/Elements/Annotations/sysID/apr2022/disease-info.csv"
masked_genes_f = "/media/euphrasie/Elements/Annotations/genesToMask/genes_a_masquer_Dec2020.txt"
clinvar_gene_f = "/media/euphrasie/Elements/Annotations/clinvar/apr2022/gene_condition_source_id"
hgmd_f = "/media/euphrasie/Elements/Annotations/HGMD/HGMD_2021.4/hgmd_pro_2021.4_"
acc_f = "/media/euphrasie/Elements/Annotations/genesACC/genesACC_Aug2019.txt"
loeuf_f = "/media/euphrasie/Elements/Annotations/LOEUF/gnomad.v2.1.1.lof_metrics.by_gene.txt"
metabo_f = "/media/euphrasie/Elements/Annotations/Metabo/metabo_and_onco_may2022.txt"
# Variant based
clinvar_variant_dir = "/media/euphrasie/Elements/Annotations/clinvar/apr2022/"
mistic_dir = "/media/euphrasie/Elements/Annotations/Mistic/"
gnomad_dir = "/media/euphrasie/Elements/Annotations/gnomAD_merge_exom2.1.1_genom3.1.2/"

########## FUNCTIONS ##########

def read_vcf(vep_vcf):
	vcf_r = vcf.Reader(open(vep_vcf, 'rb'))
	return vcf_r

def write_vcf(vcf_r, target_chrom):
	vcf_w = vcf.Writer(open('vep_and_custom_annotations_chr' + target_chrom + '.vcf', 'w'), vcf_r)
	return vcf_w

def replace_chars(string):
	new_string = string.replace('\r','').replace('\n','').replace(';','.').replace('=',':').replace(' ', '_')
	return new_string

def add_info_to_header(vep_vcf):
	with gzip.open(vep_vcf, 'rb') as source, open('header.vcf', 'w') as destination:
		for line in source:
			line = line.decode("utf-8").rstrip()
			if line.startswith('##'):
				destination.write(line + "\n")
			if line.startswith('#CHROM'):
				chrom_line = line
		for string in 'omim_inheritance', 'sysID_group', 'sysID_disease', 'sysID_synopsis', 'clinvar_gene', 'clinvar_id', 'clinvar_clnsig', 'hgmd_gene', 'hgmd_id', 'hgmd_class', 'hgmd_phen', 'flankseq', 'gnomad_flag', 'appris', 'metabo_gene', 'mistic_pred':
			destination.write("##INFO=<ID=" + string.upper().replace('_CSQ','') + ",Number=1,Type=String,Description=\"" + string + "\">\n")
		for integer in 'gene2mask', 'acc_gene', 'gnomad_homhem', 'gnomad_het', 'lab_htz', 'lab_htz_adj', 'lab_hmz', 'lab_hmz_adj', 'masked_genes':
			destination.write("##INFO=<ID=" + integer.upper().replace('_CSQ','') + ",Number=1,Type=Integer,Description=\"" + integer + "\">\n")
		for floating_nb in 'loeuf', 'cadd_phred', 'mistic_score', 'gnomad_AF':
			destination.write("##INFO=<ID=" + floating_nb.upper().replace('_CSQ','') + ",Number=1,Type=Float,Description=\"" + floating_nb + "\">\n")
		destination.write(chrom_line + "\n")

def omim_inheritance():
	print("  + OMIM")
	omim_dict = {}
	with open(omim_f, 'r') as omim:
		for line in omim:
			AD = False
			AR = False
			XL = False
			XLD = False
			XLR = False
			YL = False
			inheritance = ""
			if line.startswith("chr"):
				gene_and_aliases = line.split("\t")[6]
				description = line.split("\t")[7] + "\t" + line.split("\t")[12]
				if line.startswith("chrX"):
					if "X-linked" in description:
						if not "recessive" in description and not "dominant" in description:
								XL = True
						if "dominant" in description:
							XLD = True
						if "recessive" in description:
							XLR = True
					if XL == True and XLD == False and XLR == False:
						inheritance = "XL"
					if XLD == True and XLR == False:
						inheritance = "XLD"
					if XLD == False and XLR == True:
						inheritance = "XLR"
					if XLD == True and XLR == True:
						inheritance = "XLD/XLR"
				elif line.startswith("chrY"):
					if "Y-linked" in description:
						YL = True
					if YL == True:
						inheritance = "YL"
				else:
					if "utosomal dominant" in description:
						AD = True
					if "utosomal recessive" in description:
						AR = True
					if AD == True and AR == False:
						inheritance = "AD"
					if AD == False and AR == True:
						inheritance = "AR"
					if AD == True and AR == True:
						inheritance = "AD/AR"
				for symbol in gene_and_aliases.split(","):
					symbol = symbol.strip().upper()
					if inheritance != "":
						omim_dict[symbol.upper()] = inheritance
	return omim_dict

def sysID():
	print("  + sysID")
	sysID_dict = {}
	with open(sysid_f, 'r') as csvfile:
		csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
		# disease-info # Entrez id	Gene symbol	Gene group	Inheritance pattern	Inheritance type	Main class 	Accompanying phenotype	Limited confidence criterion	Sysid yes no	Disease subtype	Disease type	Alternative names	Additional references	Omim disease	Gene review	Haploinsufficiency yes no	Clinical synopsis	Human gene disease remark 
		for line in csv_reader:
			if "Entrez id" in line:
				symbol_idx = line.index("Gene symbol")
				group_idx = line.index("Gene group")
				# description_idx = line.index("Gene description")
				# inheritPattern_idx = line.index("Inheritance pattern")
				# inheritType_idx = line.index("Inheritance type")
				disease_idx = line.index("Disease type")
				synopsis_idx = line.index("Clinical synopsis")
			else:
				group_final = ""
				symbol = line[symbol_idx]
				group = line[group_idx]
				disease = "\"" + line[disease_idx].replace(" ","_").replace(",_",",").replace(":_",":") + "\""
				synopsis = "\"" + line[synopsis_idx].replace(" ","_").replace(",_",",").replace(":_",":") + "\""

				if group == "Negative control ID screen":
					group_final = "0"
				if "ID data freeze" in group or "urrent primary ID" in group:
					group_final = "V" # validé
				elif "ID candidate" in group:
					group_final = "C" # candidat

				# print(symbol + " | " + group + " (" + group_final + ") | Disease = " + disease + " | Synopsis = " + synopsis)

				if symbol not in sysID_dict:
					sysID_dict[symbol] = {}
					sysID_dict[symbol]['group'] = group_final
					sysID_dict[symbol]['disease'] = []
					sysID_dict[symbol]['synopsis'] = []
					sysID_dict[symbol]['disease'].append(disease)
					sysID_dict[symbol]['synopsis'].append(synopsis)
				else:
					if disease not in sysID_dict[symbol]['disease']:
						sysID_dict[symbol]['disease'].append(disease)
					if synopsis not in sysID_dict[symbol]['synopsis']:
						sysID_dict[symbol]['synopsis'].append(synopsis)
	return sysID_dict

def masked_genes():
	print("  + Masked genes")
	masked_genes_dict = {}
	with open(masked_genes_f, "r") as masked:
		for line in masked:
			line = line.rstrip()
			if line not in masked_genes_dict:
				masked_genes_dict[line] = 1
	return masked_genes_dict

def clinvar_gene():
	print("  + Clinvar genes")
	clinvar_gene_dict = {}
	with open(clinvar_gene_f, 'r') as clinvar_g:
		for line in clinvar_g:
			line = line.rstrip()
			if not line.startswith('#'):
				clinvar_disease = "\"" + line.split('\t')[4].replace(" ","_") + "\""
				if line.split('\t')[1] != "":
					clinvar_gene_symbol = line.split('\t')[1]
				else:
					clinvar_gene_symbol = line.split('\t')[2]
				if clinvar_gene_symbol not in clinvar_gene_dict:
					clinvar_gene_dict[clinvar_gene_symbol] = clinvar_disease
				else:
					if clinvar_disease not in clinvar_gene_dict[clinvar_gene_symbol]:
						clinvar_gene_dict[clinvar_gene_symbol] += "," + clinvar_disease
	return clinvar_gene_dict

def gnomad(target_chrom, vcf_string):
	gnomad_tb = tabix.open(gnomad_dir + "merge_exom_genom_chr" + target_chrom + ".vcf.gz")
	chrom = vcf_string.split('-')[0]
	pos = vcf_string.split('-')[1]
	ref = vcf_string.split('-')[2]
	alt = vcf_string.split('-')[3]
	gnomad_key = 'chr' + chrom + ":" + pos + "-" + pos
	gnomad_records = gnomad_tb.querys(gnomad_key)
	gnomad_flag = "."
	gnomad_AF = "."
	gnomad_het = "."
	gnomad_homhem = "."
	for r in gnomad_records:
		gnomad_ref = r[3]
		gnomad_alt = r[4]
		gnomad_flag = r[6]
		gnomad_info = r[7]
		if (gnomad_ref == ref) and (gnomad_alt == alt):
			for i in gnomad_info.split(';'):
				# if i.startswith("AC="): gnomad_AC = i.split('=')[1]
				# if i.startswith("AN="): gnomad_AN = i.split('=')[1]
				if i.startswith("AF="): gnomad_AF = i.split('=')[1]
				if i.startswith("het="): gnomad_het = i.split('=')[1]
				if i.startswith("homhem="): gnomad_homhem = i.split('=')[1]
			break
	return gnomad_flag, gnomad_AF, gnomad_het, gnomad_homhem

def clinvar_variant(genome, target_chrom):
	print("  + Clinvar variants")
	if genome == "hg19":
		release = "GRCh37"
	else:
		release = "GRCh38"
	clinvar_variant_f = clinvar_variant_dir + release + "/splitChr/clinvar_hg38_ascii_chr" + target_chrom + ".vcf.gz"
	clinvar_variant_dict = {}
	vcf_r = vcf.Reader(open(clinvar_variant_f, 'rb'))
	for variant in vcf_r:
		if len(variant.ALT) > 0:
			k = variant.CHROM + "-" + str(variant.end) + "-" + variant.REF + "-" + str(variant.ALT[0])
			if k not in clinvar_variant_dict:
				clinvar_variant_dict[k] = {}
				clinvar_variant_dict[k]['ID'] = variant.ID
				if "CLNSIG" in variant.INFO:
					if "Conflicting_interpretations_of_pathogenicity" in variant.INFO['CLNSIG']:
						if "CLNSIGCONF" in variant.INFO:
							clinvar_variant_dict[k]['CLNSIG'] = ','.join(variant.INFO['CLNSIG']) + ":" + ','.join(variant.INFO['CLNSIGCONF'])
						else:
							clinvar_variant_dict[k]['CLNSIG'] = ','.join(variant.INFO['CLNSIG'])
					else:
						clinvar_variant_dict[k]['CLNSIG'] = ','.join(variant.INFO['CLNSIG'])
				elif "CLNSIGINCL" in variant.INFO:
					clinvar_variant_dict[k]['CLNSIG'] = str(variant.INFO['CLNSIGINCL'])
	return clinvar_variant_dict

def hgmd(genome):
	print("  + HGMD")
	hgmd_dict_v = {}
	hgmd_dict_g = {}
	vcf_r = vcf.Reader(open(hgmd_f + genome + ".vcf", 'r'))
	for variant in vcf_r:
		vcf_string = variant.CHROM + "-" + str(variant.POS) + "-" + variant.REF + "-" + str(variant.ALT[0])
		hgmd_dict_v[vcf_string] = {}
		hgmd_dict_v[vcf_string]['ID'] = variant.ID
		if 'CLASS' in variant.INFO:
			hgmd_dict_v[vcf_string]['CLASS'] = variant.INFO['CLASS']
		else:
			hgmd_dict_v[vcf_string]['CLASS'] = '.'
		hgmd_phen = variant.INFO['PHEN'].replace("%2C",",").replace(",_",",")
		hgmd_dict_v[vcf_string]['PHEN'] = hgmd_phen

		if 'GENE' in variant.INFO:
			if variant.INFO['GENE'] not in hgmd_dict_g:
				hgmd_dict_g[variant.INFO['GENE']] = hgmd_phen
			else:
				if hgmd_phen not in hgmd_dict_g[variant.INFO['GENE']]:
					newValue = hgmd_dict_g[variant.INFO['GENE']] + "," + hgmd_phen
					hgmd_dict_g[variant.INFO['GENE']] = newValue
	return hgmd_dict_g, hgmd_dict_v

def acc():
	print("  + ACC genes")
	acc_dict = {}
	with open(acc_f, "r") as acc:
		for line in acc:
			line = line.rstrip()
			if line not in acc_dict:
				acc_dict[line] = 1
	return acc_dict

def loeuf():
	print("  + pLoF scores")
	loeuf_dict = {}
	with open(loeuf_f,'r') as lof:
		for line in lof:
			if not line.startswith("gene\ttranscript"):
				""" loeuf_dict[gene] = oe_lof_upper """
				loeuf_dict[line.split("\t")[0]] = line.split("\t")[29]
	return loeuf_dict

def metabo():
	print("  + Metabo genes")
	metabo_dict = {}
	with open("/media/euphrasie/Elements/Annotations/Metabo/metabo_and_onco_feb2022.txt", "r") as metabo_g:
		for line in metabo_g:
			line = line.rstrip()
			gene = line.split("\t")[0]
			panel = line.split("\t")[1]
			if gene not in metabo_dict:
				metabo_dict[gene] = panel
	return metabo_dict

def flank_seq(fasta_r, vcf_string):
	chrom = 'chr' + vcf_string.split('-')[0]
	pos = int(vcf_string.split('-')[1])
	ref = vcf_string.split('-')[2]
	alt = vcf_string.split('-')[3]
	flank_left, flank_right = "",""
	flank_left = fasta_r[chrom][int(pos)-13:int(pos)-1].seq
	flank_right = fasta_r[chrom][int(pos):int(pos)+12].seq

	flankseq = flank_left + "[" + ref + "/" + alt + "]" + flank_right
	return flankseq

def mistic_d(target_chrom):
	print("  + Mistic score")
	if genome == "hg19":
		release = "37"
	else:
		release = "38"
	mistic_path = mistic_dir + "/GRCh" + release + "_splitChr/MISTIC_GRCh" + release + "_chr" + target_chrom +".tsv.gz"
	# start_time = time.time()
	# print(chrom_prev)
	mistic_dict = {}
	if target_chrom != "Y":
		with gzip.open(mistic_path, 'rb') as m:
			for line in m:
				line = line.decode("utf-8").rstrip()
				if not line.startswith("#"):
					chrom = line.split("\t")[0]
					pos = line.split("\t")[1]
					ref = line.split("\t")[2]
					alt = line.split("\t")[3]
					score = line.split("\t")[4]
					pred = line.split("\t")[5]
					mistic_dict[chrom + "-" + pos + "-" + ref + "-" + alt] = score + "|" + pred
	# print("--- elapsed time: %s minutes ---" % ((time.time() - start_time) / 60.0))
	return mistic_dict

def mistic_t(genome, mistic_dir, chrom_v, pos_v, ref_v, alt_v):
	if genome == "hg19":
		release = "37"
	else:
		release = "38"
	mistic_path = mistic_dir + "MISTIC_GRCh" + release + "_sorted.tsv.gz"
	mistic_chrom_l = subprocess.check_output(['tabix', '-l', mistic_path])
	mistic_chrom_l = mistic_chrom_l.decode("utf-8").split('\n')[:-1]
	mistic_score, mistic_pred = ".", "."
	mistic_key = chrom_v + ":" + str(pos_v) + "-" + str(pos_v)
	if chrom_v in mistic_chrom_l:
		mistic_tb = tabix.open(mistic_path)
		mistic_records = mistic_tb.querys(mistic_key)
		for r in mistic_records:
			mistic_file_ref = r[2]
			mistic_file_alt = r[3]
			if mistic_file_ref == ref_v and mistic_file_alt == alt_v:
				mistic_score = r[4]
				mistic_pred = r[5]
	return mistic_score, mistic_pred

########## MAIN ##########

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", help="hg19/hg38", choices=["hg19","hg38"])
	parser.add_argument("-c", help="select chromosome", choices=["1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"])
	args = parser.parse_args()
	genome = args.g
	target_chrom = args.c
	print("\nReference Genome = " + genome)
	print("Target Chromosome = chr" + target_chrom)
	if genome == "hg19":
		genome_fasta = '/media/euphrasie/Elements/ref_genomes/hg19/hg19_std_M-rCRS_Y-PAR-mask.fa'
		chrom_fasta = '/media/euphrasie/Elements/ref_genomes/hg19/splitChr/hg19_GenDev_chr' + target_chrom + ".fa"
	else:
		genome_fasta = '/media/euphrasie/DATA/reference_genome/newhg38/GRCh38.primary_assembly.genome.fa'
		chrom_fasta = '/media/euphrasie/DATA/reference_genome/newhg38/splitChr/GRCh38.primary_assembly.chr' + target_chrom + '.fa'

	""" Copy VEP VCF file """
	subprocess.call('zgrep \"^#\" $VEP_LOCAL_VCF > header', shell=True)
	subprocess.call("zgrep ^chr" + target_chrom + "\'\t\' ${VEP_LOCAL_VCF} > vep_chr" + target_chrom + ".txt", shell=True)
	subprocess.call("cat header vep_chr" + target_chrom + ".txt | bgzip > vep_chr" + target_chrom + ".vcf.gz", shell=True)
	subprocess.call("tabix vep_chr" + target_chrom + ".vcf.gz", shell=True)
	subprocess.call("rm -f vep_chr" + target_chrom + ".txt", shell=True)


	""" Prepare VCF files """
	vep_vcf = "vep_chr" + target_chrom + ".vcf.gz"
	vcf_r = read_vcf(vep_vcf)
	vcf_w = write_vcf(vcf_r, target_chrom)
	add_info_to_header(vep_vcf)
	""" Load Annotations """
	print("\n*** Loading annotations ***")
	""" OMIM """
	omim_d = omim_inheritance()
	""" sysID """
	sysid_d = sysID()
	""" Masked Genes """
	masked_genes_d = masked_genes()
	""" Clinvar Genes """
	clinvar_gene_d = clinvar_gene()
	""" Clinvar Variants """
	clinvar_variant_d = clinvar_variant(genome, target_chrom)
	""" ACC genes """
	acc_d = acc()
	""" pLoF """
	loeuf_d = loeuf()
	""" Metabo """
	metabo_d = metabo()
	""" HGMD """
	hgmd_gene_d, hgmd_variant_d = hgmd(genome)
	""" Mistic """
	mistic_dict = mistic_d(target_chrom)
	""" Flankseq """
	print("  + Indexing genome fasta")
	fasta_r = Fasta(chrom_fasta)
	print("*** Annotations loaded ***")

	""" VEP annotations """
	print("\n## Reading variants, annotate and write new VCF file ##")
	for variant in vcf_r:
		chrom_v = variant.CHROM.replace('chr', '')
		if chrom_v == target_chrom:
			pos_v = variant.POS
			id_v = variant.ID
			ref_v = variant.REF
			alt_v = repr(variant.ALT[0])
			qual_v = variant.QUAL
			filter_v = variant.FILTER
			info_v = variant.INFO
			format_v = variant.FORMAT
			end_v = pos_v + (len(ref_v)-1)
			vcf_string = chrom_v + '-' + str(pos_v) + '-' + ref_v + '-' + alt_v
			# print(variant)

			# VEP = Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|HGVS_OFFSET|HGVSg|CLIN_SIG|SOMATIC|PHENO|PUBMED|existing_InFrame_oORFs|existing_OutOfFrame_oORFs|existing_uORFs|five_prime_UTR_variant_annotation|five_prime_UTR_variant_consequence|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|REVEL|CADD_PHRED|CADD_RAW"
			vep_csq = variant.INFO['CSQ'][0]

			vep_allele = vep_csq.split('|')[0]
			vep_consequence = vep_csq.split('|')[1]
			vep_impact = vep_csq.split('|')[2]
			vep_symbol = vep_csq.split('|')[3]
			vep_gene = vep_csq.split('|')[4]
			vep_feature_type = vep_csq.split('|')[5]
			vep_feature = vep_csq.split('|')[6]
			vep_biotype = vep_csq.split('|')[7]
			vep_exon = vep_csq.split('|')[8]
			vep_intron = vep_csq.split('|')[9]
			vep_hgvsc = vep_csq.split('|')[10]
			vep_hgvsp = vep_csq.split('|')[11]
			vep_cdna_position = vep_csq.split('|')[12]
			vep_cds_position = vep_csq.split('|')[13]
			vep_protein_position = vep_csq.split('|')[14]
			vep_amino_acids = vep_csq.split('|')[15]
			vep_codons = vep_csq.split('|')[16]
			vep_existing_variation = vep_csq.split('|')[17]
			vep_distance = vep_csq.split('|')[18]
			vep_strand = vep_csq.split('|')[19]
			vep_flags = vep_csq.split('|')[20]
			vep_variant_class = vep_csq.split('|')[21]
			vep_symbol_source = vep_csq.split('|')[22]
			vep_hgnc_id = vep_csq.split('|')[23]
			vep_canonical = vep_csq.split('|')[24]
			vep_mane_select = vep_csq.split('|')[25]
			vep_mane_plus_clinical = vep_csq.split('|')[26]
			vep_tsl = vep_csq.split('|')[27]
			vep_appris = vep_csq.split('|')[28]
			vep_ccds = vep_csq.split('|')[29]
			vep_hgvs_offset = vep_csq.split('|')[30]
			vep_hgvsg = vep_csq.split('|')[31]
			vep_clin_sig = vep_csq.split('|')[32]
			vep_somatic = vep_csq.split('|')[33]
			vep_pheno = vep_csq.split('|')[34]
			vep_pubmed = vep_csq.split('|')[35]
			vep_domains = vep_csq.split('|')[33]
			vep_existing_inframe_oorfs = vep_csq.split('|')[36]
			vep_existing_outframe_oorfs = vep_csq.split('|')[37]
			vep_existing_uorfs = vep_csq.split('|')[38]
			vep_five_prime_utr_variant_annotation = vep_csq.split('|')[39]
			vep_five_prime_utr_variant_consequence = vep_csq.split('|')[40]
			vep_spliceai_pred_dp_ag = vep_csq.split('|')[41]
			vep_spliceai_pred_dp_al = vep_csq.split('|')[42]
			vep_spliceai_pred_dp_dg = vep_csq.split('|')[43]
			vep_spliceai_pred_dp_dl = vep_csq.split('|')[44]
			vep_spliceai_pred_ds_ag = vep_csq.split('|')[45]
			vep_spliceai_pred_ds_al = vep_csq.split('|')[46]
			vep_spliceai_pred_ds_dg = vep_csq.split('|')[47]
			vep_spliceai_pred_ds_dl = vep_csq.split('|')[48]
			vep_spliceai_pred_symbol = vep_csq.split('|')[49]
			vep_revel = vep_csq.split('|')[50]
			vep_cadd_phred = vep_csq.split('|')[51]
			vep_cadd_raw = vep_csq.split('|')[52]

			""" Custom annotation for each variant """

			""" Clinvar Gene """
			clinvar_gene = "."
			if vep_symbol in clinvar_gene_d:
				clinvar_gene = clinvar_gene_d[vep_symbol]
			""" Clinvar Variant """
			clinvar_clnsig = "."
			clinvar_id = "."
			if vcf_string in clinvar_variant_d:
				clinvar_id = clinvar_variant_d[vcf_string]['ID']
				clinvar_clnsig = clinvar_variant_d[vcf_string]['CLNSIG']
			""" Transmission OMIM du gene """
			omim_inheritance = "."
			if vep_symbol in omim_d:
				omim_inheritance = omim_d[vep_symbol]
			""" Genes to mask """
			masked_genes = 0
			if vep_symbol in masked_genes_d:
				masked_genes = 1
			""" sysID """
			sysID_group, sysID_disease, sysID_synopsis = ".", ".", "."
			if vep_symbol in sysid_d:
				sysID_group = sysid_d[vep_symbol]['group']
				sysID_disease = "\"" + replace_chars(','.join(sysid_d[vep_symbol]['disease'])) + "\""
				sysID_synopsis = "\"" + replace_chars(','.join(sysid_d[vep_symbol]['synopsis'])) + "\""
			""" ACC """
			acc_gene = "0"
			if vep_symbol in acc_d:
				acc_gene = acc_d[vep_symbol]
			""" METABO """
			metabo_gene = "."
			if vep_symbol in metabo_d:
				metabo_gene = metabo_d[vep_symbol]
			""" HGMD """
			hgmd_gene = "."
			hgmd_id, hgmd_class, hgmd_phen = ".", ".", "."
			if vep_symbol in hgmd_gene_d:
				hgmd_gene = hgmd_gene_d[vep_symbol]
			if vcf_string in hgmd_variant_d:
				hgmd_id = hgmd_variant_d[vcf_string]['ID']
				hgmd_class = hgmd_variant_d[vcf_string]['CLASS']
				hgmd_phen = hgmd_variant_d[vcf_string]['PHEN']
			""" MISTIC """
			mistic_score = "."
			mistic_pred = "."
			if vcf_string in mistic_dict:
				mistic_score = mistic_dict[vcf_string].split('|')[0]
				mistic_pred = mistic_dict[vcf_string].split('|')[1]
			""" gnomAD hg38 """
			gnomad_flag, gnomad_AF, gnomad_het, gnomad_homhem = gnomad(target_chrom, vcf_string)
			""" gene pLoF """
			if vep_symbol in loeuf_d:
				if loeuf_d[vep_symbol] != "NA":
					loeuf = loeuf_d[vep_symbol]
				else:
					loeuf = ""
			""" Flankseq """
			flankseq = flank_seq(fasta_r, vcf_string)

			""" Ecriture des resultats dans un nouveau VCF """
			for i in 'omim_inheritance', 'masked_genes', 'gnomad_flag', 'gnomad_AF', 'gnomad_het', 'gnomad_homhem', 'sysID_group', 'sysID_disease', 'sysID_synopsis', 'acc_gene', 'clinvar_gene', 'clinvar_id', 'clinvar_clnsig', 'hgmd_gene', 'hgmd_id', 'hgmd_class', 'hgmd_phen', 'flankseq', 'loeuf', 'mistic_score', 'mistic_pred', 'metabo_gene': #'omim_inheritance', 'masked_genes', 'sysID_group', 'sysID_disease', 'sysID_synopsis', 'acc_gene', 'clinvar_gene', 'clinvar_id', 'clinvar_clnsig', 'hgmd_gene', 'hgmd_id', 'hgmd_class', 'hgmd_phen', 'flankseq', 'loeuf', 'mistic_score', 'mistic_pred', 'metabo_gene':
				value = replace_chars(str(locals()[i]))
				if value == "":
					value = "."
				variant.INFO[i.upper()] = value
			vcf_w.write_record(variant)

########## EOF ##########