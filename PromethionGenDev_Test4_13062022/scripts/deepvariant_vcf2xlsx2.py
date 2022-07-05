#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, argparse, pprint, xlsxwriter
from cyvcf2 import VCF

""" Author: Julien BURATTI
	julien.buratti@aphp.fr
	Hopital Pitie-Salpetriere, AP-HP """

pp = pprint.PrettyPrinter(indent=4)

parser = argparse.ArgumentParser()
#parser.add_argument("-d", help="Dad ID")
#parser.add_argument("-m", help="Mom ID")
parser.add_argument("-s", help="sample ID")

args = parser.parse_args()
#dad_id = args.d
#mom_id = args.m
#sample_id = args.c

vcf_file = os.environ["VCF"]

vcf_r = VCF(vcf_file)
samples = vcf_r.samples

#for i in dad_id, mom_id, sample_id:
#	if i == None:
#		sys.exit("Error: 3 arguments expected.")
#	elif i not in samples:
#		sys.exit("Error: Check samples names")

#dad_idx = samples.index(dad_id)
#mom_idx = samples.index(mom_id)
#sample = samples.index(sample_id)

sample = "PromethionGenDev_Test4_13062022"

# with open(sample_id + "_" + mom_id + "_" + dad_id + ".txt", "w") as output:
# Create a workbook and add a worksheet.
workbook = xlsxwriter.Workbook(sample+ '_annot.xlsx.zip')
workbook.use_zip64()

worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
three_decimals = workbook.add_format({'num_format': '0.000'})
row = 0
col = 0
header_l = ['chrom', 'pos', 'ref', 'alt', 'vcf_string', 'type', 'qual', 'AF_run', 'acc_gene', 'clinvar_clnsig', 'clinvar_gene', 'clinvar_id', 'flankseq', 'gnomad_AF', 'gnomad_flag', 'gnomad_het', 'gnomad_hom', 'hgmd_class', 'hgmd_gene', 'hgmd_id', 'hgmd_phen', 'masked_genes', 'metabo_gene', 'mistic_pred', 'mistic_score', 'omim_inheritance', 'loeuf', 'sysid_disease', 'sysid_group', 'sysid_synopsis', 'vep_consequence', 'vep_impact', 'vep_symbol', 'vep_gene', 'vep_featureType', 'vep_feature', 'vep_biotype', 'vep_exon', 'vep_intron', 'vep_hgvsc', 'vep_hgvsp', 'vep_cdna_pos', 'vep_cds_pos', 'vep_protein_pos', 'vep_aa', 'vep_codons', 'vep_existing_variation', 'vep_distance', 'vep_strand', 'vep_flags', 'vep_variant_class', 'vep_symbol_source', 'vep_hgnc_id', 'vep_canonical', 'vep_mane_select', 'vep_mane_clinical', 'vep_tsl', 'vep_appris', 'vep_ccds', 'vep_hgvs_offset', 'vep_hgvsg', 'vep_clinsig', 'vep_somatic', 'vep_pheno', 'vep_pubmed', 'vep_inframe_oorfs', 'vep_outframe_oorfs', 'vep_uorfs', 'vep_fiveprime_utr_annotation', 'vep_fiveprime_utr_csq', 'vep_spliceai_DP_AG', 'vep_spliceai_DP_AL', 'vep_spliceai_DP_DG', 'vep_spliceai_DP_DL', 'vep_spliceai_DS_AG', 'vep_spliceai_DS_AL', 'vep_spliceai_DS_DG', 'vep_spliceai_DS_DL', 'vep_spliceai_symbol', 'vep_revel', 'vep_cadd_phred', 'vep_cadd_raw', 'sample_GT_ref', 'sample_GT_alt', 'sample_DP', 'sample_AD_ref', 'sample_AD_alt', 'sample_VAF', 'sample_GQ', 'sample_RNC_ref', 'sample_RNC_alt', 'mom_GT_ref', 'mom_GT_alt', 'mom_DP', 'mom_AD_ref', 'mom_AD_alt', 'mom_VAF', 'mom_GQ', 'mom_RNC_ref', 'mom_RNC_alt', 'dad_GT_ref', 'dad_GT_alt', 'dad_DP', 'dad_AD_ref', 'dad_AD_alt', 'dad_VAF', 'dad_GQ', 'dad_RNC_ref', 'dad_RNC_alt']

for item in header_l:
	worksheet.write_string(row, col, item, bold)
	col += 1
row += 1
worksheet.freeze_panes(1, 0)
col = 0

for variant in VCF(vcf_file):
	chrom = variant.CHROM
	pos = variant.POS
	ref = variant.REF
	alt = variant.ALT[0]
	vcf_string = chrom.replace('chr', '') + '-' + str(pos) + '-' + ref + '-' + alt
	v_type = variant.var_type
	v_qual = variant.QUAL
	v_gt = variant.genotype
	#v_AF = variant.INFO['AF']
	v_acc_gene = variant.INFO['ACC_GENE']
	v_clinvar_clnsig = variant.INFO['CLINVAR_CLNSIG']
	v_clinvar_gene = variant.INFO['CLINVAR_GENE']
	v_clinvar_id = variant.INFO['CLINVAR_ID']
	v_flankseq = variant.INFO['FLANKSEQ']
	v_gnomad_AF = variant.INFO['GNOMAD_AF']
	v_gnomad_flag = variant.INFO['GNOMAD_FLAG']
	v_gnomad_het = variant.INFO['GNOMAD_HET']
	v_gnomad_hom = variant.INFO['GNOMAD_HOMHEM']
	v_hgmd_class = variant.INFO['HGMD_CLASS']
	v_hgmd_gene = variant.INFO['HGMD_GENE']
	v_hgmd_id = variant.INFO['HGMD_ID']
	v_hgmd_phen = variant.INFO['HGMD_PHEN']
	v_masked_genes = variant.INFO['MASKED_GENES']
	v_metabo_gene = variant.INFO['METABO_GENE']
	v_mistic_pred = variant.INFO['MISTIC_PRED']
	v_mistic_score = variant.INFO['MISTIC_SCORE']
	v_omim_inheritance = variant.INFO['OMIM_INHERITANCE']
	v_loeuf = variant.INFO['LOEUF']
	v_sysid_disease = variant.INFO['SYSID_DISEASE']
	v_sysid_group = variant.INFO['SYSID_GROUP']
	v_sysid_synopsis = variant.INFO['SYSID_SYNOPSIS']

	''' VEP
	##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|HGVS_OFFSET|HGVSg|CLIN_SIG|SOMATIC|PHENO|PUBMED|existing_InFrame_oORFs|existing_OutOfFrame_oORFs|existing_uORFs|five_prime_UTR_variant_annotation|five_prime_UTR_variant_consequence|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|REVEL|CADD_PHRED|CADD_RAW"> '''

	v_vep = variant.INFO['CSQ']
	v_vep_allele = v_vep.split('|')[0]
	v_vep_consequence = v_vep.split('|')[1]
	v_vep_impact = v_vep.split('|')[2]
	v_vep_symbol = v_vep.split('|')[3]
	v_vep_gene = v_vep.split('|')[4]
	v_vep_featureType = v_vep.split('|')[5]
	v_vep_feature = v_vep.split('|')[6]
	v_vep_biotype = v_vep.split('|')[7]
	v_vep_exon = v_vep.split('|')[8]
	v_vep_intron = v_vep.split('|')[9]
	v_vep_hgvsc = v_vep.split('|')[10]
	v_vep_hgvsp = v_vep.split('|')[11]
	v_vep_cdna_pos = v_vep.split('|')[12]
	v_vep_cds_pos = v_vep.split('|')[13]
	v_vep_protein_pos = v_vep.split('|')[14]
	v_vep_aa = v_vep.split('|')[15]
	v_vep_codons = v_vep.split('|')[16]
	v_vep_existing_variation = v_vep.split('|')[17]
	v_vep_distance = v_vep.split('|')[18]
	v_vep_strand = v_vep.split('|')[19]
	v_vep_flags = v_vep.split('|')[20]
	v_vep_variant_class = v_vep.split('|')[21]
	v_vep_symbol_source = v_vep.split('|')[22]
	v_vep_hgnc_id = v_vep.split('|')[23]
	v_vep_canonical = v_vep.split('|')[24]
	v_vep_mane_select = v_vep.split('|')[25]
	v_vep_mane_clinical = v_vep.split('|')[26]
	v_vep_tsl = v_vep.split('|')[27]
	v_vep_appris = v_vep.split('|')[28]
	v_vep_ccds = v_vep.split('|')[29]
	v_vep_hgvs_offset = v_vep.split('|')[30]
	v_vep_hgvsg = v_vep.split('|')[31]
	v_vep_clinsig = v_vep.split('|')[32]
	v_vep_somatic = v_vep.split('|')[33]
	v_vep_pheno = v_vep.split('|')[34]
	v_vep_pubmed = v_vep.split('|')[35]
	v_vep_inframe_oorfs = v_vep.split('|')[36]
	v_vep_outframe_oorfs = v_vep.split('|')[37]
	v_vep_uorfs = v_vep.split('|')[38]
	v_vep_fiveprime_utr_annotation = v_vep.split('|')[39]
	v_vep_fiveprime_utr_csq = v_vep.split('|')[40]
	v_vep_spliceai_DP_AG = v_vep.split('|')[41]
	v_vep_spliceai_DP_AL = v_vep.split('|')[42]
	v_vep_spliceai_DP_DG = v_vep.split('|')[43]
	v_vep_spliceai_DP_DL = v_vep.split('|')[44]
	v_vep_spliceai_DS_AG = v_vep.split('|')[45]
	v_vep_spliceai_DS_AL = v_vep.split('|')[46]
	v_vep_spliceai_DS_DG = v_vep.split('|')[47]
	v_vep_spliceai_DS_DL = v_vep.split('|')[48]
	v_vep_spliceai_symbol = v_vep.split('|')[49]
	v_vep_revel = v_vep.split('|')[50]
	v_vep_cadd_phred = v_vep.split('|')[51]
	v_vep_cadd_raw = v_vep.split('|')[52]
	sample_GT_ref = variant.genotypes[0][0]
	sample_GT_alt = variant.genotypes[0][1]
	sample_DP = variant.format('DP')[0][0]
	sample_AD_ref = variant.format('AD')[0][0]
	sample_AD_alt = variant.format('AD')[0][0]
	sample_GQ = variant.format('GQ')[0][0]
	#sample_RNC_ref = variant.format('RNC')
	#sample_RNC_alt = variant.format('RNC')
	if sample_AD_alt != -1 and sample_DP != -1 and sample_DP != 0: 
		sample_VAF = float(sample_AD_alt) / float(sample_DP)
	else:
		sample_VAF = None

	#mom_GT_ref = variant.genotypes[mom_idx][0]
	#mom_GT_alt = variant.genotypes[mom_idx][1]
	#mom_DP = variant.format('DP')[mom_idx][0]
	#mom_AD_ref = variant.format('AD')[mom_idx][0]
	#mom_AD_alt = variant.format('AD')[mom_idx][1]
	#mom_GQ = variant.format('GQ')[mom_idx][0]
	#mom_PL = variant.format('PL')[mom_idx]
	#mom_RNC_ref = variant.format('RNC')[mom_idx][0]
	#mom_RNC_alt = variant.format('RNC')[mom_idx][1]
	#if mom_AD_alt != -1 and mom_DP != -1 and mom_DP != 0: 
	#	mom_VAF = float(mom_AD_alt) / float(mom_DP)
	#else:
	#	mom_VAF = None

	#dad_GT_ref = variant.genotypes[dad_idx][0]
	#dad_GT_alt = variant.genotypes[dad_idx][1]
	#dad_DP = variant.format('DP')[dad_idx][0]
	#dad_AD_ref = variant.format('AD')[dad_idx][0]
	#dad_AD_alt = variant.format('AD')[dad_idx][1]
	#dad_GQ = variant.format('GQ')[dad_idx][0]
	#dad_PL = variant.format('PL')[dad_idx]
	#dad_RNC_ref = variant.format('RNC')[dad_idx][0]
	#dad_RNC_alt = variant.format('RNC')[dad_idx][1]
	#if dad_AD_alt != -1 and dad_DP != -1 and dad_DP != 0: 
	#	dad_VAF = float(dad_AD_alt) / float(dad_DP)
	#else:
	#	dad_VAF = None

	worksheet.write_string(row, col, chrom) # chrom (str)
	worksheet.write_number(row, col+1, pos) # pos (int)
	worksheet.write_string(row, col+2, ref) # ref (str)
	worksheet.write_string(row, col+3, alt) # alt (str)
	worksheet.write_string(row, col+4, vcf_string) # vcf_string (str)
	worksheet.write_string(row, col+5, v_type) # v_type (str)
	if v_qual == ".":
		worksheet.write_blank(row, col+6, None) 
	else:
		worksheet.write_number(row, col+6, v_qual) # v_qual (float)
	#worksheet.write_number(row, col+7, v_AF, three_decimals) # v_AF (float)
	if v_acc_gene == 0:
		worksheet.write_blank(row, col+8, None)
	else:
		worksheet.write_number(row, col+8, v_acc_gene) # v_acc_gene (int)
	if v_clinvar_clnsig == ".":
		worksheet.write_blank(row, col+9, None)
	else:
		worksheet.write_string(row, col+9, v_clinvar_clnsig) # v_clinvar_clnsig (str)
	if v_clinvar_gene == ".":
		worksheet.write_blank(row, col+10, None)
	else:
		worksheet.write_string(row, col+10, v_clinvar_gene.replace('"','')) # v_clinvar_gene (str)
	if v_clinvar_id == ".":
		worksheet.write_blank(row, col+11, None)
	else:
		worksheet.write_url(row, col+11, 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' + v_clinvar_id, string = v_clinvar_id)  # v_clinvar_id (str) 
	worksheet.write_string(row, col+12, v_flankseq) # v_flankseq (str)
	if v_gnomad_AF is None:
		worksheet.write_blank(row, col+13, None) # v_gnomad_AF (float)
	else:
		worksheet.write_number(row, col+13, v_gnomad_AF, three_decimals) # v_gnomad_AF (float)
	if v_gnomad_flag == ".":
		worksheet.write_blank(row, col+14, None)
	else:
		worksheet.write_string(row, col+14, v_gnomad_flag) # v_gnomad_flag (str)
	if v_gnomad_het is None:
		worksheet.write_blank(row, col+15, None) # v_gnomad_het (int)
	else:
		worksheet.write_number(row, col+15, v_gnomad_het) # v_gnomad_het (int)
	if v_gnomad_hom is None:
		worksheet.write_blank(row, col+16, None) # v_gnomad_hom (int)
	else:
		worksheet.write_number(row, col+16, v_gnomad_hom) # v_gnomad_hom (int)
	if v_hgmd_class == ".":
		worksheet.write_blank(row, col+17, None)
	else:
		worksheet.write_string(row, col+17, v_hgmd_class) # v_hgmd_class (str)
	if v_hgmd_gene == ".":
		worksheet.write_blank(row, col+18, None)
	else:
		worksheet.write_string(row, col+18, v_hgmd_gene.replace('"','')) # v_hgmd_gene (str)
	if v_hgmd_id == ".":
		worksheet.write_blank(row, col+19, None)
	else:
		worksheet.write_string(row, col+19, v_hgmd_id) # v_hgmd_id (str)
	if v_hgmd_phen == ".":
		worksheet.write_blank(row, col+20, None)
	else:
		worksheet.write_string(row, col+20, v_hgmd_phen) # v_hgmd_phen (str)
	worksheet.write_number(row, col+21, v_masked_genes) # v_masked_genes (int)
	if v_metabo_gene == ".":
		worksheet.write_blank(row, col+22, None)
	else:
		worksheet.write_string(row, col+22, v_metabo_gene) # v_metabo_gene (str)
	if v_mistic_pred == ".":
		worksheet.write_blank(row, col+23, None)
	else:
		worksheet.write_string(row, col+23, v_mistic_pred) # v_mistic_pred (str)
	if v_mistic_score is None:
		worksheet.write_blank(row, col+24, None) # v_mistic_score (float)
	else:
		worksheet.write_number(row, col+24, v_mistic_score) # v_mistic_score (float)
	if v_omim_inheritance == ".":
		worksheet.write_blank(row, col+25, None)
	else:
		worksheet.write_string(row, col+25, v_omim_inheritance) # v_omim_inheritance (str)
	if v_loeuf is None:
		worksheet.write_blank(row, col+26, None) # v_loeuf (float)
	else:
		worksheet.write_number(row, col+26, v_loeuf, three_decimals) # v_loeuf (float)
	if v_sysid_disease == ".":
		worksheet.write_blank(row, col+27, None)
	else:
		worksheet.write_string(row, col+27, v_sysid_disease.replace('"','')) # v_sysid_disease (str)
	if v_sysid_group == ".":
		worksheet.write_blank(row, col+28, None)
	else:
		worksheet.write_string(row, col+28, v_sysid_group) # v_sysid_group (str)
	if v_sysid_synopsis == ".":
		worksheet.write_blank(row, col+29, None)
	else:
		worksheet.write_string(row, col+29, v_sysid_synopsis.replace('"','')) # v_sysid_synopsis (str)
	worksheet.write_string(row, col+30, v_vep_consequence) # v_vep_consequence (str)
	worksheet.write_string(row, col+31, v_vep_impact) # v_vep_impact (str)
	worksheet.write_string(row, col+32, v_vep_symbol) # v_vep_symbol (str)
	worksheet.write_string(row, col+33, v_vep_gene) # v_vep_gene (str)
	worksheet.write_string(row, col+34, v_vep_featureType) # v_vep_featureType (str)
	worksheet.write_string(row, col+35, v_vep_feature) # v_vep_feature (str)
	worksheet.write_string(row, col+36, v_vep_biotype) # v_vep_biotype (str)
	worksheet.write_string(row, col+37, v_vep_exon ) # v_vep_exon (str)
	worksheet.write_string(row, col+38, v_vep_intron) # v_vep_intron (str)
	worksheet.write_string(row, col+39, v_vep_hgvsc) # v_vep_hgvsc (str)
	worksheet.write_string(row, col+40, v_vep_hgvsp) # v_vep_hgvsp (str)
	worksheet.write_string(row, col+41, v_vep_cdna_pos) # v_vep_cdna_pos (str)
	worksheet.write_string(row, col+42, v_vep_cds_pos) # v_vep_cds_pos (str)
	worksheet.write_string(row, col+43, v_vep_protein_pos) # v_vep_protein_pos (str)
	worksheet.write_string(row, col+44, v_vep_aa) # v_vep_aa (str)
	worksheet.write_string(row, col+45, v_vep_codons) # v_vep_codons (str)
	worksheet.write_string(row, col+46, v_vep_existing_variation) # v_vep_existing_variation (str)
	worksheet.write_string(row, col+47, v_vep_distance) # v_vep_distance (str)
	worksheet.write_string(row, col+48, v_vep_strand) # v_vep_strand (str)
	worksheet.write_string(row, col+49, v_vep_flags) # v_vep_flags (str)
	worksheet.write_string(row, col+50, v_vep_variant_class) # v_vep_variant_class (str)
	worksheet.write_string(row, col+51, v_vep_symbol_source) # v_vep_symbol_source (str)
	worksheet.write_string(row, col+52, v_vep_hgnc_id) # v_vep_hgnc_id (str)
	worksheet.write_string(row, col+53, v_vep_canonical) # v_vep_canonical (str)
	worksheet.write_string(row, col+54, v_vep_mane_select) # v_vep_mane_select (str)
	worksheet.write_string(row, col+55, v_vep_mane_clinical) # v_vep_mane_clinical (str)
	worksheet.write_string(row, col+56, v_vep_tsl) # v_vep_tsl (str)
	worksheet.write_string(row, col+57, v_vep_appris) # v_vep_appris (str)
	worksheet.write_string(row, col+58, v_vep_ccds) # v_vep_ccds (str)
	worksheet.write_string(row, col+59, v_vep_hgvs_offset) # v_vep_hgvs_offset (str)
	worksheet.write_string(row, col+60, v_vep_hgvsg) # v_vep_hgvsg (str)
	worksheet.write_string(row, col+61, v_vep_clinsig) # v_vep_clinsig (str)
	worksheet.write_string(row, col+62, v_vep_somatic) # v_vep_somatic (str)
	worksheet.write_string(row, col+63, v_vep_pheno) # v_vep_pheno (str)
	worksheet.write_string(row, col+64, v_vep_pubmed) # v_vep_pubmed (str)
	worksheet.write_string(row, col+65, v_vep_inframe_oorfs) # v_vep_inframe_oorfs (str)
	worksheet.write_string(row, col+66, v_vep_outframe_oorfs) # v_vep_outframe_oorfs (str)
	worksheet.write_string(row, col+67, v_vep_uorfs) # v_vep_uorfs (str)
	worksheet.write_string(row, col+68, v_vep_fiveprime_utr_annotation) # v_vep_fiveprime_utr_annotation (str)
	worksheet.write_string(row, col+69, v_vep_fiveprime_utr_csq) # v_vep_fiveprime_utr_csq (str)
	worksheet.write_string(row, col+70, v_vep_spliceai_DP_AG) # v_vep_spliceai_DP_AG (str)
	worksheet.write_string(row, col+71, v_vep_spliceai_DP_AL) # v_vep_spliceai_DP_AL (str)
	worksheet.write_string(row, col+72, v_vep_spliceai_DP_DG) # v_vep_spliceai_DP_DG (str)
	worksheet.write_string(row, col+73, v_vep_spliceai_DP_DL) # v_vep_spliceai_DP_DL (str)
	worksheet.write_string(row, col+74, v_vep_spliceai_DS_AG) # v_vep_spliceai_DS_AG (str)
	worksheet.write_string(row, col+75, v_vep_spliceai_DS_AL) # v_vep_spliceai_DS_AL (str)
	worksheet.write_string(row, col+76, v_vep_spliceai_DS_DG) # v_vep_spliceai_DS_DG (str)
	worksheet.write_string(row, col+77, v_vep_spliceai_DS_DL) # v_vep_spliceai_DS_DL (str)
	worksheet.write_string(row, col+78, v_vep_spliceai_symbol) # v_vep_spliceai_symbol (str)
	worksheet.write_string(row, col+79, v_vep_revel) # v_vep_revel (str)
	worksheet.write_string(row, col+80, v_vep_cadd_phred, three_decimals) # v_vep_cadd_phred (str)
	worksheet.write_string(row, col+81, v_vep_cadd_raw, three_decimals) # v_vep_cadd_raw (str)
	worksheet.write_number(row, col+82, sample_GT_ref) # sample_GT_ref (int)
	worksheet.write_number(row, col+83, sample_GT_alt) # sample_GT_alt (int)
	worksheet.write_number(row, col+84, sample_DP) # sample_DP (int)
	worksheet.write_number(row, col+85, sample_AD_ref) # sample_AD_ref (int)
	if sample_AD_alt < 0:
		worksheet.write_number(row, col+86, -1)
	else:
		worksheet.write_number(row, col+86, sample_AD_alt) # sample_AD_alt (int)

	if sample_VAF is None:
		worksheet.write_blank(row, col+87, None) # sample_VAF (float)
	else:
		worksheet.write_number(row, col+87, sample_VAF, three_decimals) # sample_VAF (float)

	worksheet.write_number(row, col+88, sample_GQ) # sample_GQ (int)
	#if sample_RNC_ref == ".":
	#	worksheet.write_blank(row, col+89, None)
	#else:
	#	worksheet.write_string(row, col+89, sample_RNC_ref) # sample_RNC_ref (str)
	#if sample_RNC_alt == ".":
	#	worksheet.write_blank(row, col+90, None)
	#else:
	#	worksheet.write_string(row, col+90, sample_RNC_alt) # sample_RNC_alt (str)
	#worksheet.write_number(row, col+91, mom_GT_ref) # mom_GT_ref (int)
	#worksheet.write_number(row, col+92, mom_GT_alt) # mom_GT_alt (int)
	#worksheet.write_number(row, col+93, mom_DP) # mom_DP (int)
	#worksheet.write_number(row, col+94, mom_AD_ref) # mom_AD_ref (int)
	#if mom_AD_alt < 0:
	#	worksheet.write_number(row, col+95, -1)
	#else:
	#	worksheet.write_number(row, col+95, mom_AD_alt) # mom_AD_alt (int)
	#if mom_VAF is None:
	#	worksheet.write_blank(row, col+96, None) # mom_VAF (float)
	#else:
	#	worksheet.write_number(row, col+96, mom_VAF, three_decimals) # mom_VAF (float)
	#worksheet.write_number(row, col+97, mom_GQ) # mom_GQ (int)
	#if mom_RNC_ref == ".":
	#	worksheet.write_blank(row, col+98, None)
	#else:
	#	worksheet.write_string(row, col+98, mom_RNC_ref) # mom_RNC_ref (str)
	#if mom_RNC_alt == ".":
	#	worksheet.write_blank(row, col+99, None)
	#else:
	#	worksheet.write_string(row, col+99, mom_RNC_alt) # mom_RNC_alt (str)
	#worksheet.write_number(row, col+100, dad_GT_ref) # dad_GT_ref (int)
	#worksheet.write_number(row, col+101, dad_GT_alt) # dad_GT_alt (int)
	#worksheet.write_number(row, col+102, dad_DP) # dad_DP (int)
	#worksheet.write_number(row, col+103, dad_AD_ref) # dad_AD_ref (int)
	#if dad_AD_alt < 0:
	#	worksheet.write_number(row, col+104, -1)
	#else:
	#	worksheet.write_number(row, col+104, dad_AD_alt) # dad_AD_alt (int)

	#if dad_VAF is None:
	#	worksheet.write_blank(row, col+105, None) # dad_VAF (float)
	#else:
	#	worksheet.write_number(row, col+105, dad_VAF, three_decimals) # dad_VAF (float)

	#worksheet.write_number(row, col+106, dad_GQ) # dad_GQ (int)
	#if dad_RNC_ref == ".":
	#	worksheet.write_blank(row, col+107, None)
	#else:
	#	worksheet.write_string(row, col+107, dad_RNC_ref) # dad_RNC_ref (str)
	#if dad_RNC_alt == ".":
	#	worksheet.write_blank(row, col+108, None)
	#else:
	#	worksheet.write_string(row, col+108, dad_RNC_alt) # dad_RNC_alt (str)

	row += 1

worksheet.autofilter(0, 0, row-1, 108)
workbook.close()
