#!/usr/bin/env python

# This is code takes input as your interested STR regions and SNP loci and check if any read in the a perticular STR region is enclosing the read and the same time the Snp loci (matching the snp allele with the reference or alternative allele) finally it gives the variation occured in STR region which is associating with the Snp allele.

# packages

import gzip
from tqdm import tqdm
from collections import Counter
import pysam
import sys



def Determine_PofO(NHsnps_regions, Aligned_file_loc):
    
    samfile = pysam.AlignmentFile(f'{Aligned_file_loc}', "rc")

    # list containing determined PofO, (there are still few cases with ambiguity, which will be passed to ambiguity module to clarify them)
    Determined_PofO = []
    
    for repeat_region in tqdm(NHsnps_regions):
        
        # info
        start =  int(repeat_region[1])  
        end =    int(repeat_region[2])
        chrom = repeat_region[0]
        rgn = f'{chrom}-{start}'

        
        snp_pos = int(repeat_region[-5])-1 # comparing with read.get_aligned_pairs(with_seq=True)
        snp_gts = [repeat_region[-4], repeat_region[-3]]
        snp_gts_representation = [repeat_region[-8], repeat_region[-7], repeat_region[-6]]
        denovo_gts = repeat_region[5].split('/') # ***************************************************************change the var name : denovo_gts to offspring alleles
        repeat_len = int(repeat_region[3])

        # here we are calling the function to store the snp allele parent of origin.
        snp_allele_origin = snp_allele_origin_classification(snp_gts, snp_gts_representation)

        # dictionary to maintian the snp_loc and repeat region enclosing the reads.
        read_info = {}

        
        # Here we will itterate the reads which enclose the repeat region and also NHsnp.
        for read in samfile.fetch(chrom, snp_pos, snp_pos):

            #
            if read.mapping_quality < 1: continue

            #
            read_query = read.query_name
            read_start = read.reference_start
            read_end = read.reference_end

            # checking if the ( read is enclosed) and ( covering the interested snp position in the flank). { R1:[s,e],[reference pos],[read base at reference pos] }
            if not read_query in read_info: 
                read_info[read_query] = [] # [STR_alleles, NHsnp]

                #
                base_info = sum(read.get_aligned_pairs(with_seq=True), ())
                snp_on_read_info = snp_check(read, snp_pos, snp_gts, base_info)           
                read_info[read_query].append(snp_on_read_info)


        # # Here we will itterate the reads which enclose the NHsnp.
 
    # Here we are returning list of determined regions.
    return Determined_PofO




# **** Function 1: Snp_AssociatedSTR_GT()
# ****************************************************************************************************************************************************************************************
# Function which returns the Denovo genotype assiciating with the Snp genotype bases which further can be related to Denovo allele parent of origin.

def Snp_AssociatedSTR_GT(read_info, snp_gts, denovo_gts, start, end, snp_pos):
    
    # A copy of read_info dictionary to delete the read id's which are not usefull. ()
    read_info_copy = read_info.copy()
    Offspring_A12 = {denovo_gts[0]:[], denovo_gts[1]:[]}
    
    # Before entering into the current repeat region read_info dictionary, we must ignore if ( read is not enclosed the repeat region | No snp on any read | or snp on R1 or R2 but no read overlapping the repeat region.
    temp = [Offspring_A12[str(read_info_copy[read_id][1])].append(read_info_copy[read_id][0]) for read_id in read_info_copy  if len(read_info_copy[read_id]) == 2 and str(read_info_copy[read_id][1]) in denovo_gts]
    
    A12 = {denovo_gts[0]:[], denovo_gts[1]:[]}    

        # 
    return NHsnp_asso_offspring_alleles
            



# Function2 : snp_check()
# ****************************************************************************************************************************************************************************************
# Creating function to check the snp_pos overlapping with which read and matching with which snp_allele

def snp_check(read, snp_pos, snp_gts, base_info):
    
    # base info = [read_pos, refn_pos, read_base]
    reference_pos = list(base_info)[1::3]
    read_bases = list(base_info)[2::3]
    
    # check if read is overalpped with snp_pos and base.



# Function 3 : STR_variation_from_read()
# ****************************************************************************************************************************************************************************************
# Considering only Repeat overlapped CIGAR variants.

def STR_variation_from_read(read, start, end, repeat_len):

        
    # Finding out the way to get the CIGAR string pos in read, reference & type of variation.
    # Then figuring out the way to include the way to fecth the info in our set of reads.

    repeat_start  = start
    repeat_end    = end


    cigar_s = read.cigarstring
    cigar_t = sum(read.cigartuples, ())

    num_bases = list(cigar_t)[1::2]
    alignment_type = list(cigar_t)[::2]

    # We update while itterating through the read positions to cover the entire read positions with any len and capture the variations which encounter within the repeat region.
    return len_diff_in_read





# Function 4 : snp_allele_origin_classification()
# ****************************************************************************************************************************************************************************************
# Here we relate snp alleles to the maternal and paternal side.

def snp_allele_origin_classification(snp_gts, snp_gts_representation):
    
    Snp_origin = []
    
    # We have to 0 or 1 being unique in list of GT combinations.
    offspring = snp_gts_representation[0].split('|')
    parents = snp_gts_representation[1].split('|') + snp_gts_representation[2].split('|')
    GT_bases = {'0':snp_gts[0], '1':snp_gts[1]}
    
    # Based on the '0' and '1' counts, either being count 1 where that allele indicates it is coming from criteria satisfied condition.
    p1 =  parents.count(offspring[0])
    p2 = parents.count(offspring[1])
    
    #
    paternal = snp_gts_representation[1]#.split('|')
    maternal = snp_gts_representation[2]#.split('|')

 
    return Snp_origin


