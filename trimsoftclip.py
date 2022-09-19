#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 15:54:44 2022

@author: mobilab5
"""

import os
import pandas as pd
import pysam
import argparse
import string
# FUNCTION
def make_subcigar(cigar):    
    """ 
    Break cigar string in .sam file into parts.
    
    Args:
    ----------
    cigar_S : STR 
        cigar string in .sam file
        
    Returns:
    ----------
    subcigar_L : LIST
        list of subcigar
        
    Example:
    ----------
    >>> subcigar = make_subcigar('150M2D15I27M')
    >>> subcigar
    ['150M', '2D', '15I', '27M']
    """
    
    alphabet_list = list(string.ascii_uppercase)
    subcigar_L = []
    cigar_len = len(cigar)
    i = 0
    while i < cigar_len:
        sub_cigar = ""
        while cigar[i] not in alphabet_list:
            sub_cigar = sub_cigar + cigar[i]
            i += 1
        sub_cigar = sub_cigar + cigar[i]
        subcigar_L.append(sub_cigar)
        i += 1
    return subcigar_L


# PARSE ARGUMENTS
__author__      = "TA-Nguyen"
__version__     = "1.0.0"
__copyright__   = "Copyright 2022, MolBioLab"

parser = argparse.ArgumentParser(description = 'Trim soft-clipped part from aligned read of SAM file format',
					usage = 'Use "python %(prog)s --help" for more information \n Simple usage: python3 trimsoftclip.py -i <in.sam> -o <out.sam>',
        			formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i','--input-file',  nargs='+', dest="inFile", action = 'store',
                    help='a SAM file contain soft-clip part in alignment report')

parser.add_argument('-o','--output-file',  nargs='+', dest="outFile", action = 'store',
                    help='a out SAM file with trimmed soft-clip part in alignment report')

args = parser.parse_args()

# RUN
if len(args.inFile) != 1:
    exit()
    
in_samfile = pysam.AlignmentFile(args.inFile[0], "rb", check_sq = False)

# CREATE NEW SAM file
with open(args.outFile[0], "+a") as out:
    out.write(in_samfile.text)
    
# TRIM SOFT-CLIP
for read in in_samfile.fetch():
    cigar_T = read.cigar
    
    # A. Create new CIGAR field
    original_cigar = read.to_dict()["cigar"]
    original_subcigar_L = make_subcigar(original_cigar)
    new_subcigar_L = [subCigar for subCigar in original_subcigar_L if "S" not in subCigar]
    new_cigar = "".join(new_subcigar_L)    
    
    
    # B. Create new POS field
    new_Startpos = read.pos
    
        # Check first soft-clipped
    cigarID_first, cigarLen_first = cigar_T[0]
    if cigarID_first == 4: 
        new_Startpos = new_Startpos + cigarLen_first
        
        
    # C. Create new SEQ & QUAL field
    new_seq = read.seq
    new_qual = read.qual
    
        # 1. Check first soft-clipped
    cigarID_first, cigarLen_first = cigar_T[0]
    if cigarID_first == 4: 
        new_seq = new_seq[cigarLen_first:]
        new_qual = new_qual[cigarLen_first:]
        
        # 2. Check last soft-clipped
    cigarID_last, cigarLen_last = cigar_T[-1]
    if cigarID_last == 4: 
        new_seq = new_seq[:-cigarLen_last]
        new_qual = new_qual[:-cigarLen_last]
    
    # D. Create new field: LW:i:<Length of sequence Without soft-clipped part>    
    new_seqLen = len(read.seq)
    for cigarID, cigarLen in cigar_T:
        if cigarID == 4:
            new_seqLen = new_seqLen - int(cigarLen)
    
    
    # E. Write new SAM file
    with open(args.outFile[0], "+a") as out:
        aligned_info = read.to_dict()
        fields = tuple(aligned_info.keys())
        
        # Change some fields info
        aligned_info["ref_pos"] = new_Startpos
        aligned_info["cigar"] = new_cigar
        aligned_info["seq"] = new_seq
        aligned_info["qual"] = new_qual
        aligned_info["tags"].append("LW:i:" + str(new_seqLen))
        aligned_info["tags"] = "\t".join(aligned_info["tags"])
        
        
        newLine_L = [str(aligned_info[field]) for field in fields]
        newLine = '\t'.join(newLine_L)   
        out.write(newLine + "\n")
    