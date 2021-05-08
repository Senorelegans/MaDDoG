#!/usr/bin/env python
import argparse
import csv
import numpy as np
import os
import shlex
import shutil
import subprocess
import sys
import pandas as pd
import string
import glob

# Running in parallel
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Part II: Calculating run in/on for sense and antisense strands.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--annotation_file_with_path', action= 'store', metavar='annotation_file_with_path', help='ensembl or wormbase annotation file for your organism.  The file may be downloaded from http://www.ensembl.org/info/data/ftp/index.html for ensembl or ftp://ftp.wormbase.org/pub/wormbase/releases/ for wormbase')
    parser.add_argument('--control_BAM_list', nargs="*", action= 'store', metavar='input_path_BAM_folder', help='List of files for your control samples')
    parser.add_argument('--treatment_BAM_list', nargs="*", action= 'store', metavar='input_path_BAM_folder', help='List of files for treated samples')
    parser.add_argument('--input_R_template_file', action= 'store', metavar='input_R_template_file', help='input_R_template_file. Example "R_template/R_template_new.R" ')    
    parser.add_argument('--cpus', action= 'store', dest='cpus', metavar='cpus', default= 1, type=int, help='Enter available cpus per node.  The more cpus the faster Dogcatcher performs. Default: "1"')
    parser.add_argument('--padj', action= 'store', dest='padj', metavar='padj', default= 0.05, type=float, help='DSeq2 padj value cutoff to select non-significant genes. Default: "0.05"')

    args = parser.parse_args()

    annotation_file_with_path = args.annotation_file_with_path

    gtf = annotation_file_with_path

    annotation_file = os.path.basename(annotation_file_with_path)
    annotation_path = os.path.dirname(annotation_file_with_path)
    control_BAM_list = args.control_BAM_list
    treatment_BAM_list = args.treatment_BAM_list

    input_R_template_file = args.input_R_template_file
    padj = args.padj
    cpus = args.cpus



    bamlist = str(control_BAM_list + treatment_BAM_list).strip("[").strip("]")

    ###Generate col data file
    control_BAM_list = [  "X"+bam.replace("/",".") + "\tC" for bam in control_BAM_list ]
    treatment_BAM_list = [  "X"+bam.replace("/",".") + "\tT" for bam in treatment_BAM_list ]
    control_BAM_list = [ bam.replace("..",".") for bam in control_BAM_list ]
    treatment_BAM_list = [ bam.replace("..",".") for bam in treatment_BAM_list ]
    with open("col_data.txt", "w") as outfile:
        outfile.write("\tcondition" + "\n")
        for bam in control_BAM_list:
            outfile.write(bam + "\n")
        for bam in treatment_BAM_list:
            outfile.write(bam + "\n")


    def create_R_script():
        with open(input_R_template_file, "r") as infile, open("Rsubread_DESeq2_initial.R", "w") as outfile:
            for line in infile:
                line = line.replace("COL_DATA", "col_data.txt")
                line = line.replace("bamlist", bamlist)
                line = line.replace("initial_annotation_file", str(gtf))
                line = line.replace("cpus", str(cpus))
                outfile.write(line)
    create_R_script()
