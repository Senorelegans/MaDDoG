
#################################################
### THIS FILE WAS AUTOGENERATED! DO NOT EDIT! ###
#################################################
# file to edit: dev_nb/20_Dogcatcher.ipynb
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
import multiprocessing as mp
import time
from nb_00_HelperFunctions import *

import numpy as np
import scipy.stats as stats



def makeFolders(folder_list):
    for directory in folder_list:
        if not directory: # Make sure not an empty string ""
            continue
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)

def removeFile(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

def removeFolder(foldername):
    try:
        shutil.rmtree(foldername)
    except OSError:
        print(f"{foldername} not found")
        pass

def match_bed_gtf(df_gtf, df_bed, printSTUFF=False):
    #print("Ensuring chr names in bed match annotation")
    gtf_names=df_gtf['chr'].unique().tolist()

    gtf_only = df_gtf["chr"].unique().tolist()
    bed_only = df_bed["chr"].unique().tolist()

    df_gtf = df_gtf[df_gtf["chr"].isin(df_bed["chr"])]
    df_bed = df_bed[df_bed["chr"].isin(df_gtf["chr"])]

    chr_names = df_gtf['chr'].unique().tolist()

    if printSTUFF:
    #     if (len(df_gtf) == 0) & (len(df_bed)==0):
    #         print("Amount/Names of Chromosomes in annotation or bed file don't match ")
    #         print("annotation chr names detected \n ", gtf_only)
    #         print("bed chr names detected \n ", bed_only)

        print("Detected chr names in annotation and bed file are ", chr_names)
    return df_gtf, df_bed, chr_names

def get_bg_df(f1, df_gtf):
    my_col = ["chr", "start", "end","count"]
    df_bed = pd.read_csv(f1, sep="\t", header=None, names=my_col, comment="#", dtype={"chr": str, "start": int, "end": int, "count": float}, low_memory=False)
    df_bed = df_bed[~df_bed["chr"].str.contains("\.") ]    # Take out patches

    df_bed["length"] = df_bed["end"] - df_bed["start"]
    df_gtf, df_bed, chr_names = match_bed_gtf(df_gtf, df_bed)
    df_bed["count"] = df_bed["count"].abs()
    df_b = df_bed[ df_bed["count"] > 0.0 ]
    return df_bed, df_b, df_gtf, chr_names

def writeBGOverZero(f1):
    my_col = ["chr", "start", "end","count"]
    df = pd.read_csv(f1, sep="\t", header=None, names=my_col, comment="#", dtype={"chr": str, "start": int, "end": int, "count": float}, low_memory=False)
    df = df[~df["chr"].str.contains("\.") ]    # Take out patches
    df = df[ df["count"] > 0.0 ]
    df.to_csv(f1[:-4]+".BedGraph",sep="\t",header=None,index=None)


def strandDF(df):
    df_plu = df[df["strand"]=="+"]
    df_min = df[df["strand"]=="-"]
    return df_plu, df_min

def get_chrm_dic(df):
    df["length"] = df["end"] - df["start"]
    names=df['chr'].unique().tolist()
    d = {chrom : df.loc[df.chr==chrom] for chrom in names}
    return d, names

def get_cov_perc(df, chr, start, end):
    """this function will take start and end and return the percent coverage over the area.
    for a BedGraph file for only one chromosome"""
    cov_per = 0.0
    rows_df = df[ (df["end"] > start) & (df["start"] < end) ]
    if len(rows_df) == 0:
        return cov_per
    else:
        start_no_cov_length = 0
        if start > rows_df["start"].iloc[0]:
            start_no_cov_length = start - rows_df["start"].iloc[0]

        end_no_cov_length = 0
        if end < rows_df["end"].iloc[-1]:
            end_no_cov_length = rows_df["end"].iloc[-1] - end

        cov_sum = rows_df["length"].sum() - (start_no_cov_length + end_no_cov_length)    #add new start and end distance and all rows inbetween
        cov_length = end - start
        if cov_length == 0:
            cov_per = 0.0
            return cov_per
        else:
            cov_per = cov_sum / cov_length
            return cov_per
    return 0.0

def get_right_DOG(df, df_bed_chr, w, cp, plu_min):
    l1 = len(df)
    if plu_min == "plu":
        D = ""
    if plu_min == "min":
        D = "A"

    df["DOG_length"] = 0
    df["cov_perc"] = 0.0
    df.sort_values(by=["start","chr"], inplace=True, ascending=True)
    df_right_overlap = df[ df["start"].shift(-1) < df["end"] ]

    df["right_int_length"] = df["start"].shift(-1) - df["end"]
    df["right_int_length"] = df[['right_int_length']].fillna(value=100_000)

    cols = list(df)

    df_final = pd.DataFrame(columns=cols)
    df2 = df [ df['right_int_length'] <= w]
    if len(df2) > 0:
        df2['cov_perc'] = df2.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + row["right_int_length"]), axis=1)
        df2["DOG_length"] = np.where(df2['cov_perc'] > cp, df2['right_int_length'], 0 )
        df_final = df_final.append(df2)


    df = df [ df['right_int_length'] > w]

    c = w
    while len(df) > 0:
        print("df: ", len(df))
        df['cov_perc'] = df.apply(lambda row: get_cov_perc(df_bed_chr, row['chr'], row['end'], row['end'] + w), axis=1)

        dfu = df [ df['cov_perc'] < cp ]
        df = df [ df['cov_perc'] >= cp]
        if len(dfu) > 0:
            dfu["DOG_length"] = dfu["DOG_length"] + w - c
            df_final = df_final.append(dfu)

        w = w + c
    df_final.sort_values(by=["start","chr"], inplace=True, ascending=True)
    l2 = len(df_final)

    return df_final

def CalcCoverage(df, start, end):
    #(df_bed_chr, row['chr'], row['end'], row['end'] + row["right_int_length"]), axis=1)

    start = int(start)
    end = int(end)
    cov = 0
    rows_df = df[ (start < df["end"] ) & (df["start"] < end) ]
    if len(rows_df) == 0:
        return cov
    else:
        if rows_df["start"].iloc[0] < start:
            rows_df["start"].iloc[0] = start

        if end < rows_df["end"].iloc[-1]:
            rows_df["end"].iloc[-1] = end

        rows_df["length"] = rows_df["end"] - rows_df["start"]
        rows_df["count_by_length"] = rows_df["count"] * rows_df["length"]
        cov_sum = rows_df["count_by_length"].sum()
        cov_length = end - start
        if cov_length == 0:
            cov = 0
        else:
            cov = cov_sum / cov_length
    return cov

def ReturnReads(df, start, end):

    rows_df = df[ (start < df["end"] ) & (df["start"] < end) ]
    rows_df = rows_df[["start","end","count"]]

    if len(rows_df) > 0:
        rows_df["start"].iloc[0] = start
        rows_df["end"].iloc[-1] = end
        d = rows_df.to_dict("list")
        return d
    else:
        return "NA"

def SORT(df):
    df.sort_values(by=["chr","start"], inplace=True, ascending=True)
    return df

def setInt(df):
    df["DOG_start_local"] = df["DOG_start_local"].astype(int)
    df["DOG_start_meta"] = df["DOG_start_local"].astype(int)
    df["DOG_end_local"] = df["DOG_end_local"].astype(int)
    df["DOG_end_meta"] = df["DOG_end_meta"].astype(int)
    df["DOG_length"] = df["DOG_length"].astype(int)
    return df

def Rcoverage(df, df_bed_chr, w, dog_gene_perc_cov, region, region_coverage_min,InsideGeneSize):
    dog_gene_perc_cov = dog_gene_perc_cov / 100
    df["DOG_length"] = 0
    df = SORT(df)

    df["right_int_length"] = df["start"].shift(-1) - df["end"]
    df_right_overlap = df[ df["start"].shift(-1) < df["end"] ]
    df = df[~df["name"].isin(df_right_overlap["name"])]

    df["right_int_length"] = df[['right_int_length']].fillna(value=100_000)


    if region == "LAST_EXON":
        df["region_start"] = df["LA_start"]
        df["region_end"] = df["LA_end"]
    if region == "GENE":
        df["region_start"] = df["start"]
        df["region_end"] = df["end"]
    if region == "WINDOW":
        df["region_start"] = df["end"] - w
        df["region_start"] = np.where(df["region_start"] < df["start"], df["start"], df["region_start"] )
        df["region_end"] = df["end"]

    df["region_cov"] = df.apply(lambda row: CalcCoverage(df_bed_chr, row["region_start"],row["region_end"]), axis=1)
    dfr = df[df["region_cov"]<=region_coverage_min]
    df = df[df["region_cov"]>region_coverage_min]


    df_final = pd.DataFrame(columns=list(df))
    df2 = df [ df['right_int_length'] <= w]
    if len(df2) > 0:
        df2['DOG_cov'] = df2.apply(lambda row: CalcCoverage(df_bed_chr, row['end'], row['end'] + row["right_int_length"]), axis=1)
        df2["DOG_region_perc_cov"] =  df2["DOG_cov"] / df2["region_cov"]
        df2["DOG_length"] = np.where(df2["DOG_region_perc_cov"] > dog_gene_perc_cov, df2['right_int_length'], 0 )
        df_final = df_final.append(df2)

    df = df [ df['right_int_length'] > w]


    c1 = 0
    c2 = w
    while len(df) > 0:

        df['DOG_cov'] = df.apply(lambda row: CalcCoverage(df_bed_chr, row['end']+ c1, row['end']+ c2), axis=1)
        df["DOG_region_perc_cov"] =  df["DOG_cov"] / df["region_cov"]

        dfu = df [ df['DOG_region_perc_cov'] < dog_gene_perc_cov ]
        df = df [ df['DOG_region_perc_cov'] >= dog_gene_perc_cov]
        if len(dfu) > 0:
            dfu["DOG_length"] = c1
            df_final = df_final.append(dfu)
        c1 = c1 + w
        c2 = c2 + w

    df = pd.concat([df_final,df_right_overlap,dfr])
    df = SORT(df)
    df["TYPE"] = "DOG"
    df["DOG_start_local"] = df["end"]
    df["DOG_start_meta"] = df["end"]
    df["DOG_end_meta"] = df["DOG_start_local"] + df["DOG_length"]
    df["DOG_into_downstream_gene"] = ( df["DOG_length"] >= df['right_int_length'] ) & ( df["DOG_length"] > 0 )
    df["DOG_end_local"] = np.where(df["DOG_into_downstream_gene"] == True , df["end"]+df["right_int_length"], df["DOG_end_meta"] )
    df["DOG_from_upstream_gene"] = np.where(df["DOG_into_downstream_gene"].shift(1) == True , True, False)
    df = df[df["DOG_length"]>0]
    df = setInt(df)
    #ReturnReads(df, start, end)
    if len(df)> 0:
        df['DOG_reads_local'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['DOG_start_local'], row["DOG_end_local"]), axis=1)
        df['DOG_reads_meta'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['DOG_start_meta'], row["DOG_end_meta"]), axis=1)
        df['ALL_reads_local'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['region_start'], row["DOG_end_local"]), axis=1)
        df['ALL_reads_meta'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['region_start'], row["DOG_end_meta"]), axis=1)
    return df

def Lcoverage(df, df_bed_chr, w, dog_gene_perc_cov, region, region_coverage_min,InsideGeneSize):

    dog_gene_perc_cov = dog_gene_perc_cov / 100
    df["DOG_length"] = 0
    df = SORT(df)

    df_left_overlap = df[df["start"] < df["end"].shift(1)]
    df['left_int_length'] = (df["end"].shift(1) - df["start"]) * (-1)
    df["left_int_length"] = df[['left_int_length']].fillna(value=100_000)

    df = df[~df["name"].isin(df_left_overlap["name"])]

    if region == "LAST_EXON":
        df["region_start"] = df["FA_start"]
        df["region_end"] = df["FA_end"]
    if region == "GENE":
        df["region_start"] = df["start"]
        df["region_end"] = df["end"]
    if region == "WINDOW":
        df["region_end"] = df["start"] + InsideGeneSize
        df["region_end"] = np.where(df["region_end"] > df["end"], df["end"], df["region_end"] )
        df["region_start"] = df["start"]

    df["region_length"] = df["region_end"] - df["region_start"]
    df["region_cov"] = df.apply(lambda row: CalcCoverage(df_bed_chr, row["region_start"],row["region_end"]), axis=1)
    dfr = df[df["region_cov"]<=region_coverage_min]
    df = df[df["region_cov"]>region_coverage_min]



    df_final = pd.DataFrame(columns=list(df))
    df2 = df [ df['left_int_length'] <= w]
    if len(df2) > 0:
        df2['DOG_cov'] = df2.apply(lambda row: CalcCoverage(df_bed_chr, row['start'] - row["left_int_length"], row["start"]), axis=1)
        df2["DOG_region_perc_cov"] =  df2["DOG_cov"] / df2["region_cov"]
        df2["DOG_length"] = np.where(df2["DOG_region_perc_cov"] > dog_gene_perc_cov, df2['left_int_length'], 0 )
        df_final = df_final.append(df2)

    df = df [ df['left_int_length'] > w]

    c1 = 0
    c2 = w
    while len(df) > 0:
        df['DOG_cov'] = df.apply(lambda row: CalcCoverage(df_bed_chr, row['start'] - c2, row['start'] - c1), axis=1)
        df["DOG_region_perc_cov"] =  df["DOG_cov"] / df["region_cov"]

        dfu = df [ df['DOG_region_perc_cov'] < dog_gene_perc_cov ]
        df = df [ df['DOG_region_perc_cov'] >= dog_gene_perc_cov]
        if len(dfu) > 0:
            dfu["DOG_length"] = c1
            df_final = df_final.append(dfu)
        c1 = c1 + w
        c2 = c2 + w

    df = pd.concat([df_final,df_left_overlap,dfr])
    df = SORT(df)
    df["TYPE"] = "DOG"
    df["DOG_end_local"] = df["start"]
    df["DOG_end_meta"] = df["start"]
    df["DOG_start_meta"] = df["start"] - df['DOG_length']
    df["DOG_into_downstream_gene"] = (df['DOG_length'] >= df['left_int_length']) & (df['DOG_length'] > 0)
    df["DOG_start_local"] = np.where(df["DOG_into_downstream_gene"] == True, df["start"]-df['left_int_length'],df["DOG_start_meta"])
    df["DOG_from_upstream_gene"] = np.where(df["DOG_into_downstream_gene"].shift(1) == True, True,False)
    df = df[df["DOG_length"]>0]
    df = setInt(df)
    if len(df)> 0:
        df['DOG_reads_local'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['DOG_start_local'], row["DOG_end_local"]), axis=1)
        df['DOG_reads_meta'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['DOG_start_meta'], row["DOG_end_meta"]), axis=1)
        df['ALL_reads_local'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['DOG_start_local'], row["region_end"]), axis=1)
        df['ALL_reads_meta'] = df.apply(lambda row: ReturnReads(df_bed_chr, row['DOG_start_meta'], row["region_end"]), axis=1)

    return df

class DOG():
    def __init__(self, OutFolder, Title, fname,Annotation, f1_plu, f1_min,cpus=4, RunDOG=True,SlidingWindow=1000, DogGenePercCov=1 ,GeneRegion='WINDOW', GeneRegionCoverageMin=3, InsideGeneSize=1000, CoveragePerc=5,splitBed=False):
        pd.set_option('mode.chained_assignment',None)
        self.fname = fname
        self.OutFolder=OutFolder
        self.DogGenePercCov=DogGenePercCov
        self.Region = GeneRegion
        self.RunDOG = RunDOG
        self.RegionCoverageMin = GeneRegionCoverageMin
        self.Title = Title
        self.InsideGeneSize = InsideGeneSize

        self.splitBed = splitBed
        self.cpus = cpus
        self.SlidingWindow = SlidingWindow
        self.CoveragePerc = CoveragePerc / 100
        gtf = f"{Annotation[:-4]}_flat.txt"
        self.f1_plu = f1_plu
        self.f_plu = os.path.basename(f1_plu)
        self.f1_min = f1_min
        self.f_min = os.path.basename(f1_min)

        t = fExist(self.f1_plu)
        t = fExist(self.f1_min)

        f1s = self.f1_min.split("_min")
        self.f1 = f1s[0]
        self.ext = f1s[1]
        self.f = os.path.basename(self.f1)

        makeFolders([OutFolder,f"{OutFolder}/beds",f"{OutFolder}/csv"])

        self.plu_min_gtf = "plu"


        self.df_gtf = pd.read_csv(gtf,sep="\t",dtype={"chr": str, "start": int,
                                  "end": int,"name":str,"strand":str,
                                  "FA_start":str,"FA_end":str,"LA_start":str,"LA_end":str},low_memory=False)
        self.df_gtf_plu, self.df_gtf_min = strandDF(self.df_gtf)

        self.chr_names = self.df_gtf['chr'].unique().tolist()
        self.chr_names_plu = self.df_gtf_plu['chr'].unique().tolist()
        self.chr_names_min = self.df_gtf_min['chr'].unique().tolist()


        self.d_gtf_plu, self.chr_names_plu = get_chrm_dic(self.df_gtf_plu)
        self.d_gtf_min, self.chr_names_min = get_chrm_dic(self.df_gtf_min)

        self.makeBeds()
        self.multiCOV()
        self.parseIGV()

    def setStrand(self):
        if self.plu_min == "plu":
            self.strand = "+"
        if self.plu_min == "min":
            self.strand = "-"
        if self.plu_min_gtf == "plu":
            self.strand_gtf = "+"
        if self.plu_min_gtf == "min":
            self.strand_gtf = "-"


    def tempBed(self, chrom):
            fout = f"{self.fol}/{self.f}_{str(chrom)}_{self.plu_min}.BedGraph"
            print(f"chrom: {chrom} worker: {os.getpid()}")
            with open(self.fin, "r", newline="\n") as infile, open(fout, "w", newline="\n") as outfile:
                for line in infile:
                    tab1 = line.index("\t")
                    line_chr = line[:tab1]
                    if line_chr == chrom:
                        outfile.write(line)

    def makeBeds(self):
        self.d = "Dogcatcher_temporary_BedGraph"
        self.fol = f"{self.d}/{self.f}"
        makeFolders([self.d])

        if fExist(self.fol):
                    print(f"Not splitting... because BedGraphs detected: {self.fol} or splitBed is False")
        if (fExist(self.fol) == False) or (self.splitBed == True):
            print(f"Splitting BedGraphs to {self.fol}")
            makeFolders([self.fol])

            for self.fin in [self.f1_plu, self.f1_min]:
                if self.fin == self.f1_plu:
                    self.plu_min = "plu"
                if self.fin == self.f1_min:
                    self.plu_min = "min"
                #Start multiprocess with temp bedgraph
                pool    = mp.Pool(processes=self.cpus)
                multiple_results = [pool.apply_async(self.tempBed, args=(chrom,)) for chrom in self.chr_names]
                pool.close()
                pool.join()  # block at this line until all processes are done
                print(f"Finished Multiprocess: {self.fol}")

    def COV(self, chrom):
        print(f"chrom: {chrom} worker: {os.getpid()}")
        fout = f"{self.fol}/{self.f}_{str(chrom)}_{self.plu_min}.BedGraph"

        df_gtf = self.df_gtf[self.df_gtf["strand"]==self.strand]
        df_bed, df_b, df_gtf, chr_names = get_bg_df(fout, df_gtf)
        if self.RunDOG:
            if self.plu_min == "plu":
                df = Rcoverage(df_gtf, df_bed, self.SlidingWindow, self.DogGenePercCov, self.Region, self.RegionCoverageMin,self.InsideGeneSize)
            if self.plu_min == "min":
                df = Lcoverage(df_gtf, df_bed, self.SlidingWindow, self.DogGenePercCov, self.Region, self.RegionCoverageMin,self.InsideGeneSize)

        if self.RunDOG != True:
                df = df_gtf

                df['reads_local'] = df.apply(lambda row: ReturnReads(df_bed, row['start'], row["end"]), axis=1)
                df["region_start"] = df["start"]
                df["region_end"] = df["end"]

        if self.RunDOG:
                df["DOG_length"] = df["DOG_end_local"] - df["DOG_start_local"] - self.InsideGeneSize

        print(self.f)
        df["sample_name"] = self.f
        df["region"] = self.Region

        df["DogGenePercCov"] = self.DogGenePercCov
        df["RegionCoverageMin"] = self.RegionCoverageMin
        df["SlidingWindow"] = self.SlidingWindow
        df["InsideGeneSize"] = self.InsideGeneSize
        return df



    def multiCOV(self):
        df = pd.DataFrame()
        for self.plu_min in ["plu","min"]:
            self.setStrand()
            print(f"Processing {self.plu_min} strand DOG", "*-*"*10)
            pool = mp.Pool(processes=self.cpus)
            # Dont forget to take out later
#             self.chr_names = ["chr4"]
            multiple_results = [pool.apply_async(self.COV, args=(chrom,)) for chrom in self.chr_names]
            pool.close()
            pool.join()
            print("All multiprocess workers done... ")
            dfs = ([ res.get() for res in multiple_results])

            for df2 in dfs:
                df = pd.concat([df,df2])

        df.to_csv(f"{self.OutFolder}/csv/{self.Title}_{self.f}.csv",sep="\t",index=None)
        self.df = df
        df2 = df.copy()
        del df2["DOG_reads_local"]
        del df2["DOG_reads_meta"]
        del df2["ALL_reads_local"]
        del df2["ALL_reads_meta"]
        df2.to_csv(f"{self.OutFolder}/csv/{self.Title}_{self.f}_NoReads.csv",sep="\t",index=None)


    def parseIGV(self):
        purple = "157,50,238"
        light_green = "144,238,144"
        df = self.df
        df["RGB"] = np.where(df["strand"]=="+",light_green,purple)
        df["start2"] = df["start"]
        df["stop2"] = df["end"]
        df["zero"] = "0"

        df = df[["chr","DOG_start_local","DOG_end_local","name","zero","strand","start2","stop2","RGB"]]
        fout = f"{self.OutFolder}/beds/{self.Title}_{self.f}.bed"
        open(fout, 'w').write('track name="'+self.f+'" description="'+self.f+'" visibility=2 itemRgb="On"\n')
        df.to_csv(fout, sep="\t", index=None, columns=None, header=None)





import argparse
def parse_arguments():
        parser = argparse.ArgumentParser(description='Flatten gtf or bed to first and last exon file. Options in currently are ENSEMBL, BED')
        parser.add_argument('--annotation_in', action= 'store', metavar='annotation_in')
        parser.add_argument('--f1_plu', action= 'store', metavar='bg_plu')
        parser.add_argument('--f1_min', action= 'store', metavar='bg_min')
        parser.add_argument('--cpus', action= 'store', metavar='cpus',default=4)
        parser.add_argument('--GeneRegion', action= 'store', metavar='GeneRegion')
        parser.add_argument('--RunDOG', action= 'store', metavar='RunDOG')
        parser.add_argument('--DogGenePercCov', action= 'store', metavar='DogGenePercCov')
        parser.add_argument('--RegionCoverageMin', action= 'store', metavar='RegionCoverageMin')
        parser.add_argument('--InsideGeneSize', action= 'store', metavar='InsideGeneSize')
        parser.add_argument('--SlidingWindow', action= 'store', metavar='SlidingWindow')
        parser.add_argument('--OutFolder', action= 'store', metavar='OutFolder')


        args = parser.parse_args()
        return args

if __name__=="__main__":
    args = parse_arguments()
    Annotation = args.annotation_in
    f1_plu = args.f1_plu
    f1_min = args.f1_min
    cpus = int(args.cpus)
    GeneRegion = args.GeneRegion
    RunDOG = args.RunDOG
    DogGenePercCov = int(args.DogGenePercCov)
    RegionCoverageMin = args.RegionCoverageMin
    InsideGeneSize = args.InsideGeneSize
    SlidingWindow = args.SlidingWindow
    Title = f"Region_{Region}_SlidingWindow_{SlidingWindow}_DogGenePercCov_{DogGenePercCov}"
    OutFolder=f"{args.OutFolder}_{Title}"
    D = DOG(OutFolder=OutFolder, Title=Title, Annotation=Annotation, f1_plu=f1_plu,
            f1_min=f1_min,cpus=4, RunDOG=True,SlidingWindow=SlidingWindow, DogGenePercCov=DogGenePercCov,
            GeneRegion=GeneRegion, GeneRegionCoverageMin=GeneRegionCoverageMin, InsideGeneSize=InsideGeneSize)
