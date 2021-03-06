
import sys
import os

def main(Bedgraphfile, flagstatfile, readcountcorrBedgraphfile):
    total_reads = open(flagstatfile).readlines()[0]
    with open(Bedgraphfile) as infile, open(readcountcorrBedgraphfile, "w") as outfile:
        for line in infile:
            line = line.strip("\n")
            line = line.split("\t")
            num_of_reads = line[-1]
            mp = str(float(num_of_reads) / (int(total_reads)/1_000_000))
            line[-1] = mp
            outfile.write("\t".join(map(str,line))+"\n")


if __name__=="__main__":
    Bedgraphfile, flagstatfile, readcountcorrBedgraphfile = sys.argv[1],sys.argv[2], sys.argv[3]
    main(Bedgraphfile, flagstatfile, readcountcorrBedgraphfile)
