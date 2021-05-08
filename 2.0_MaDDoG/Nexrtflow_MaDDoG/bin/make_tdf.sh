1. Split the bam file into 4 files, two forward strand and reverse strand
    A. the two "forward strand" files which I'd call the "forward strand RNA" files
        i.second read on forward strand
        ii. first read on reverse strand
    B. the two "reverse strand" files which I'd call the "reverse strand RNA" files
        i. second read on reverse strand
        ii. first read on forward strand
2. Combine the forward strand bams, and the reverse strand bams...
   forward.bam (all reads from forward strand RNAs, regardless of the read strand)
   reverse.bam (all reads from reverse strand RNAs, regardless of the read strand)
3. bedtools genomecov each file (without the strand flag!)
    A. for the reverse.bam, flip the counts to a negative number
4. Cat the two bedgraph files to sample.bedgraph
5. Sort the sample.bedgraph file
6. run millionsmapped correction on this.
7. run igvtools on the millionsmapped corrected bedgraph


I fix the problem a different way.

1. Split the bam file into first and second reads
2. genomecov each strand of both bed files while flipping the negative strand counts to a negative number (which is read 1 on the positive strand and read2 on the negative strand)
3. Create 2 bedgraph-like files, one for the positive strand transcripts and one for the negative strand transcripts. These files will have counts in columns 4 and 5 (read 1 counts are in column 4 and read 2 counts are in column 5)
4. Sum the 4th and 5th columns
5. Cat the positive and negative summed files
6. Sort that file.
7. run millionsmapped correction on this
8. run igvtools on the millionsmapped corrected bedgraph

scripts for my method below

Either of those two ways work to fix the issue. 



#http://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_only_the_first_read_from_paired_end_BAM_files
${samtoolsdir}samtools view -h -b -f 0x0040 ${od}${bamroot}.sorted.bam > ${bedgraphfortdfdir}${bamroot}.pairfirst.bam
# 0x0040 is hexadecimal for 64 (i.e. 16 * 4), which is binary for 1000000, corresponding to the read in the first read pair.


#need to know the flag for the second strand
# Jess gave me this https://broadinstitute.github.io/picard/explain-flags.html
#128 means second in pair
#128 in hexadecimal is 0x0080
${samtoolsdir}samtools view -h -b -f 0x0080 ${od}${bamroot}.sorted.bam > ${bedgraphfortdfdir}${bamroot}.pairsecond.bam

#That should get me the two parts of pair separate
#Then I need to run genomecoverage and swap the strand info on the second pair


#Then I need to run genomecoverage on each of them
${bedtoolsdir}genomeCoverageBed -bg -split -strand - -ibam ${bedgraphfortdfdir}${bamroot}.pairfirst.bam -g $genome >$outdir/genomecoveragebedpair/${entry}.pairfirst.pos.bed

${bedtoolsdir}genomeCoverageBed -bg  -split -strand + -ibam ${bedgraphfortdfdir}${bamroot}.pairfirst.bam -g $genome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> $outdir/genomecoveragebedpair/${entry}.pairfirst.neg.bed

${bedtoolsdir}genomeCoverageBed -bg -split -strand + -ibam  ${bedgraphfortdfdir}${bamroot}.pairsecond.bam -g $genome > $outdir/genomecoveragebedpair/${entry}.pairsecond.pos.bed

${bedtoolsdir}genomeCoverageBed -bg -split -strand - -ibam ${bedgraphfortdfdir}${bamroot}.pairsecond.bam -g $genome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> $outdir/genomecoveragebedpair/${entry}.pairsecond.neg.bed

#first I need to sort the Bedgraphs
${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairfirst.pos.bed > ${bedgraphfortdfdir}${bamroot}.pairfirst.pos.BedGraph.sort

${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairfirst.neg.bed > ${bedgraphfortdfdir}${bamroot}.pairfirst.neg.BedGraph.sort


${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairsecond.pos.bed > ${bedgraphfortdfdir}${bamroot}.pairsecond.pos.BedGraph.sort

${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairsecond.neg.bed > ${bedgraphfortdfdir}${bamroot}.pairsecond.neg.BedGraph.sort

#Then I need to add two Bedgraphs

#this should put the values in columns 4 and 5


${bedtoolsdir}unionBedGraphs -i ${bedgraphfortdfdir}${bamroot}.pairfirst.pos.BedGraph.sort ${bedgraphfortdfdir}${bamroot}.pairsecond.pos.BedGraph.sort >${bedgraphfortdfdir}${bamroot}.pos.Bedgraphcol

${bedtoolsdir}unionBedGraphs -i ${bedgraphfortdfdir}${bamroot}.pairfirst.neg.BedGraph.sort ${bedgraphfortdfdir}${bamroot}.pairsecond.neg.BedGraph.sort >${bedgraphfortdfdir}${bamroot}.neg.Bedgraphcol

#then I need to sum cols 4 and 5

awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' ${bedgraphfortdfdir}${bamroot}.pos.Bedgraphcol >${bedgraphfortdfdir}${bamroot}.pos.Bedgraph

awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' ${bedgraphfortdfdir}${bamroot}.neg.Bedgraphcol >${bedgraphfortdfdir}${bamroot}.neg.Bedgraph

#then I need to cat the two Bedgraphs

cat $outdir/genomecoveragebedpair/${entry}.pos.Bedgraph ${bedgraphfortdfdir}${bamroot}.neg.Bedgraph >${bedgraphfortdfdir}${bamroot}.bed

#then I need to sort the final Bedgraph so it can be divided by millions mapped and converted into tdf
${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.bed >${bedgraphfortdfdir}${bamroot}.Bedgraph