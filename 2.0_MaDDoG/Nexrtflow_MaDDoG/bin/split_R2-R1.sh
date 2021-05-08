set -ue



# Read this http://thegenomefactory.blogspot.com/2013/08/paired-end-read-confusion-library.html
# then this https://davetang.org/muse/2014/07/24/getting-started-paired-end-reads/
# split script https://www.biostars.org/p/92935/

# Get the bam file from the command line
DATA=$1

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 $DATA > fwd1.bam
samtools index fwd1.bam

samtools view -b -f 80 $DATA > fwd2.bam
samtools index fwd2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f fwd.bam fwd1.bam fwd2.bam
samtools index fwd.bam

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 $DATA > rev1.bam
samtools index rev1.bam

samtools view -b -f 64 -F 16 $DATA > rev2.bam
samtools index rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f rev.bam rev1.bam rev2.bam
samtools index rev.bam