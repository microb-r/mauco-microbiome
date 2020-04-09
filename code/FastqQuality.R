library(seqTools)

R1_reads <- dir(path="../data/raw/fastq/", pattern="*R1_001.fastq.gz")
R2_reads <- dir(path="../data/raw/fastq/", pattern="*R2_001.fastq.gz")

R1_quality <- fastqq(file.path("../data/raw/fastq", R1_reads), k=6)
R2_quality <- fastqq(file.path("../data/raw/fastq", R2_reads), k=6)

plotMergedPhredQuant(R1_quality)
plotMergedPhredQuant(R2_quality)
