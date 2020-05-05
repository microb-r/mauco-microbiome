# Libraries

library(dada2)
library(tidyverse)
library(DECIPHER)

# Set working directory
setwd(".")

# First step, sequence filtering

# Fastq path
fastq_path <- "../data/raw/fastq"

# reads
fnFs <- sort(list.files(fastq_path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(fastq_path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Get the sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Path to the clean fastq
clean_fastq <- "../data/processed/clean_fastq"

filtFs <- file.path(clean_fastq, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(clean_fastq, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Clean reads
clean_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, minLen=50,
                     trimLeft = c(19,20))

# Calculate error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Run Dada2
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Generate seqtab
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Generate output table with sequences
getN <- function(x) sum(getUniques(x))
track <- cbind(clean_out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names


# Save objects
write.table(track, file= "../results/tables/SequenceProcessing.tsv",
            sep="\t", row.names=TRUE, quote = FALSE)
saveRDS(seqtab.nochim, "../data/processed/seqtab.nochim.RDS")

### Assign taxonomy

# Get the reference information
load("../data/references/SILVA_SSU_r138_2019.RData")

# Obtain DNA strings from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

# Assign taxonomy using IDTAXA
ids <- IdTaxa(dna, trainingSet, strand="both", processors=20, verbose=TRUE)

#Ranks to store
ranks <- c("domain", "phylum", "class", "order", "family", "genus")

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
colnames(taxid) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Add species information to the table
taxid <- addSpecies(taxid, "../data/references/silva_species_assignment_v138.fa.gz", 
                    tryRC=TRUE, verbose=TRUE, allowMultiple = TRUE)

# Save the taxid object
saveRDS(taxid, "../data/processed/taxid.rds")

##### 5- Generate Phylogenetic Tree

## Generate phylogenetic tree

seqs <- getSequences(seqtab.nochim)
write(seqs, "../data/processed/seqtab_nochim.list")
system("Align_Tree_Dada2otus.py -i ../data/processed/seqtab_nochim.list -o ../data/processed/seqtab_nochim")








