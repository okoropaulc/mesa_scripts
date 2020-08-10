library(dada2) # version '1.11.3' #obtained with packageVersion("dada2")
path <- "Z:/paul/mets61/"
#list.files(path)

# Forward and reverse fastq filenames have format: SampleID_F.fastq and SampleID_R.fastq
fnFs <- sort(list.files(path, pattern="_F.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R.fastq", full.names = TRUE))

# Extract sample IDs, assuming filenames have format: SampleID_F.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Examine quality profiles of forward and reverse reads
#plotQualityProfile(fnFs[1:4])
#plotQualityProfile(fnRs[1:4])

mypath <- "/Users/Okoro/OneDrive/Desktop/Thesis_Papers/mets61" #just to save the filtered data to my laptop

#Perform filtering and trimming
filt_path <- file.path(mypath, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,175), maxN=0, matchIDs = TRUE, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
#View how many reads passed the filter
head(out)

##Examine quality profiles of filtered forward and reverse reads
#pdf("/home/paul/metstry1/Second/qplotTrimmedF.pdf")
plotQualityProfile(filtFs[1:4])
#dev.off()
#pdf("/home/paul/metstry1/Second/qplotTrimmedR.pdf")
plotQualityProfile(filtRs[1:4])
#dev.off()

#Learn the Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Depreplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
#We are now ready to apply the core sequence-variant inference algorithm to the dereplicated data.
#Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#Inspecting the dada-class object returned by dada:
dadaFs[[1]]
head(dadaFs[[1]]$clustering)
#head(dadaRs[[1]]$clustering)

#Merge Paired reads
#Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate = TRUE, verbose=TRUE)
#mergePairs(ddF, drpF, ddR, drpR, minOverlap=4)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct the sequence table
#The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Remove Chimeric Sequence
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Track reads through the pipeline
#As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
#write.table(track, file ="/home/paul/metstry1/Third/mergersReport.txt", sep="\t", quote = F)

#Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/paul/silva_nr_v128_train_set.fa.gz", multithread=TRUE)

#Because I used justConcatenate=TRUE (which do not allow species to be assigned), I have to write a code to specially addSpecie to the sequences
asvSeqs <- gsub("N{10}[A-z]*", "", getSequences(seqtab.nochim))
asgnSpecs <- assignSpecies(asvSeqs, "/home/paul/silva_species_assignment_v128.fa.gz")

# transfer species annotation to same genus ASVs else "NA"
same.genus <- asgnSpecs[,"Genus"] == taxa[,"Genus"]
same.genus[is.na(same.genus)] <- FALSE
taxa <- cbind(taxa, Species = rep(NA, nrow(taxa)))
taxa[same.genus, "Species"] <- asgnSpecs[same.genus, "Species"]
#head(taxa[same.genus, "Species"])

taxa <- addSpecies(taxa, "/home/paul/silva_species_assignment_v128.fa.gz")
#Lets inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) #Prints the taxonomy with species


#Mock Validation
#unqs.mock <- seqtab.nochim["GSF856-Dugas-101-QH2O",]
#unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
#cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

library(phyloseq)
library(ggplot2)

#Phyloseq Tryout
#read sample phenotype information
met61 <- read.table(file = "/home/paul/mets61obesity.txt")
#Order the ID
met61 <- met61[order(met61$id),]
#drop the x column of met62
#met62$X <- NULL

#extract the rownames of OTU table - seqtab.nochim
samples.out <- rownames(seqtab.nochim)
#assign the extracted rownames to be the rownames of met61
rownames(met61) <- samples.out

#Save the the sorted met61 
write.table(met61, file = "/home/paul/mets61/mets61obesity_sorted.txt")

#Ordination
# making our phyloseq object with transformed table
OTUcount <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
met61_sam <- sample_data(met61)
taxonomy <- tax_table(taxa)
#make phyloseq object
ps <- phyloseq(OTUcount, met61_sam, taxonomy)
#visualize aplha-diversity
shannon <- plot_richness(ps, x="id", measures=c("Shannon"), color="Site")
write.table(shannon$data, file = "/home/paul/mets61/mets61_alpha_index_shannon.txt") #save Shannon index values
plot_richness(ps, x="id", measures=c("Shannon"), color="Site", title = "METS Alpha Shannon Index Plot", shape = "Obesity_status")

simpson <- plot_richness(ps, x="id", measures=c("Simpson"), color="Site")
write.table(simpson$data, file = "/home/paul/mets61/mets61_alpha_index_simpson.txt") #save Simpson index values
plot_richness(ps, x="id", measures=c("Simpson"), color="Site", title = "METS Alpha Simpson Index Plot", shape = "Obesity_status")

observed <- plot_richness(ps, x="id", measures=c("Observed"), color="Site")
write.table(observed$data, file = "/home/paul/mets61/mets61_alpha_index_observed.txt") #save Observed index values
plot_richness(ps, x="id", measures=c("Observed"), color="Site", title = "METS Alpha Observed Index Plot", shape = "Obesity_status")

chao1 <- plot_richness(ps, x="id", measures=c("Chao1"), color="Site")
write.table(chao1$data, file = "/home/paul/mets61/mets61_alpha_index_chao1.txt") #save chao1 index values
plot_richness(ps, x="id", measures=c("Simpson"), color="Site", title = "METS Alpha Chao1 Index Plot", shape = "Obesity_status")

ace <- plot_richness(ps, x="id", measures=c("ACE"), color="Site")
write.table(ace$data, file = "/home/paul/mets61/mets61_alpha_index_ace.txt") #save ACE index values
plot_richness(ps, x="id", measures=c("ACE"), color="Site", title = "METS Alpha ACE Index Plot", shape = "Obesity_status")

InvSimpson <- plot_richness(ps, x="id", measures=c("InvSimpson"), color="Site")
write.table(InvSimpson$data, file = "/home/paul/mets61/mets61_alpha_index_InvSimpson.txt") #save InvSimpson index values
plot_richness(ps, x="id", measures=c("InvSimpson"), color="Site", title = "METS Alpha InvSimpson Index Plot", shape = "Obesity_status")

fisher <- plot_richness(ps, x="id", measures=c("Fisher"), color="Site")
write.table(fisher$data, file = "/home/paul/mets61/mets61_alpha_index_fisher.txt") #save fisher index values
plot_richness(ps, x="id", measures=c("Fisher"), color="Site", title = "METS Alpha Fisher Index Plot", shape = "Obesity_status")

#plot_richness(ps, measures=c("Shannon", "Simpson"), color="Site")
#plot_richness(ps, measures=c("Shannon", "Simpson"), color="Obesity_status")
#plot_richness(ps, measures=c("Shannon", "Simpson"), color = "Obesity_status") + aes(shape="Site")

#Ordination
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Site", title="METS Bray NMDS", shape = "Obesity_status")
#plot_ordination(ps.prop, ord.nmds.bray, color="Obesity_status", title = "Bray NMDS")


#Bar Plot for top20 sequences
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Family") + facet_wrap(~Obesity_status, scales="free_x")
plot_bar(ps.top20, fill="Family") + facet_wrap(~Site, scales="free_x")

#write setab.nochim and taxa to csv files
write.csv(seqtab.nochim, "/home/paul/mets61/met61_OTU_table.csv")
write.csv(taxa, "/home/paul/mets61/met61_taxonomy_table.csv")
write.table(seqtab.nochim, "/home/paul/mets61/met61_OTU_table.txt")
write.table(taxa, "/home/paul/mets61/met61_taxonomy_table.txt")





#Extracting the result from R to generate a Fasta file of ASV, a count table and a taxonomy table
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/paul/metstry1/Third/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "/home/paul/mets61/mets61_ASVs_counts.txt", sep="\t", quote=F)

# taxa table with no species assignment:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "/home/paul/mets61/mets61_ASVs_taxonomy_with_no_species.txt", sep="\t", quote=F)

# taxa table with species:
asv_tax <- taxa.print
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "/home/paul/mets61/mets61_ASVs_taxonomy_with_Species.txt", sep="\t", quote=F)

