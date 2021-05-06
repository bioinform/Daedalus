library(data.table)

dbFiles <- list.files("/sc1/groups/pls-redbfx/users/dannebar/immunoPETE/immunoDB/database_files", full=TRUE)

dbFiles <- dbFiles[grepl("gdna", dbFiles)]
Vdb <- dbFiles[grepl("_V_", dbFiles)]
Jdb <- dbFiles[grepl("_J_", dbFiles)]
Vgenes <- do.call(rbind, lapply(Vdb, fread))
Jgenes <- do.call(rbind, lapply(Jdb, fread))

write.csv(Vgenes, "ImmunoDB_V_gdna_reference.csv", quote=FALSE, row.names=FALSE)
write.csv(Jgenes, "ImmunoDB_J_gdna_reference.csv", quote=FALSE, row.names=FALSE)

#####
#####

Vp <- fread("Vgene_V2_primers.csv")
Jp <- fread("Vgene_V2_primers.csv")


Vp

Vprim <- data.table(sequence_id = Vp$symbol,
                    sequence = Vp$primer_seq)
Jprim <- data.table(sequence_id = Jp$symbol,
                    sequence = Jp$primer_seq)

write.csv(Vprim, "Vgene_V2_primer_reference.csv", quote=FALSE, row.names=FALSE)
write.csv(Jprim, "Jgene_V2_primer_reference.csv", quote=FALSE, row.names=FALSE)

#####
#####

samples <- fread("Expt56_samples.csv")
params <- fread("Ipete_params.csv")


params$umi1 <- ""
params$umi2 <- "NNNNNNNNN"
params$umi_type <- "R2"

params$subsample <- "100000"
params$subsample <- "200000"

params$trim_primers <- 'true'
params$checkSpikein <- 'false'
params$umi_mode <- 'true'

params$primerRef <- '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/V2.2_primers.csv'
params$geneRef <- '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/immunoDB_cdna.csv'
params$spikeRef <- '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/synthetic_seqs_fixed.fasta'
params$vPrimerRef <- '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/Vprimers.fasta'
params$jPrimerRef <- '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/Jprimers.fasta'
params$vGeneRef <- '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/Vgenes_cdna.fasta'
params$jGeneRef <- '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/Jgenes_cdna.fasta'

samples$sample_sheet <-  "/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/iPETEV2_Expt56_LibPool11_NexSeq_Run7_031320_MANIFEST.csv"
samples$project <- "OmicsBricks_ipeteV2_Expt56_testRun"

samples$run_folder <- "/sc1/raw/illumina/nextseq/200313_NB551443_0103_AH3TJ5AFX2"
manifest <- data.table(samples, params)
write.csv(manifest, "Expt56_V2_manifest.csv", quote=FALSE, row.names=FALSE)

manifest

library(data.table)
library(Biostrings)
J <- fread("Jgene_V2_primers.csv")
V <- fread("Vgene_V2_primers.csv")

Vprim <- DNAStringSet(V$primer_seq)
names(Vprim) <- V$symbol
Vprim <- Vprim[!duplicated(Vprim)]
writeXStringSet(Vprim, "Vprimers.fasta")
Jprim <- DNAStringSet(J$primer_seq)
names(Jprim) <- J$symbol
Jprim <- Jprim[!duplicated(Jprim)]
writeXStringSet(Jprim, "Jprimers.fasta")


primers <- rbind(V, J)
primers <- primers[!duplicated(primer_seq)]



primerData <- data.table(sequence_id = primers$symbol, sequence = primers$primer_seq)
write.csv(primerData, "V2.2_primers.csv", quote=FALSE, row.names=FALSE)



write.csv(primers, "VJprimers_V2.csv")


##sample,run_folder,project,sample_sheet,checkSpikeIn


#####
igv <- fread("IG_V_cdna_reference.csv")
trv <- fread("TR_V_cdna_reference.csv")
vgenes <- rbind(igv, trv)

Vseqs <- DNAStringSet(vgenes$sequence)
names(Vseqs) <- vgenes$sequence_id
writeXStringSet(Vseqs, "Vgenes_cdna.fasta")

igj <- fread("IG_J_cdna_reference.csv")
trj <- fread("TR_J_cdna_reference.csv")
jgenes <- rbind(igj, trj)

Jseqs <- DNAStringSet(jgenes$sequence)
names(Jseqs) <- jgenes$sequence_id
writeXStringSet(Jseqs, "Jgenes_cdna.fasta")

