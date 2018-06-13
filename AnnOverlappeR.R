#!/usr/bin/env Rscript

# libraries 1
suppressMessages(library(GetoptLong))

# optional VARs
numcores <- 7
download <- TRUE
tool_path <- "~/projects/ensembl2ncbi/"
output_path <- "./"
VERSION = "1.1"
latest_NCBI = FALSE
latest_name = "none"

if(length(commandArgs(TRUE)) == 0){
  stop("No arguments found: Rscript AnnOverlappeR.R --help for details\n")
}

GetoptLong(
  "species=s", "Species in scientific name (e.g.: homo_sapiens), string, mandatory option",
  "output_path=s", "output_path, string, optional option",
  "numcores=i", "Number of cpus , integer, optional (default 4)",
  "download!", "--no-download for no downloads, boolean, mandatory option (default TRUE)",
  "latest_ncbi!", "--no-latest_NCBI for no downloads, boolean, mandatory option (default FALSE)",
  "latest_name=s", "latest name found for ncbi latest assembly, string, optional option",
  "verbose",  "print messages",
  head = 'Rscript AnnOverlappeR.R ',
  foot = 'Please contact jochen.bick@usys.ethz.ch for comments'
)


# libraries 2
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(XML))
suppressMessages(library(pbmcapply))

# check if tool_path has "/" at the end
output_path <- paste0(output_path, "/")
output_path <- gsub(pattern = "//", "/", output_path)
output_path <- normalizePath(output_path)

#download  <- FALSE
cat("#######################\n")
cat("###::AnnOverlappeR::###\n")
cat("#######################\n\n")

cat ("INFO::\t\t", "Download boolean set to: ", download, "\n")
cat ("INFO::\t\t", "Number of used cores: ", numcores, "\n")
cat ("INFO::\t\t", "Files will be saved to: ", output_path, "\n")

source(paste0(tool_path, "functions_1.1.R"))

# EXAMPLE:
# species:
# species <- "ovis_aries"
# species <- "sus_scrofa"
# species <- "homo_sapiens"
# species <- "bos_taurus"
# species <- "equus_caballus"
# species <- "canis_familiaris"
# species <- "oryctolagus_cuniculus"
# latest_NCBI <- T
# latest_name <- "Canis_lupus"
# output_path <- "~/data/AnnOverlapper/2017-09-13/"
# download  <- FALSE
# download  <- TRUE
# numcores <- 7



species <- gsub("(.*)", "\\L\\1", species, perl=T) # lower case it
cat("INFO::\t\t", "Species: ", species, "\n\n")

# read refseq and ensembl files (not yet in gff or gtf format)
# change dir for download of gtf and gff of spezific species
# setwd(paste0("data/",species,"/"))
species_path <- paste0(output_path, "/data/", species, "/")
dir.create(path=species_path, recursive=T, showWarnings = FALSE)

# retrieve ensembl data

if(download){
  ens_species.gz <- paste0(species, ".gtf.gz")
  ens_url <- paste0("ftp://ftp.ensembl.org/../pub/current_gtf/", species, "/*[0-9].gtf.gz")
  cat("STATUS::\t", "Download GTF files from ensembl:", ens_url, " ...")
  system(paste0("wget ", ens_url, " -q -O ", species_path, ens_species.gz))
  system(paste0("gunzip -f ", species_path, "*.gtf.gz"))
  cat("done\n")
}
#check number of files in data/species folder for import
ens_file <-  paste0(species_path, species, ".gtf")


cat("STATUS::\t", "Translate chromosome names...")
system(paste0(tool_path, " perl src/chromosome_translater.pl -gff ", ens_file, " -out ", 
              species_path, species, "_ens_mod.gtf"))
cat("done\n")

# clean
system(paste0("rm -f ", ens_file))
system(paste0("rm -f ", species_path, "*.gz"))

ens_species <- paste0(species, "_ens_mod.gtf") # update
ens_file <- paste0(species_path, ens_species)

# shorten ens file
system(paste0("grep -P '\t[CDSexonge]+\t\' ", ens_file , " >", species_path, "data.gtf"))
system(paste0("grep -P '^#' ", ens_file , " >", species_path, "head.gtf"))
system(paste0("cat ", species_path, "head.gtf ", species_path, "data.gtf >", ens_file))

# clean
system(paste0("rm -f ", species_path, "data.gtf"))
system(paste0("rm -f ", species_path, "head.gtf"))

cat("STATUS::\t", "Import GTF file...")

ens_all.gtf <- import.gff(ens_file, format="gtf") 
ens.gtf <- ens_all.gtf

# clean
# system(paste0("rm -f ", ens_file))

cat("done\n")

# subset gtf
ens.gtf@elementMetadata <- ens.gtf@elementMetadata[c("type","gene_id", "gene_name", "gene_biotype", "transcript_id", "protein_id")]#, "source", transcript_id" )] ### cutoff


##### retrieve ncbi data

UC_species <- gsub("(^\\w)(.*)", "\\U\\1\\L\\2", species, perl=T)
if(latest_name != "none"){
  UC_species <- latest_name
}

if(download){
  ncbi_species.gz <- paste0(species, ".gff3.gz")
  if(latest_NCBI){
    ncbi_url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/", UC_species, "/latest_assembly_versions/")
    system(paste0("wget ", ncbi_url, " -q -O ", species_path, "index.html")) # may specifiy a download name
    rawHTML <- paste(readLines(paste0(species_path, "index.html")))
    assemblies <- grep(pattern = "GCF", value = TRUE, x = sub(pattern = ".*latest_assembly_versions/(.*)\".*", replacement = "\\1", x = rawHTML))[1]
    system(paste0("rm -f ", paste0(species_path, "index.html")))
    ncbi_url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/", UC_species, "/latest_assembly_versions/", assemblies, "/*.gff.gz")
  }else{
    ncbi_url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/", UC_species, "/GFF/ref*_top_level.gff3.gz")
    
  }
  cat("STATUS::\t", "Download GFF files from ncbi:", ncbi_url, " ...")
  system(paste0("wget ", ncbi_url, " -q -O ", species_path, ncbi_species.gz)) # may specifiy a download name
  system(paste0("gunzip -f ", species_path, "*.gff3.gz"))
}

ncbi_species <- paste0(species, ".gff3")

# wget chomosome translate table
chr_path <- paste0("/tmp/", species, "/chrtranslation/")
dir.create(path=chr_path, recursive=T, showWarnings = FALSE)

UC_species <- gsub("(^\\w)(.*)", "\\U\\1\\L\\2", species, perl=T)

if(download){
  ncbi_chr_url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/", UC_species, "/Assembled_chromosomes/*")
  system(paste0("wget ", ncbi_chr_url, " -N -q -P ", chr_path))
  cat("done\n")
}
chr_list <- list.files(path=chr_path, pattern="*")
chr_files <- gsub(", ", ",", toString(paste0(chr_path, chr_list))) # take all!!


cat("STATUS::\t", "Translate chromosome names...")

ncbi_file <- paste0(species_path, ncbi_species)

system(paste0(tool_path, " perl src/chromosome_translater.pl -gff ", ncbi_file, " -chr ", chr_files,
       " -out ", species_path, species, "_ncbi_mod.gff3"))

cat("done\n")


# clean
system(paste0("rm -f ", ncbi_file))
system(paste0("rm -f ", species_path, "*.gz"))

ncbi_species <- paste0(species, "_ncbi_mod.gff3") # update
ncbi_file <- paste0(species_path, ncbi_species)

# shorten ncbi file
system(paste0("grep -P '\t[CDSexonge]+\t\' ", ncbi_file , " >", species_path, "data.gff"))
system(paste0("head -n 20 ", ncbi_file ," | grep -P '^#'", " >", species_path, "head.gff"))
system(paste0("cat ", species_path, "head.gff ", species_path, "data.gff >", ncbi_file))

# clean
system(paste0("rm -f ", species_path, "data.gff"))
system(paste0("rm -f ", species_path, "head.gff"))

cat("STATUS::\t", "Import GFF file...")

ncbi_all.gff <- import.gff(ncbi_file, format="gff3")

ncbi.gff  <- ncbi_all.gff

ncbi.gff@elementMetadata <- ncbi.gff@elementMetadata[c("type","geneid", "gene", "transcript_id", "protein_id")] ### do not take description or product because its only for mRNAs


# clean
#system(paste0("rm -f ", ncbi_file))

cat("done\n")


GG1 <- ncbi.gff[ncbi.gff$type %in% "gene", ]
GG1@elementMetadata <- GG1@elementMetadata[c("type","geneid", "gene")] # make it short again
GG1_all <- ncbi_all.gff[ncbi_all.gff$type %in% "gene", ]
GG2 <- ens.gtf[ens.gtf$type %in% "gene", ]

cat(paste0("INFO::\t\t","AnnOverlappeR retrieved: ", length(GG1), " genes from NCBI and ", length(GG2), " genes from Ensembl!\n"))

save(ens_all.gtf, file=paste0(species_path, "ens_all.gtf.RData"))
save(ncbi_all.gff, file=paste0(species_path, "ncbi_all.gff.RData"))

save(GG1, file=paste0(species_path, "GG1.RData"))
save(GG2, file=paste0(species_path, "GG2.RData"))



###################
# ncbi VS ensembl #
###################

cat("STATUS::\t", "Check gene overlapping...")

results <- makegenehits(GG1, GG2)

result_genes <- results$finalresult
no_result_genes <- results$no_result



final.per.order <- result_genes[order(result_genes$percent, decreasing=T, result_genes$geneid1), ]

write.table(x=final.per.order, file=paste0(species_path, "final.per.order.xls"), row.names=F, sep="\t")
save(final.per.order, file=paste0(species_path, "final.per.order.RData"))

cat("done\n")

cat(paste0("INFO::\t\t","AnnOverlappeR calculated: ", nrow(result_genes), " overlapping genes!\n"))


##### exons

cat("STATUS::\t", "Check exon overlapping...")
# read in exons
eGG1 <- ncbi.gff[ncbi.gff$type %in% "exon", ]
eGG1@elementMetadata <- eGG1@elementMetadata[c("type","geneid", "gene", "transcript_id")] # make it short again
eGG2 <- ens.gtf[ens.gtf$type %in% "exon", ]

# static VARs
# exon1=eGG1
# exon2=eGG2

save(eGG1, file=paste0(species_path, "eGG1.RData"))
save(eGG2, file=paste0(species_path, "eGG2.RData"))

# CDS checker (new)
cGG1 <- ncbi.gff[ncbi.gff$type %in% "CDS", ]
cGG1@elementMetadata <- cGG1@elementMetadata[c("type","geneid", "gene", "protein_id")] # make it short again
cGG2 <- ens.gtf[ens.gtf$type %in% "CDS", ]

save(cGG1, file=paste0(species_path, "cGG1.RData"))
save(cGG2, file=paste0(species_path, "cGG2.RData"))


cat(" starting multicore...")


# static vars
ref_name = final.per.order$geneid1
ens_name = final.per.order$gene_id2

result_exon_all <- pbmclapply(seq_along(ref_name), function(x){
  result_exon <- exoncheck(ref_name[x], ens_name[x])#, F) # changed because of overhead problems
}, mc.cores = numcores)

# prepare no_result data
no_result_genes$less_splices_ex <- no_result_genes$noofhits <- no_result_genes$percent_ss <- no_result_genes$percent.mean_exon <- no_result_genes$percent.mean_cds <- 0

# save no result genes
save(no_result_genes ,file=paste0(species_path, "no_result_genes.RData"))
write.table(x=no_result_genes , file=paste0(species_path, "no_result_genes.xls"), row.names=F, quote = F, sep="\t")

save(result_exon_all, file=paste0(species_path, "result_exon_all.RData"))

exon.df <- do.call(rbind, result_exon_all)


cat("done\n")

final_ge.df <- cbind(final.per.order, exon.df) 
write.table(x=final_ge.df , file=paste0(species_path, "final-gene_exon.xls"), row.names=F, sep="\t")
gid_gs.df <- data.frame(geneid=GG1_all@elementMetadata$geneid, genesymbol=GG1_all@elementMetadata$gene)
write.table(x=gid_gs.df , file=paste0(species_path, "gid_gs.xls"), row.names=F, col.names=F, quote = F, sep="\t")

save(final_ge.df, file=paste0(species_path, "final_ge.RData"))



cat("STATUS::\t", "Filter results and save them for DB...")


filter_easy50 <- filter_easy(results.df = final_ge.df, cutoff = 50)
save(filter_easy50, file=paste0(species_path, "filter_easy50.RData"))
write.table(x=filter_easy50 , file=paste0(species_path, "filter_easy50.xls"), row.names=F, quote = F, sep="\t")


final_ge.df_filter50 <- filter_df(results.df = final_ge.df, cutoff = 50)
write.table(x=final_ge.df_filter50 , file=paste0(species_path, "fin_filter50.xls"), row.names=F, quote = F, sep="\t")

filter_50ge <- final_ge.df_filter50

filter_50ge$unique <- paste0(filter_50ge$geneid1, ":", filter_50ge$gene_id2) 

# check order
filter_50ge <- filter_50ge[order(filter_50ge$unique, decreasing=F), ]

filter_50ge_nodup <- filter_50ge[!duplicated(filter_50ge$unique), ] # get rid of duplications

write.table(x=filter_50ge_nodup , quote=F, file=paste0(species_path, "filter_50ge_nodup.txt"), row.names=F, sep="\t", col.names=F)

save(filter_50ge_nodup, file=paste0(species_path, "filter_50ge_nodup.RData"))

# change digets rounding
filter_50ge_nodup$percent  <- round(filter_50ge_nodup$percent, digits = 3)
filter_50ge_nodup$percent.mean_exon  <- round(filter_50ge_nodup$percent.mean_exon, digits = 3)
filter_50ge_nodup$percent_ss  <- round(filter_50ge_nodup$percent_ss, digits = 3)
filter_50ge_nodup$percent.mean_cds  <- round(filter_50ge_nodup$percent.mean_cds, digits = 3)

cat("done!\n")

cat(paste0("INFO::\t\t", nrow(result_genes)-nrow(filter_50ge_nodup), " did not pass filter criteria!\n"))

ent2ens <- filter_50ge_nodup[c("geneid1", "gene_id2", "gene1", "percent", "percent.mean_exon", "less_splices_ex", "percent_ss")]

save(ent2ens, file=paste0(species_path, "ent2ens.RData"))

write.table(x=filter_50ge_nodup[c("geneid1", "gene_id2")] , quote=F, file=paste0(species_path, "ent2ens.txt"), 
            row.names=F, sep="\t", col.names=F)

write.table(x=filter_50ge_nodup[c("gene_id2", "geneid1", "gene1", "gene_biotype2", "percent",
                                  "percent.mean_exon", "less_splices_ex", "percent_ss",
                                  "percent.mean_cds")] ,
            quote=F, file=paste0(species_path, "ensembl2entrez_for_db_all.txt"), row.names=F, sep="\t", col.names=F)

cat(paste0("INFO::\t\t","AnnOverlappeR calculated: ", nrow(filter_50ge_nodup), " identifier pairs!\n"))

## easy 50 part
filter_easy50ge <- filter_easy50

filter_easy50ge$unique <- paste0(filter_easy50ge$geneid1, ":", filter_easy50ge$gene_id2) 

# check order
filter_easy50ge <- filter_easy50ge[order(filter_easy50ge$unique, decreasing=F), ]

filter_easy50ge_nodup <- filter_easy50ge[!duplicated(filter_easy50ge$unique), ] ### get rid of duplications

cat("INFO::\t\t", nrow(filter_easy50ge)-nrow(filter_easy50ge_nodup), " we found to be duplicated!\n")

# change digets rounding
filter_easy50ge_nodup$percent  <- round(filter_easy50ge_nodup$percent, digits = 3)
filter_easy50ge_nodup$percent.mean_exon  <- round(filter_easy50ge_nodup$percent.mean_exon, digits = 3)
filter_easy50ge_nodup$percent_ss  <- round(filter_easy50ge_nodup$percent_ss, digits = 3)
filter_easy50ge_nodup$percent.mean_cds  <- round(filter_easy50ge_nodup$percent.mean_cds, digits = 3)
#names(filter_50ge_nodup)
cat("done!\n")

write.table(x=filter_easy50ge_nodup[c("gene_id2", "geneid1", "gene1", "gene_biotype2",
                                      "percent", "percent.mean_exon", "less_splices_ex",
                                      "percent_ss", "percent.mean_cds")] ,
            quote=F, file=paste0(species_path, "ensembl2entrez_easy50_for_db_all.txt"), row.names=F, sep="\t", col.names=F)

cat(paste0("INFO::\t\t","AnnOverlappeR calculated: ", nrow(filter_easy50ge_nodup), " identifier pairs with the easy 50 filter!\n"))


cat("STATUS::\t", "COMPLETED successfully!!!!\n")



################################
############# END ##############
################################

