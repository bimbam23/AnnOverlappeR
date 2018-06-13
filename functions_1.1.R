  
#####
# FUNCTIONS
#####


filter_easy <- function (results.df = final_ge.df, cutoff = 50){
  results.df$percent.mean_exon <- as.double(sub(pattern=NaN,replacement=0L,x=results.df$percent.mean_exon))
  results.df$percent.mean_cds <- as.double(sub(pattern=NaN,replacement=0L,x=results.df$percent.mean_cds))

  no_rows_unfilt <- nrow(results.df)
  cat("INFO:: There are", no_rows_unfilt, "unfiltered rows\n")
  
  good_genes <- results.df[results.df$percent >= cutoff, ] # keep genes that have a greater overlap than 50% anyway without looking at exon levels
  good_genes$group <- "50<=gOL"
  no_rows_good_genes <- nrow(good_genes)
  cat("INFO:: There are", no_rows_good_genes, "good genes\n")

  subresults.df <- results.df[results.df$percent < cutoff, ]

  # exon and CDS OL check
  subresults1_in.df <- subresults.df[subresults.df$percent.mean_exon >= cutoff | subresults.df$percent.mean_cds >= cutoff, ] # true exon overlap bigger cutoff
  subresults1_in.df$group <- "50>gOL_50<=ecOL"
  no_rows_subresults1_in.df <- nrow(subresults1_in.df)
  cat("INFO:: There are", no_rows_subresults1_in.df, "rows with 50% exon or CDS OL\n")
  
  subresults1_out.df <- subresults.df[subresults.df$percent.mean_exon < cutoff & subresults.df$percent.mean_cds < cutoff, ]

  filterEasy_50ge <- rbind(good_genes, subresults1_in.df)
  no_rows_filterEasy_50ge <- nrow(filterEasy_50ge)
  cat("INFO:: There are", no_rows_filterEasy_50ge, "rows left after filtering\n")
  
  return(filterEasy_50ge)
}

makegenehits <- function (gGR1=GG1, gGR2=GG2){
  
  no_hit_pos  <- seq_len(length(gGR1))
  hits <- suppressWarnings(findOverlaps(gGR1, gGR2)) # do the gene operlap
  no_hit_pos  <- no_hit_pos[-queryHits(hits)]
  innercut_wd <- suppressWarnings(width(pintersect(gGR1[queryHits(hits)], gGR2[subjectHits(hits)]))) ## width of overlapregion of both genes
  outercut_wd <- suppressWarnings(width(punion(gGR1[queryHits(hits)],gGR2[subjectHits(hits)]))) ## width of both gene regions together
  percentoverlap <- (innercut_wd/outercut_wd)*100 # percent of overlap

  # filtering on outerOL of two genes
  wGR1 <- width(gGR1[queryHits(hits)])
  wGR2 <- width(gGR2[subjectHits(hits)])
  
  wGR <- lapply(seq_along(hits), function(x){
    data.frame(max=max(wGR1[x], wGR2[x]), min=min(wGR1[x], wGR2[x])) # sort and use min an max
  })
  wGR <- do.call(rbind, wGR)
  
  no_innerOL <- outercut_wd > wGR$max # no. genes which are included inside each other
  ol_to_minwGR <- -(outercut_wd[no_innerOL]-wGR$min[no_innerOL]-wGR$max[no_innerOL])/wGR$min[no_innerOL] # percentage of smaller gene overlapping bigger gene
  
  ol_to_minwGR_filtered <- ol_to_minwGR<0.25 # cutoff 25% min overlapp of smaller gene with bigger gene
  no_innerOL[which(no_innerOL)] <- ol_to_minwGR_filtered
  filter <- !no_innerOL 
  
  result <- gGR1[queryHits(hits)]
  no_result <- gGR1[no_hit_pos]
  result <- as.data.frame(result)
  names(result) <- paste0(names(result), 1) # rename first GR
  result$percent <- percentoverlap # set percent gene overlap
  result2 <- gGR2[subjectHits(hits)] 
  result2 <- as.data.frame(result2)
  
  names(result2) <- paste0(names(result2), 2) # rename second GR
  finalresult <- cbind(result, result2)
  # filter step
  finalresult_f <- finalresult[filter, ]

  # add data.frame for no_results
  if(length(no_result) != 0){ 
    no_result <- as.data.frame(no_result)
    names(no_result) <- paste0(names(no_result), 1)
    lnr <- nrow(no_result)
    no_result2 <- data.frame(percent = rep(0, lnr), seqnames2 = rep("-", lnr), start2 = rep(0, lnr),
                             end2 = rep(0, lnr), width2 = rep(0, lnr), strand2 = rep(0, lnr), type2 = rep(0, lnr), 
                             gene_id2 = rep("-", lnr))
    no_result <- cbind(no_result, no_result2)
    return(list(finalresult=finalresult, no_result=no_result, filtered_out=finalresult[!filter, ]))
  }else{
    return(list(finalresult=finalresult, no_result=data.frame(), filtered_out=finalresult[!filter, ]))
  }
}


exoncheck <- function (ref_name, ens_name){
  exon1_sub <- eGG1[eGG1$geneid %in% ref_name]
  exon2_sub <- eGG2[eGG2$gene_id %in% ens_name]

  cds1_sub <- cGG1[cGG1$geneid %in% ref_name]
  cds2_sub <- cGG2[cGG2$gene_id %in% ens_name]
  
  # splice side overlap check
  splices_ex1 <- unique(c(start(ranges(exon1_sub)), end(ranges(exon1_sub))))
  splices_ex2 <- unique(c(start(ranges(exon2_sub)), end(ranges(exon2_sub))))
  
  ss_hits <- as.data.frame(table(table(c(splices_ex1, splices_ex2))))
  less_splices_ex <- min(length(splices_ex1), length(splices_ex2))
  noofhits <- ss_hits[ss_hits$Var1 == "2" , "Freq"]
  percent_ss <- (noofhits/less_splices_ex)*100
  if(identical(noofhits, integer(0))){
    noofhits <- 0L
    percent_ss <- 0L
  }
  percent.mean_exon <- suppressWarnings(mean_ol(GGRange1 = exon1_sub, GGRange2 = exon2_sub))
  percent.mean_cds <- suppressWarnings(mean_ol(cds1_sub, cds2_sub))
  data.frame(less_splices_ex, noofhits, percent_ss, percent.mean_exon, percent.mean_cds)
}


mean_ol <- function(GGRange1, GGRange2){
  hits <- findOverlaps(GGRange1, GGRange2)
  if(length(hits) != 0){
    innercut_wd <- width(pintersect(GGRange1[queryHits(hits)], GGRange2[subjectHits(hits)])) ## width of Overlapregion
    outercut_wd <- width(punion(GGRange1[queryHits(hits)], GGRange2[subjectHits(hits)])) ## width of two generegions
  
    percentoverlap <- (innercut_wd/outercut_wd)*100
    
    inds <- split( 1:length(hits), queryHits(hits) )
    shits.mark <- data.frame(queryHits=names(inds), 
                             percOL=sapply(inds, function(u) pmin(sum(percentoverlap[u]),100)),
                             subjectHits=sapply(inds, function(u) subjectHits(hits)[u[1]]))
 
    hits.df <- shits.mark
    return (mean(hits.df$percOL))
  }else{
    return(0)
  }
  
}

# optional extra filtering functions
filter_df <- function (results.df=final_ge.df, cutoff=50){
  results.df$percent.mean_exon <- as.double(sub(pattern=NaN,replacement=0L,x=results.df$percent.mean_exon))
  results.df$percent.mean_cds <- as.double(sub(pattern=NaN,replacement=0L,x=results.df$percent.mean_cds))
  
  no_rows_unfilt <- nrow(results.df)
  cat("INFO:: There are", no_rows_unfilt, "unfiltered rows\n")
  
  good_genes <- results.df[results.df$percent >= 50, ] # keep genes that have a greater overlap than 50% anyway without looking at exon levels
  good_genes$group <- "50<=gOL"
  no_rows_good_genes <- nrow(good_genes)
  cat("INFO:: There are", no_rows_good_genes, "good genes\n")
  
  results.df <- results.df[results.df$percent < 50, ] # rest!
  
  # number of spliceslite hits
  subresults1.df <- results.df[results.df$noofhits == 0, ] # no splicesite hit -> check for general exon overlap
  subresults2.df <- results.df[results.df$noofhits > 0, ] # splicesite hit
  
  subresults1_in.df <- subresults1.df[subresults1.df$percent.mean_exon >= cutoff, ] # true exon overlap bigger cutoff
  subresults1_in.df$group <- "50>gOL_0==nss_50<=eOL"
  no_rows_subresults1_in.df <- nrow(subresults1_in.df)
  cat("INFO:: There are", no_rows_subresults1_in.df, "rows with this condtion:50>gOL_0==nss_50<=eOL\n")
  
  subresults1_out.df <- subresults1.df[subresults1.df$percent.mean_exon < cutoff, ]
  subresults1b_in.df <- subresults1_out.df[subresults1_out.df$percent.mean_cds >= cutoff, ] # cds filter
  subresults1b_in.df$group <- "50>gOL_0==nss_50<=cOL"
  no_rows_subresults1b_in.df <- nrow(subresults1b_in.df)
  cat("INFO:: There are", no_rows_subresults1b_in.df, "rows with this condtion:50>gOL_0==nss_50<=cOL\n")
  subresults1b_out.df <- subresults1_out.df[subresults1_out.df$percent.mean_cds < cutoff, ]
  
  subresults2_in.df <- subresults2.df[subresults2.df$percent_ss >= cutoff, ] # keep it if 50% of the splicesites are overlapping
  subresults2_in.df$group <- "50>gOL_0<nss_50<=pssOL"
  no_rows_subresults2_in.df <- nrow(subresults2_in.df)
  cat("INFO:: There are", no_rows_subresults2_in.df, "rows with this condtion:50>gOL_0<nss_50<=pssOL\n") 
  subresults2_out.df <- subresults2.df[subresults2.df$percent_ss < cutoff, ] # maybe set this to cutoff var
  
  # second filtering on cutoff 
  subresults3_in.df <- subresults2_out.df[subresults2_out.df$percent.mean_exon >= cutoff, ]
  subresults3_in.df$group <- "50>gOL_0<nss_50>pssOL_50<=eOL"
  no_rows_subresults3_in.df <- nrow(subresults3_in.df)
  cat("INFO:: There are", no_rows_subresults3_in.df, "rows with this condtion:50>gOL_0<nss_50>pssOL_50<=eOL\n")
  subresults3_out.df <- subresults2_out.df[subresults2_out.df$percent.mean_exon < cutoff, ]
  subresults3b_in.df <- subresults3_out.df[subresults3_out.df$percent.mean_cds >= cutoff, ] # cds filter
  subresults3b_in.df$group <- "50>gOL_0<nss_50>pssOL_50<=cOL"
  no_rows_subresults3b_in.df <- nrow(subresults3b_in.df)
  cat("INFO:: There are", no_rows_subresults3b_in.df, "rows with this condtion:50>gOL_0<nss_50>pssOL_50<=cOL\n")
  subresults3b_out.df <- subresults3_out.df[subresults3_out.df$percent.mean_cds < cutoff, ]
  
  filter_ge <- rbind(good_genes, subresults1_in.df, subresults1b_in.df, subresults2_in.df, subresults3_in.df, subresults3b_in.df) #, out)
  no_rows_filter_ge <- nrow(filter_ge)
  cat("INFO:: There are", no_rows_filter_ge, "rows left after filtering\n")
  return(filter_ge)
  lost_genes <- rbind(subresults1_out.df, subresults1b_out.df, subresults2_out.df, subresults3_out.df, subresults3b_out.df) #, out)
  
  filtergroups.list <- list(good_genes = paste0(good_genes$geneid1, ":", good_genes$gene_id2), 
                            no_ss_exon = paste0(subresults1_in.df$geneid1, ":", subresults1_in.df$gene_id2),
                            no_ss_cds = paste0(subresults1b_in.df$geneid1, ":", subresults1b_in.df$gene_id2),
                            bigger50p_ss_cds = paste0(subresults2_in.df$geneid1, ":", subresults2_in.df$gene_id2),
                            smaller50_ss_exon = paste0(subresults3_in.df$geneid1, ":", subresults3_in.df$gene_id2),
                            smaller50_ss_cds = paste0(subresults3b_in.df$geneid1, ":", subresults3b_in.df$gene_id2))
  save(filtergroups.list, file = paste0(species_path, "filtergroup.list.RData"))
  
}

filter_advanced <- function (results.df = final_ge.df, geneCutoff = 50, SSexonCDSCutoff = 50, highOLCutoff = 80, mediumOLCutoff = 70){
  results.df$percent.mean_exon <- as.double(sub(pattern=NaN,replacement=0L,x=results.df$percent.mean_exon))
  results.df$percent.mean_cds <- as.double(sub(pattern=NaN,replacement=0L,x=results.df$percent.mean_cds))
  
  no_rows_unfilt <- nrow(results.df)
  cat("INFO:: There are", no_rows_unfilt, "unfiltered rows\n")
  
  good_genes <- results.df[results.df$percent >= geneCutoff, ] # keep genes that have a greater overlap than 50% anyway without looking at exon levels
  good_genes$group <- "50<=gOL"
  no_rows_good_genes <- nrow(good_genes)
  cat("INFO:: There are", no_rows_good_genes, "good genes\n")
  
  sub_results.df <- results.df[results.df$percent < geneCutoff, ]
  
  # number of spliceslite hits
  subresults1.df <- sub_results.df[sub_results.df$noofhits == 0, ] # no splicesite hit -> check for general exon overlap
  no_rows_subresults1.df <- nrow(subresults1.df)
  cat("INFO:: There are", no_rows_subresults1.df, "rows with no ssOL\n")
  subresults2.df <- sub_results.df[sub_results.df$noofhits > 0, ] # splicesite hit
  no_rows_subresults2.df <- nrow(subresults2.df)
  cat("INFO:: There are", no_rows_subresults2.df, "rows with more than one ssOL\n")
  
  
  subresults1_in.df <- subresults1.df[subresults1.df$percent.mean_exon >= highOLCutoff | subresults1.df$percent.mean_cds >= highOLCutoff, ] # true exon overlap bigger cutoff
  subresults1_in.df$group <- "50>gOL_0==nss_80<=ecOL"
  no_rows_subresults1_in.df <- nrow(subresults1_in.df)
  cat("INFO:: There are", no_rows_subresults1_in.df, "rows with 80% exon or CDS OL\n")
  
  subresults1_out.df <- subresults1.df[subresults1.df$percent.mean_exon < highOLCutoff & subresults1.df$percent.mean_cds < highOLCutoff, ]
  nrow(subresults1_out.df)
  
  subresults2_in.df <- subresults2.df[subresults2.df$percent_ss >= SSexonCDSCutoff, ] # keep it if 50% of the splicesites are overlapping
  subresults2_in.df$group <- "50>gOL_0<nss_50<=pssOL"
  no_rows_subresults2_in.df <- nrow(subresults2_in.df)
  cat("INFO:: There are", no_rows_subresults2_in.df, "rows with 50% splice side OL\n")
  
  subresults2b.df <- subresults2.df[subresults2.df$percent_ss < SSexonCDSCutoff, ]
  
  # second filtering on cutoff
  subresults3_in.df <- subresults2b.df[subresults2b.df$percent.mean_exon >= mediumOLCutoff | subresults2b.df$percent.mean_cds >= mediumOLCutoff, ]
  subresults3_in.df$group <- "50>gOL_0<nss_50>pssOL_50<=ecOL"
  no_rows_subresults3_in.df <- nrow(subresults3_in.df)
  cat("INFO:: There are", no_rows_subresults3_in.df, "rows with 70% exon or CDS OL\n")
  
  subresults3_out.df <- subresults2b.df[subresults2b.df$percent.mean_exon < mediumOLCutoff & subresults2b.df$percent.mean_cds < mediumOLCutoff, ]
  
  
  filter_60ge <- rbind(good_genes, subresults1_in.df, subresults2_in.df, subresults3_in.df)
  filter_60ge
  no_rows_filter_60ge <- nrow(filter_60ge)
  cat("INFO:: There are", no_rows_filter_60ge, "rows left after filtering\n")
  
  lost_genes <- rbind(subresults1_out.df, subresults3_out.df) #, out)
  filtergroups.list <- list(good_genes = paste0(good_genes$geneid1, ":", good_genes$gene_id2), 
                            no_ss_exon = paste0(subresults1_in.df$geneid1, ":", subresults1_in.df$gene_id2),
                            no_ss_cds = paste0(subresults1b_in.df$geneid1, ":", subresults1b_in.df$gene_id2),
                            bigger50p_ss_cds = paste0(subresults2_in.df$geneid1, ":", subresults2_in.df$gene_id2),
                            smaller50_ss_exon = paste0(subresults3_in.df$geneid1, ":", subresults3_in.df$gene_id2),
                            smaller50_ss_cds = paste0(subresults3b_in.df$geneid1, ":", subresults3b_in.df$gene_id2))
  save(filtergroups.list, file = paste0(species_path, "filtergroup.list.RData"))
}
