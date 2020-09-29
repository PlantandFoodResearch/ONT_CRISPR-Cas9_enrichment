#!/usr/bin/echo source me dont execute me

#'
#'#' Rsamtoools unit tests: https://github.com/Bioconductor/Rsamtools/tree/master/inst/unitTests
#'

#'
#'#' bamCoverage by `which` region 
#'
bamCoverage <- function(index, which, bamFile, orient="all", ...) {
  which_element <- which[index]
  param <- ScanBamParam(what=c('flag', 'pos', 'qwidth'), which=which_element,
                        flag=scanBamFlag(isUnmappedQuery=FALSE))
  x <- scanBam(bamFile, ..., param=param)
  if(orient == "all") {
    ir <- IRanges(x[[1]][["pos"]], width=x[[1]][["qwidth"]])
  } else if (orient == "fwd") {
    ind <- which(x[[1]][["flag"]]==0)
    ir <- IRanges(x[[1]][["pos"]][ind], width=x[[1]][["qwidth"]][ind])
  } else if(orient == "revcomp") {
    ind <- which(x[[1]][["flag"]]==16)
    ir <- IRanges(x[[1]][["pos"]][ind], width=x[[1]][["qwidth"]][ind])
  } else {
    stop(" Arg 'orient' options: all | fwd | revcomp")
  }
  if(length(ir)) {
    seqname <- as.character(GenomicRanges::seqnames(which_element))
    start   <- min(start(ir))
    end     <- max(end(ir))
    return(list(seqname=seqname, start=start, end=end, cov=coverage(ir, shift= -start+1), ir=ir))
  } else {
    return(NULL)
  }
}

#'
#' Extract ensembl coord system
#'
which_coords <- function(which){
  paste0(GenomicRanges::seqnames(which), ":", GenomicRanges::start(which),"-", GenomicRanges::end(which))
}

#'
#' Split ensembl coords
#'
split_coords <- function(x){
  parts <- strsplit(x, c("[:-]"))[[1]]
  y     <- list(seqname=parts[1], start=as.integer(parts[2]), end=as.integer(parts[3]))
  return(y)
}

#'
#'#' KRISPR rna oligos
#'
krispr_rna_oligos <- function(){
  rna <- c(crRNA_RF_1_F = "GTCATATCTAAGGACCCGCGTGG",
           crRNA_RF_2_F = "TCTGTACTCCGTCTGTCGGTCGG",  ## SNP at end of guide decreases efficiency - also interferes with Cas9 RE
           crRNA_RF_3_R = "AGAAGACTGTCAATCCCGAGTGG",
           crRNA_RF_4_F = "TGTCTGGAAAGTTTCTAACGCGG"   ## SNP at beginning of guide descreases efficiency
  )
  return(DNAStringSet(rna))
}

#'
#'#' Median quality scores
#'
filter_quality <- function(index, FastqQual, fun=mean) {
  if(missing(index)){
    ## Process all if index missing
    index <- seq(length(FastqQual))
  }
  FastqQual_sub <- FastqQual[index]
  metric <- rep(NA, length(FastqQual_sub))

  for(i in seq(length(FastqQual_sub))) {
    metric[i] <- fun(as(FastqQual_sub[i], "numeric"))
  }

  return(metric)
}

split_ids <- function(x) {
  unname(sapply(as.character(x),
                         function(x)strsplit(x, " ")[[1]][1]))
}

ont_duplicates <- function(x) {
  stopifnot(class(x)=="ShortReadQ")

  ont_ids <- split_ids(ShortRead::id(x))
  dups    <- Biostrings::duplicated(ShortRead::id(x))

  out     <- list()
  out$all_ids <- ont_ids

  if(any(dups)) {
    out$dup_ind  <- Biostrings::duplicated(ShortRead::id(x))
    out$dup_ids  <- ont_ids[out$dup_ind]
  }

  return(out)
}

basecall_stats_versions <- function(){

  #
  ## Checking if statistics summaries are identical between runs - _no they are not_
  #
  albacore_summary_hramwd <- read.table("/workspace/hramwd/github/analysis-workflows/Malus/Red_Flesh_ON/albacore2/Red_flesh_ON_run1_Cas9/all_summary_run1.txt", header=TRUE, as.is = TRUE)
  albacore_summary_hrpelg <- read.table("/workspace/hrpelg/Red_Flesh_ON/albacore2/Red_flesh_ON_run1_Cas9/all_summary_run1.txt", header=TRUE, as.is = TRUE)

  cat("[ hramwd run albacore dimensions ]\n")
  print(dim(albacore_summary_hramwd))
  cat("[ hrpelg run albacore dimensions ]\n")
  print(dim(albacore_summary_hrpelg))

  expect_equal(names(albacore_summary_hramwd), names(albacore_summary_hrpelg))

  ## Returned row order is different
  o1 <- order(albacore_summary_hramwd$read_id)
  o2 <- order(albacore_summary_hrpelg$read_id)

  plot(o1, o2, xlab="hramwd run order", ylab="hrpelg run order", main="Albacore runs (log scale)", pch=16, cex=0.3, log="xy")
}
