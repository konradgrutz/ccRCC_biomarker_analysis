 
# filter a set of bam files 
# only allow reads that have at least n x M in cigar string in sum
# e.g. n=30  12S14M2I16M -> 14+16M
# AND at most k edit distance:
#  nM tag contains edit distance, but only counts mismatches for STAR bams
#  add deletions and insertions from CIGAR (STAR) contains: D,I,M,N,S -> N is huge gap -> omit = split read

# WARNING: only uniquely mapped reads, others unsure

library(Rsamtools) # read/write/filter bam
library(stringr) # process cigar string
library(parallel)

args <- commandArgs(trailingOnly = T)

if(length(args) <3) {
  cat("ERROR: at least 3 parameters needed: 
			number of CPUs, 
        least number of M in total in CIGAR for each read (e.g. 30), 
        max percent of InDels/soft-clipping (e.g. 10),
        new suffix for bam files (e.g. -30MCIGAR), 
        bam files\n")
  q("no")
}

nbCPUs <- as.numeric(args[1]) # n.b. rsamtools needs a lot of RAM
minM <- as.numeric(args[2]) # 25
maxNonMPercent <- as.numeric(args[3]) # max 10% of nucleotides may be non M
suffix <- args[4] # "-25M10percIDS-uniqMap"

bamFiles <- args[5:length(args)]
#bamFiles <- list.files(path = ".",patter="Aligned.sortedByCoord.out-25M-CIGAR-uniqMap.bam$")
#bamFiles <- list.files(path = ".",pattern = "Aligned.sortedByCoord.out.bam$",recursive = T)
bamFiles[1:min(length(bamFiles),20)]
bamFiles <- grep("\\/Alig",bamFiles,value=T)
param <- ScanBamParam(what=c("cigar","qwidth"),tag=c("nM","NH"))
bamTagFilter(param) <- list(NH = 1) # unique mapped reads, does not work as filter in filterBam() function



## filter to a file
filter <- FilterRules(list(cigarM = function(x) {
  m <- str_extract_all(x$cigar,"\\d+M") 
  IDS <- str_extract_all(x$cigar,"\\d+[IDS]")
  nm <- as.numeric(x$nM)
  cat("m",length(m)," IDS",length(IDS)," nm",length(nm),"\n")
  # n.b. unlike in scanBam, x is a DataFrame with cigar, qwidth, nM, NH as columns
  # n.b. x$cigar is a vector of all cigars (in chunks) of all sequences
  # hence, m is also a vector, hence, use sapply
  # min 25M & max k percent D/I/S cigar letters
  (sapply(m,function(x) sum(as.numeric(gsub("M","",x)))) >= minM) & 
    (sapply(IDS,function(x) sum(as.numeric(gsub("(\\d+)\\w+","\\1",x)))) + nm) <= x$qwidth*maxNonMPercent/100
  #sapply(m,function(x) sum(as.numeric(gsub("M","",x)))) >= minM   # min 25M
  # min 25M & nonM < x%: 
  # sapply(m,function(x) sum(as.numeric(gsub("M","",x)))) >= minM & 
  #   sapply(nonM,function(x) sum(as.numeric(gsub("(\\d+)\\w+","\\1",x)))) <= x$qwidth*maxNonMPercent/100
}))

filter1bam <- function(bamFile) {
  outFile <- gsub("\\.bam",paste0(suffix,".bam"),bamFile)
  if(file.exists(outFile)) {
    cat(outFile,"already exists, not overwriting\n")
  } else {
    dest <- filterBam(bamFile, outFile, param=param, filter=filter)
  }
}
# outFile <- gsub("\\.bam",paste0(suffix,".bam"),bamFiles[1])
# dest <- filterBam(bamFiles[1], outFile, param=param, filter=filter)

mclapply(X = bamFiles,function(x) filter1bam(x), mc.cores = nbCPUs)
# 27min, 15cpus, 96 Samples, matches + mismatches + InDels + soft-cl.
