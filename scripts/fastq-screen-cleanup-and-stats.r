# 
# after calling fastq_screen_parallel.r
# apply this script to clean and make statistics over all samples and runs
#
# clean: for each given run folder remove *.tagged.fastq.gz files, move *.tagged_filter.fastq.gz to *.fastq.gz, remove *png and *html
# statistics: for each run and sample, read in sample-id_screen.txt; print min/max/mean of % of reads kept after filtering
#

args <- commandArgs(trailingOnly = T)
if(length(args) <2) {
  cat("ERROR: at least 2 parameters needed: 
      outfile base name (e.g. fastq_screen_overview),
      run folders of trimmed_fastq_screen (e.g. trimmed_uncontaminated/run-01 trimmed_uncontaminated/run-02)\n")
  q("no")
}

fastqScreenOutBase <- args[1]
outDirs <- args[2:length(args)]
cat("outDirs (",length(outDirs),")",outDirs)
outDir <- outDirs[1]
percentagesNoHitAllRuns <- list()
allPercentages <- NULL # percentage of each run, each sample
htmlOutText <- "<html>\n<body width=\"1000\">\n"

for (outDir in outDirs ) {
  cat("directory:",outDir,"\n")
  htmlOutText <- c(htmlOutText,"<h1>",unlist(strsplit(outDir,"\\/"))[2],"</h1>\n")

  # remove tagged.fastq.gz files, these contain all reads:
#  files <- list.files(path = outDir,pattern = ".+tagged\\..*gz", full.names = T)
#  if(length(files)>0) {
#    cat("removing",files,"\n")
#    file.remove(files)
#  }
  
  # rename tagged_filter.fastq.gz to original file name:
  currentFileNames <- list.files(path = outDir,pattern = ".+tagged_filter\\..*gz", full.names = T)
  desiredFileNames <- gsub("\\.tagged_filter","",currentFileNames)
  if(length(currentFileNames)>0) {
    if(all(currentFileNames != desiredFileNames)) {
      file.rename(from = currentFileNames, to = desiredFileNames)
    }
  }
  
  # remove html
  files <- list.files(path = outDir,pattern=".*_screen.html", full.names = T)
  if(length(files)>0) {
    cat("removing",files,"\n")
    file.remove(files)
  }
  
  # list png files to build a huge html with them:
  files <- list.files(path = outDir,pattern=".*_screen.png", full.names = T)
  if(length(files)>0) {
    htmlOutText <- c(htmlOutText,paste0("<img width=\"45%\" src=\"",files,"\">\n"))
  }
  
  
  # read in and make stats:
  files <- list.files(path=outDir, pattern=".*_screen.txt", full.names = T)
  cat(files,"\n")
  percentages <- NULL
  f <- files[1]
  for(f in files) {
    dat <- read.csv(file = f,header = F,sep = "\t",skip=1)
    percentages <- c(percentages,as.numeric(gsub("%Hit_no_genomes: ","",dat$V1[nrow(dat)])))
  }
  names(percentages) <- sapply(files,function(f) {x <- unlist(strsplit(f,"\\/")); x <- x[length(x)]; gsub("_screen.txt","",x)})
  
  percentagesNoHitAllRuns <- c(percentagesNoHitAllRuns,list(percentages))
  names(percentagesNoHitAllRuns)[length(percentagesNoHitAllRuns)] <- outDir

  if(!is.null(allPercentages)) {
    newnames <- setdiff(names(percentages),rownames(allPercentages))
    if(length(newnames)>0) { # add new samples in rows
      # 
      allPercentages[nrow(allPercentages)+length(newnames),] <- NA
      rownames(allPercentages)[(nrow(allPercentages)-length(newnames)+1):nrow(allPercentages)] <- newnames
    } 
    allPercentages[names(percentages),outDir] <- percentages
  } else {
    allPercentages <- data.frame("one"=percentages)
    rownames(allPercentages) <- names(percentages)
    colnames(allPercentages)[1] <- outDir
  }
}
htmlOutText <- c(htmlOutText,"</body></html>\n")
writeLines(text = htmlOutText,con = paste0(fastqScreenOutBase,".html"))

stats <- sapply(percentagesNoHitAllRuns, function(p) c(min(p, na.rm = T),max(p, na.rm = T),mean(p, na.rm = T)))
rownames(stats) <- c("min","max","average")
colnames(stats) <- gsub("trimmed_uncontaminated/","",colnames(stats))
stats <- round(stats,1)
stats <- as.data.frame(stats)
stats <- cbind(stat=rownames(stats),stats)
cat("percent of reads without a hit/kept:\n")
write.table(x = stats,file = paste0(fastqScreenOutBase,".csv"),quote = F,sep = "\t",row.names = F)

# add all percentages to the same file:
allPercentages <- cbind("sample"=rownames(allPercentages),allPercentages)
colnames(allPercentages) <- gsub("trimmed_uncontaminated/","",colnames(allPercentages))

write.table(x = allPercentages,file = paste0(fastqScreenOutBase,".csv"),quote = F,sep = "\t",row.names = F, append=T)

cat("convert html to pdf:",paste0(fastqScreenOutBase,".html"),"\n")
cat("then remove png files:\n",paste0("rm ",outDirs,"/*.png\n"))
