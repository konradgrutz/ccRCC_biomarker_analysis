# 
# call fastq_screen for several samples
# 
# limitations: 
#  so far only single end seq. supported
#  only using bowtie2
#  only -nohits option used


args <- commandArgs(trailingOnly = T)
if(length(args) <4) {
  cat("ERROR: at least 4 parameters needed: 
  out dir (e.g. trimmed_uncontaminated/run-01/), 
  rna_screen configuration file (e.g. fastq_screen.conf.ecoli_ncRNA),
  number of threads to use,
   single end sequenced fastq.gz files\n")
  q("no")
}

outDir <- args[1]
configFile <- args[2]
nbCpus <- as.numeric(args[3])
infiles <- args[4:length(args)]

i <- 0
for (infile in infiles) {
  i <- i+1
  cat(i,"of",length(infiles)," ",round(i/length(infiles)*100),"%\n")
  cmd <- paste0("fastq_screen -conf ",configFile," -threads ",nbCpus," -nohits -aligner bowtie2 -force -outdir ",
                outDir, " ",infile)
  cat(cmd,"\n")
  system(cmd)
}

cat("\nDONE\n")

