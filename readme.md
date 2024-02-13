# Data analysis pipeline

for the publication

### Identification and validation of novel snoRNA-based biomarkers for clear cell renal cell carcinoma from urine-derived extracellular vesicles

Konrad Grützmann, Karsten Salomo, Alexander Krüger, Andrea Lohse-Fischer, Kati Erdmann, Michael Seifert, Gustavo Baretton, Daniela Aust, Doreen William, Evelin Schröck, Christian Thomas, Susanne Füssel


## directory structure

    - annotation/
      - gene annotations (Homo_sapiens.GRCh37.75-chrRename-noHaplo.RData)
      - chromosome lengths (hg19.chrom.sizes)
      - sample annotation (annotation-63.csv|xls, annotation-78.csv|xls)
		- patient obesity and hypertension risk factors (annotation-obesity-hypertension.xlsx)
      - fungi-base-urls.txt for fastq_screen
      - further gene annotation of the signif regions (signif-regions-hand-curated-annotation.txt)
      - exon intron structure of SNHG1 for plot (UCSC-genes-hg19-SNHG1-structure.gtf)
    - conda /
      - the required conda environment will be installed there (installation see below)
    - data/
      - coverage/ for genomic region coverage for each sample
      - genomes/ to download and index genomes for fastq_screen
      - PCR/  data we measured to validate the gene candidates
      - star/ STAR alignments for each sample
      - trimmed/ where the trimmed reads will be
      - trimmed_uncontaminated/ trimmed reads that are uncontaminated after fastq_screen
    - FiguresTables/
      - Figures and Tables of the paper produced with the pipeline will be stored there
    - scripts
       - all R scripts and jupyter notebooks of the analysis explained below 

    - download all these data from zenodo: https://doi.org/10.5281/zenodo.10654767

## Setting up conda environments


- use this conda call to setup the full environment

    ```
    conda create --prefix conda -c conda-forge -c bioconda r-base=4 r-irkernel notebook \
      r-circlize bioconductor-complexheatmap r-ggrepel bioconductor-deseq2 r-ggplot2 \
      r-gplots r-scales r-gridextra r-proc r-xlsx r-readxl bioconductor-rsamtools \
      r-reshape2 r-stringr parallel bowtie2 fastq-screen samtools STAR bedtools
    ```

- if you do not have the fastq reads (we cannot provide them for our cohort for privacy reasons), you will start the analyses with clustering (see below). This requires only a smaller conda environment without STAR, fastq_screen, bowtie, samtools, bedtools:

    ```
    conda create --prefix conda -c conda-forge -c bioconda r-base=4 r-irkernel notebook \
      r-circlize bioconductor-complexheatmap r-ggrepel bioconductor-deseq2 r-ggplot2 \
      r-gplots r-scales r-gridextra r-proc r-xlsx r-readxl bioconductor-rsamtools \
      r-reshape2 r-stringr parallel 
    ```
   
- then activate the conda environment and start working

    ```
    conda activate conda/
    ```
- for the jupyter notebooks start jupyter locally:

    ```
    jupyter notebook
    ```
    
- or work remotely. On the server:

    ```
    jupyter notebook --port=9000 --no-browser
    ```
- on your local machine:

    ```
    ssh -N -L 9000:localhost:9000 user@remote-server.domain.com
    ```
- then open a internet browser with 
    ```
    localhost:9000/
    ```
    

## Full workflow

- For privacy reason, we cannot publish read data, but we give instructions on how we handled the fastq read data. This includes: trimming, fastq_screen, mapping with STAR, and filtering bam files by cigar string
- without the fastq or bam files you start at the point "Reduced workflow" with clustering found below

### read filtering

- trim reads after sequencing with cutadapt with the parameters -q 15,15 -l 50 -m 15 -u 3 -a AAAAAAAAAA

#### fastq_screen

- to filter reads for possible contamination
- search and download the genomes of the species that may show up in the urine samples
- NCBI assembly:
  - 328 archea
    query: Search all[filter] AND archaea[filter] AND "latest refseq"[filter] AND "complete genome"[filter] AND  ( all[filter] NOT "derived from surveillance project"[filter] AND all[filter] NOT anomalous[filter] )   Sort by: SORTORDER
  - 1529 representative bacteria 
    query: Search all[filter] AND bacteria[filter] AND "latest refseq"[filter] AND "complete genome"[filter] AND "representative genome"[filter] AND (all[filter] NOT "derived from surveillance project"[filter] AND all[filter] NOT anomalous[filter] )
- 280 fungal species from NCBI RefSeq, all urls are in annotation/fungi-base-urls.txt, use e.g. wget to download these
- from Ensembl:
  - barley (Hordeum_vulgare.IBSC_v)2, wheat (Triticum_aestivum.IWGSC), and from NCBI: corn (Zea mays), rice (Oryza sativa), sorghum, soybean (glycine_max), wine (vitis_vinifera), wheat (Triticum_aestivum.IWGSC.dna_rm.toplevel.fa.gz)
  - soybean: GCF_000004515.5_Glycine_max_v2.1_genomic.fna.gz, sorghum GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz, zea mays: GCF_000005005.2_B73_RefGen_v4_genomic.fna.gz, barley: Hordeum_vulgare.IBSC_v2.dna_rm.toplevel.fa.gz, wine: GCF_000003745.3_12X_genomic.fna.gz, rice: GCF_001433935.1_IRGSP-1.0_genomic.fna.gz
  
- then index all genomes using bowtie2 with standard parameters, e.g. 

    ```
    bowtie2-build --threads 8 -f GCF_000003745.3_12X_genomic.fna wine
    ```
    
- make a fastq_screen config file that includes all standard databases and the above downloaded and indexed genomes: fastq_screen_Adapter_Lambda_Mitoch_PhiX_rRNA_Vectors_my_genomes.conf
- call fastq_screen in parallel for the trimmed fastq.gz files:
  - at least 4 parameters needed: output directory, rna_screen configuration file, number of threads to use, single end sequenced fastq.gz files
  
    ```
    Rscript scripts/fastq_screen_parallel.r data/trimmed_uncontaminated/run-1/ \ 
        data/fastq_screen_Adapter_Lambda_Mitoch_PhiX_rRNA_Vectors_my_genomes.conf 6 data/trimmed/run-1/*gz
      ```

  - call the above fastq-screen-parallel.r for every NGS run separately

- cleanup and make statistics:

    ```
    Rscript scripts/fastq-screen-cleanup-and-stats.r fastq_screen_overview data/trimmed_uncontaminated/run-*
    ```
    - produces
      - fastq_screen_overview.html with all plots of all samples
      - fastq_screen_overview.csv with min/max/mean % of reads kept after filtering
    
- concatenate all fastq.gz of each run to one fastq.gz for each sample
  - important: don't concatenate *tagged.fastq.gz, only sampleName.fastq.gz, i.e. the filtered ones
  - all reads land in directory merged_runs_fastq/
  - in data/trimmed_uncontaminated/ call:
  
    ```
    ls run-*/*.fastq.gz|cut -d"/" -f2|cut -d"." -f1|cut -d"_" -f1|sort|uniq > sample-ids.csv
    for n in $(cat sample-ids.csv); do echo "gzip -cd run-*/$n*1.fastq.gz > \ 
      merged_runs_fastq/$n.fastq"; gzip -cd run-*/$n*1.fastq.gz > merged_runs_fastq/$n.fastq ; done
    ```
  - you should gzip the fastq files and remove the unpacked ones:

    ```
    for n in $(cat sample-ids.csv); do echo "gzip -c merged_runs_fastq/$n.fastq > \
      merged_runs_fastq/$n.fastq.gz "; gzip -c merged_runs_fastq/$n.fastq > merged_runs_fastq/$n.fastq.gz; done
    ```

### read mapping and filtering

- call STAR for each sample with default parameters, except these: --alignIntronMax 1 --outFilterMismatchNmax 1 --outFilterMatchNmin 16 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterMismatchNoverLmax 0.03
  - merge all splice junctions after first mapping round: cat */SJ.out.tab | cut -f1-4 | sort | uniq > SJ.all
  - then call STAR again adding the parameter --sjdbFileChrStartEnd SJ.all
  - this produces files Aligned.sortedByCoord.out.bam for each sample
- mark duplicated reads with sambamba:

    ```
    sambamba markdup -t 4 -l 9 -p Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.mkdup.bam
    ```
  
- reduce to uniquely mapped reads for each sample:

    ```
    samtools view -h sample/Aligned.sortedByCoord.out.mkdup.bam|grep --perl "NH:i:1[^\d]|^@" | \
      samtools view -h -b -o sample/Aligned.sortedByCoord.out.mkdup.uniqMap.bam - 
    ```

- index all bam files: 

    ```
    samtools index Aligned.sortedByCoord.out.mkdup.uniqMap.bam
    ```

- filter bam alignments for matched bases and soft-clipping
  - at least 25M in CIGAR string, at most 10% may have InDel or Soft-clip
  - 3 or more parameters needed:  number of CPUs, the least number of M in total in CIGAR for each read, max percent of InDels/soft-clipping, new suffix for bam files, bam files
  - call within each sample subdirectory:
  
    ```
    scripts/filter-bams-by-cigar-matches-soft-clip.r 4 25 10 -uniqMap25M10percIDS \
        Aligned.sortedByCoord.out.mkdup.uniqMap.bam
    ```
  
  - produces files Aligned.sortedByCoord.out-uniqMap25M10percIDS.bam

- generate chromosome regions/tiles:
  - scripts/make-chromosome-tiles.r.ipynb
  - yields annotation/chromosome-tiles.bed to determine genomic region coverage with bedtools 

### coverage variance

- scripts/calc-coverage-variance.r.ipynb
- slices the genome into 1 Mb sized tiles (coverage-variance-chromosome-tiles.bed)
- calculate the read coverage for each tile for each sample
- then, calculate the mean, s.d. and s.d./mean ("coverage variance") over all tiles for each sample
- produces coverage-variance.csv which is part of Table S1

    
### genomic region coverage

- use "bedtools coverage" to get the expression of each genomic region
- call for each sample 

    ```
    bedtools coverage -counts -a annotation/chromosome-tiles.bed \
        -b data/star/sample/Aligned.sortedByCoord.out-uniqMap25M10percIDS.bam > \
        data/coverage/sample-bedcov-uniqMap25M10percIDS.csv
    ```

### find expressed regions

- scripts/find-expressed-regions.r.ipynb
- uses the coverage files produced above (*-bedcov-uniqMap25M10percIDS.csv) to find regions/tiles that are expressed with at least 5 reads in at least 20% of the cohort (12 of 63 patients)
- generates read counts of the expressed regions: data/expressedTiles-featureCounts-25M10percIDS-min5reads20percent.RData
    
### nucleotide-wise RNA coverage

- for Figure 4 - nucleotide-wise coverage plot of SNHG1 and SNORD22, SNORD26

    ```
    samtools depth -H -g SECONDARY,QCFAIL,DUP -a \
        -b data/SNHG1-extended.bed data/star/*/Aligned.sortedByCoord.out-uniqMap25M10percIDS.bam \
        > data/sam-depth-25MpercIDS-uniqMap-wFlags-SNHG1-extended.txt
    ```
- for Suppl. Figures S5-S12

    ```
    samtools depth -H -g SECONDARY,QCFAIL,DUP -a \
        -b data/signif-regions-extended-for-covg-plot.bed \
        data/star/*/Aligned.sortedByCoord.out-uniqMap25M10percIDS.bam \
        > data/sam-depth-25MpercIDS-uniqMap-wFlags-all-signif-regions-extended.txt &
    ```
- see notes below to produce the figures from the data extracted here


## Reduced workflow

- you start from here if you do not have the fastq or bam files
- the full workflow continues here


### clustering

- scripts/clustering-and-PCA.r.ipynb
- principal component analysis and sample expression distance heatmap
- produces in the directory Figures-and-Tables/
  - Figure 2 - principle component analysis (PCA)
  - Suppl Figure S3 - PCA with higher principal components
  - Suppl Figure S4 - Heatmap with clustering of expression data
  - gene-expression.xls - expression table with total reads, counts, TPMs and rlog(TPMs)
  

### Differential region expression

- scripts/deseq2-withAdjustment.r.ipynb
- produces in the directory Figures-and-Tables/
  - Supplementary Table S4 with log-fold-changes, pvalues, etc., including counts, TPMs, rlog(TPMs)
  - Figure 3 - volcano plot


### Plot regions with nucleotide-wise read coverage

- uses the files produced above ("nucleotide-wise RNA coverage") by samtools

- scripts/plot-UCSC-gene-candidates.r.ipynb
    - produces Figure-4-nucl-wise-coverage.png in subfolder Figures-and-Tables/
    
- scripts/make-wig-UCSC-of-sign-regions.r.ipynb
    - script for Suppl. Figures S5-S12
    - produces wig files in subfolder Figures-and-Tables/ for upload into UCSC 
        - signif-region-SNORD99.wig, signif-region-tRNA-Met.wig, signif-region-RNU2-2P.wig, signif-region-SNORD22_26.wig, signif-region-mascRNA.wig, signif-region-SNORA50C.wig, signif-region-SNORA81.wig, signif-region-SNORD50B.wig
    - load these into UCSC: My Data -> Custom Tracks -> upload and submit
    - right click the UCSC wig track and select full instead of dense view
    - also produces Figures-and-Tables/differentially-expressed-tiles.bed for UCSC browser. These are the vertical light blue bars in the Figures. They must be added by hand: zoom into the given region, mark the region by drag and drop inside the upper part of UCSC browser (chromosome coordinates), select "Add Highlight" in the popup menu.
    - note, all coordinates are based on hg19/GRCh37 

  
### PCR validation and diagnostic models

- scripts/evaluate-PCR-predictions.r.ipynb
    - evaluates the qPCR validation measurements
    - produces
        - Table 4 with classification measurements (sensitivity, specificity, accuracy, ...) with and without adjustment for obesity and hypertension
        - Figure 5 with ROC curves
        - Table S5 that shows the performance when using 2-4 genes for classification and with or without adjustment for obesity and hypertension

### Crossvalidation of the diagnostic models

- scripts/model-crossvalidation.r.ipynb
    - does the fivefold repeated crossvalidation of the models with 1, 2, 3, 4 genes
    - produces Suppl. Table S6 with the mean and standard deviations of the model metrics

