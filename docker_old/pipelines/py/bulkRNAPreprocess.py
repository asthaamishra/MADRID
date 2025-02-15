#!/usr/bin/python3
from bioservices import BioDBNet
import pandas as pd
from project import configs
import re
import os, time, sys
import getopt
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
#from rpy2.robjects import r, pandas2ri
#import rpy2.robjects as ro
#from rpy2.robjects.conversion import localconverter
#from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

#pandas2ri.activate()

#limma = importr("limma")
tidyverse = importr("tidyverse")
#edgeR = importr("edgeR")
#genefilter = importr("genefilter")
#biomaRt = importr("biomaRt")
#sjmisc = importr("sjmisc")

# automatically convert ryp2 dataframe to Pandas dataframe
string="""
library(tidyverse)

organizeFiles <- function(data_dir, technique) {
  SampMetrics <- list()
  countFiles <- list()
  fragmentFiles <- list()
  N_samples <- list()
  Names <- list()
  count_dir_i <- paste(c(data_dir, "/geneCounts"),collapse="")
  tab_dir <- list.dirs(path = count_dir_i, full.names = TRUE, recursive = FALSE)
  frag_dir_i <- paste(c(data_dir, "/fragLengths"),collapse="")
  frag_dir <- list.dirs(path = frag_dir_i, full.names = TRUE, recursive = FALSE)
  for (j in 1:length(tab_dir) ) {
    s <- tab_dir[j]
    f <- frag_dir[j]
    sname <- unlist(strsplit(s,"/"))
    sname <- sname[length(sname)]
    entry <- list()
    if ( technique=="zFPKM" ) {
      fragGlob <- paste(c(f, "/*.frag"), collapse="")
      fragFiles <- Sys.glob(fragGlob)
      entry[["FragmentFiles"]] <- fragFiles
    }
    else {
      entry[["FragmentFiles"]] <- NA
    }
    fileGlob <- paste(c(s, "/*.tab"), collapse="")
    cntFiles <- Sys.glob(fileGlob)
    n_samples <- length(cntFiles)
    sample_names <- rep(0,n_samples)
    for (i in 1:n_samples) {
      samp_file <- str_match(cntFiles[i], "geneCounts/\\\\s*(.*?)\\\\s*.tab")[,2]
      sample_names[i] <- unlist(strsplit(samp_file,"/"))[2]
    }
    entry[["CountFiles"]] <- cntFiles
    entry[["NumSamples"]] <- n_samples
    entry[["SampleNames"]] <- sample_names
    SampMetrics[[sname]] <- entry
  }
  return(SampMetrics)
}

prepSampCnts <- function(cntFile,files) {
  if ( grepl("R\\\\d+r1", cntFile, ignore.case=FALSE) ) {
    basename <- unlist(strsplit(cntFile, "r1"))[1]
    basename <- unlist(strsplit(basename, "/"))
    basename <- basename[length(basename)]
    search_str <- paste(c(basename, "r\\\\d+"),collapse="")
    run_files <- files[grepl(search_str, files)]
    samp_count <- NULL
    for ( f in run_files) {
      run_count <- read.delim(f)
      run_len <- length(run_count[,1])
      run_count <- run_count[4:run_len,]
      run_count <- na.omit(run_count)
      genes <- run_count[,1]
      cl <- which.max(colSums(run_count[,2:4])) + 1
      run_count <- data.frame(cbind(run_count[,1], run_count[,cl]))
      run_count[,1] <- genes # force rewrite genes bc rpy2 does not handle previous step properly
      colnames(run_count) <- c("genes", "counts")
      if ( is.null(samp_count) ) {
        samp_count <- run_count
      }
      else {
        samp_count <- merge(samp_count, run_count, by="genes", all=TRUE)
        samp_count[is.na(samp_count)] <- 0
      }
    }
    genes <- samp_count["genes"]
    samp_count["genes"] <- NULL
    samp_count <- sapply(samp_count, as.numeric)
    samp_sum <- rowSums(samp_count)
    samp_count <- data.frame(cbind(genes, samp_sum))
    samp_count[,1] <- genes # rewrite bc rpy2 is weird
    colnames(samp_count)[2] <- "counts"
    return(samp_count)
  }
  else if ( grepl("R\\\\d+r", cntFile)==TRUE ) {
    return("skip")
  }
  else {
    samp_count <- read.delim(cntFile)
    samp_len <- length(samp_count[,1])
    samp_count <- samp_count[4:samp_len,]
    samp_count <- na.omit(samp_count)
    genes <- samp_count[,1]
    cl <- which.max(colSums(samp_count[,2:4])) + 1
    samp_count <- data.frame(cbind(samp_count[,1], samp_count[,cl]))
    samp_count[,1] <- genes # force rewrite genes bc rpy2 does not handle previous step properly
    colnames(samp_count) <- c("genes", "counts")
    return(samp_count)
  }
}

createCountMatrix <- function(files, sample_names, n_samples,
                              fragFiles, insertSizes, technique) {
                              
  if ( technique=="zFPKM" ) {
    if ( grepl("r1",fragFiles[1],ignore.case=FALSE ) ) {
      groupSize <- c(fragFiles[1])
    }
    else {
      groupSize <- c()
      re_insertSizes <- c(insertSizes[1])
    }
  }
  i_adjust <- 0
  counts <- prepSampCnts(files[1])
  colnames(counts)[2] <- sample_names[1]
  for ( i in 2:n_samples) {
    new_cnts <- prepSampCnts(files[i], files)
    options(warn=-1)
    if ( new_cnts=="skip" ) {
      i_adjust <- i_adjust+1
      if ( technique=="zFPKM" ) {
        groupSize <- c(groupSize, insertSizes[i])
      }
      next
    }
    counts <- merge(counts, new_cnts, by="genes", all=TRUE)
    counts[is.na(counts)] <- 0
    if ( grepl("R\\\\d+r1", sample_names[i], ignore.case=FALSE) ) {
      samp_name <- unlist(strsplit(sample_names[i], "r"))[1]
      if ( technique=="zFPKM" ) {
        if ( length(groupSize) > 1 ) {
          re_insertSizes[i-i_adjust] <- mean(groupSize)
          groupSize <- c(insertSizes[i])
        }
        else {
          re_insertSizes[i-i_adjust] <- insertSizes[i]
        }
      }
    }
    else {
      samp_name <- sample_names[i]
      if ( technique=="zFPKM" ) {
        if ( length(groupSize) > 1 ) {
          re_insertSizes[i_adjust] <- mean(groupSize)
          groupSize <- c()
        }
        else {
          re_insertSizes[i-i_adjust] <- insertSizes[i]
        }
      }
    }
    colnames(counts)[i+1-i_adjust] <- samp_name
  }
  if ( technique=="zFPKM" ) {
    if ( length(groupSize) > 1 ) {
      re_insertSizes[i-i_adjust] <- mean(groupSize)
      groupSize <- c()
    }
  }
  if ( technique=="zFPKM" ) {
    res <-list("counts"=counts, "insertSizes"=re_insertSizes)
  }
  else {
    res <- counts
  }
  return(res)
}

genCountMatrix_main <- function(data_dir, out_dir, technique="quantile") {
  
  SampMetrics <- organizeFiles(data_dir, technique)
  
  
  for ( i in 1:length(SampMetrics) ) {
    if ( technique=="zFPKM" ) {
      j<-0
      insertSizes<-rep(0, SampMetrics[[i]][["NumSamples"]])
      for ( file in SampMetrics[[i]][["FragmentFiles"]] ) {
        j<-j+1
        lines <- readLines(file, n= 50)
        read_idx <- grep("METRICS CLASS", lines)
        size <- as.numeric(unlist(strsplit(lines[read_idx+2],"\t"))[6])
        insertSizes[j] <- as.numeric(unlist(strsplit(lines[read_idx+2],"\t"))[6])
      }
      SampMetrics[[i]][["InsertSizes"]] <- insertSizes
    }
    else {
      SampMetrics[[i]][["InsertSizes"]] <- NA
    }
  }
  for ( i in 1:length(SampMetrics) ) {
    res <- createCountMatrix(SampMetrics[[i]][["CountFiles"]],
                             SampMetrics[[i]][["SampleNames"]],
                             SampMetrics[[i]][["NumSamples"]],
                             SampMetrics[[i]][["FragmentFiles"]],
                             SampMetrics[[i]][["InsertSizes"]],
                             technique)
    if ( technique=="zFPKM" ) {
      SampMetrics[[i]][["CountMatrix"]] <- res$counts
      SampMetrics[[i]][["InsertSizes"]] <- res$insertSizes
    }
    else {
      SampMetrics[[i]][["CountMatrix"]] <- res
    }
    SampMetrics[[i]][["NumSamples"]] <- ncol(SampMetrics[[i]][["CountMatrix"]])
    if ( i == 1 ) {
      #full_count_matrix <- data.frame(SampMetrics[[i]][["CountMatrix"]])
      full_count_matrix <- SampMetrics[[i]][["CountMatrix"]]
    } else {
      #add_mat <-data.frame(SampMetrics[[i]][["CountMatrix"]])
      add_mat <- SampMetrics[[i]][["CountMatrix"]]
      full_count_matrix <- merge(full_count_matrix, add_mat,
                                 by="genes", all=TRUE)
      #row.names(full_count_matrix) <- full_count_matrix$Row.names
      #full_count_matrix["Row.names"] <- NULL
    }
  }
  file_split <- unlist(strsplit(data_dir, "/"))
  file_name <- paste(c(
    out_dir, "/", "BulkRNAseqDataMatrix_", file_split[length(file_split)], ".csv"),collapse="") 
  write.csv(full_count_matrix, file_name, row.names=FALSE)
  cat("Count Matrix written at ", file_name, "\n")
}
"""
genCountMatrixio = SignatureTranslatedAnonymousPackage(string, "genCountMatrixio")

def fetch_gene_info(input_values, input_db='Ensembl Gene ID',
                    output_db=['Gene Symbol','Gene ID','Chromosomal Location'],
                    delay=15):
 
    s = BioDBNet()   
    # input_db = 'Agilent ID'
    # input_values = df_results.ProbeName.tolist()

    df_maps = pd.DataFrame([],columns=output_db)
    df_maps.index.name=input_db
    i = 0
    # for i in range(0,len(input_values),500):
    while i < len(input_values):
        print('retrieve {}:{}'.format(i,min(i+500,len(input_values))))
        df_test = s.db2db(input_db, output_db, input_values[i:min(i+500,len(input_values))], 9606)
        if isinstance(df_test, pd.DataFrame):
            df_maps = pd.concat([df_maps, df_test], sort=False)
        elif df_test == '414':
            print("bioDBnet busy, try again in {} seconds".format(delay))
            time.sleep(delay)
            continue
        i += 500
    return df_maps

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hn:c:f:t:",
                                   ["tissue_name=",
                                    "create_counts_matrix=",
                                    "gene_format=",
                                    "technique="])
    except getopt.GetoptError:
        print('python3 bulkRNAPreprocess.py -n <tissue_name> -c <create_counts_matrix> -f <gene_format> -t <technique>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 bulkRNAPreprocess.py -n <tissue_name> -c <create_counts_matrix> -f <gene_format> -t <technique>')
            sys.exit()
        elif opt in ("-n", "--tissue_name"):
            tissue_name = arg
        elif opt in ("-c", "--create_counts_matrix"):
            make_matrix = arg
        elif opt in ("-f", "--gene_format"):
            gene_format = arg  
        elif opt in ("-t", "--technique"):
            technique = arg
    input_dir = os.path.join(configs.rootdir, 'data', 'bulkData', tissue_name)
    output_dir = os.path.join(configs.rootdir, 'data')
    print('Input directory is "{}"'.format(input_dir))
    print('Output directory is "{}"'.format(output_dir))
    print('Active gene determination technique is "{}"'.format(technique))
    if make_matrix:
        print("Creating Counts Matrix")
        genCountMatrixio.genCountMatrix_main(input_dir, output_dir, technique)
    geneCountFile = os.path.join(output_dir, ("BulkRNAseqDataMatrix_"+tissue_name+".csv"))
    print('Fetching gene info using genes in "{}"'.format(geneCountFile))
    genes = pd.read_csv(geneCountFile)['genes'].to_list()
    output_db=['Ensembl Gene ID', 'Gene Symbol', 'Gene ID', 'Chromosomal Location']               
    if gene_format.upper()=="ENSEMBL":
         form = "Ensembl Gene ID"                     
    elif gene_format.upper()=="ENTREZ":
         form = "Gene ID"
    elif gene_format.upper()=="SYMBOL":
          form = "Gene Symbol"                       
    output_db.remove(form)  
    gene_info = fetch_gene_info(genes, input_db=form, output_db=output_db)
    gene_info['start_position'] = gene_info['Chromosomal Location'].str.extract("chr_start: (\d+)")
    gene_info['end_position'] = gene_info['Chromosomal Location'].str.extract("chr_end: (\d+)")
    gene_info.index.rename("ensembl_gene_id", inplace=True)
    gene_info.rename(columns={"Gene Symbol": "hgnc_symbol", "Gene ID": "entrezgene_id"}, inplace=True)
    gene_info.drop(['Chromosomal Location'], axis=1, inplace=True)
    gene_info_file = os.path.join(output_dir, ("GeneInfo_"+tissue_name+".csv"))
    gene_info.to_csv(gene_info_file)
    print('Gene Info file written at "{}"'.format(gene_info_file))


if __name__ == "__main__":
    print(sys.argv)
    main(sys.argv[1:])
