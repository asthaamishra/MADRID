# This file will cluster transcriptomic and proteomic sources
# It will perform the following functions
#   1. reads the zFPKMs (or TPMs or CPMs) for each replicate
#               *binarizes them if given that arg (see RNAseq.R using the kOverA function)
#   2. clusters all replicates for all batches/studies across all contexts/cell types within one data source at a time
#   3. clusters all replicates for all batches/studies across all data sources within one cell type at a time
#   4. clusters all batches/studies across all contexts/cell types within one data source at a time
#   5. clusters all batches/studies across all data sources within one cell type at a time (edited)

library(stringr)
library(tools)
library(uwot)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)


get_study_value <- function(file_path) {
    # This function will get the S## value from the file path
    
    # Find "S##.csv", but use a lookahead to avoid "taking" the ".csv" portion
    match <- stringr::str_extract(string = toupper(file_path), pattern = "S\\d{1,2}(?=\\.CSV)")

    # If match is not NA, uppercase the value
    if (!is.na(match))
        match <- toupper(match)

    return(match)
}

get_replicate_files <- function(results_directory, context_names, source_type, use_trna, use_mrna) {
    # This function will get the file paths of zFPKMs for each replicate

    # Create a list to hold all cell types
    all_context_files <- list()

    # Iterate through each context/cell type given
     for (context_name in context_names) {
        lower_context_name <- tolower(context_name)
        source_type <- tolower(source_type)
        
        current_context_files <- list.files(file.path(results_directory, lower_context_name), full.names = TRUE, recursive = TRUE)
        context_files <- c()
        
        for (file in current_context_files) {
            file_name <- tolower(file)
            
            # Check the current file meets our criteria
            if (
              tools::file_ext(file_name) == "csv" &&                      # Ensure the current file is a CSV file
              grepl(lower_context_name, file_name) &&  # Test if the current file is part of the current context (i.e., naiveB, immNK)
              grepl(source_type, file_name) &&   # Test if the current file has source type (zFPKM, TPM, CPM)
              # Create a "group" that works if either trna or mrna is TRUE
              (
                (use_trna && grepl("total", file_name)) ||         # Test if the current file is a total-rna file
                (use_mrna && grepl("mrna", file_name))             # Test if the current file is an mRNA (polyA) file
              )
            ) {
                context_files <- append(context_files, file)
            }
        }
        
        # Only append new list if context_files has at least one item
        if (length(context_files) > 0)
            all_context_files[[context_name]] <- context_files
    }
    
    # Return list if it has at least one item, otherwise return "NA"
    if (length(all_context_files) > 0) {
        return(all_context_files)
    } else {
        return(NA)
    }
}

read_matrix_values <- function(context_files) {
    # This function is responsible for reading in the matrix values found within the replicate files
    # It takes the list of replicate files and returns a list of lists of matrix values
    replicate_dataframes <- list()
    context_names <- names(context_files)

    # Iterate through each context/cell type
    for (context in context_names) {

        # Define the index of replicate_dataframes we will be appending to
        index <- 1

        # Collect our current context data
        current_context_files <- context_files[[context]]
        context_dataframe <- c()
        for (file in current_context_files) {
            dataframe <- read.csv(file, header = TRUE)

            # Get the S## value from the file path
            study <- get_study_value(file)
            context_dataframe[[study]] <- dataframe

            index <- index + 1
        }

        # Only append new list if context_dataframe has at least one item
        if (length(context_dataframe) > 0) {
            replicate_dataframes[[context]] <- context_dataframe
        }
    }

    # Return list if it has at least one item, otherwise return "NA"
    if (length(replicate_dataframes) > 0) {
        return(replicate_dataframes)
    } else {
        return(NA)
    }
}


binarize_matrix_values <- function(context_dataframes, default_bin) {
    # study_dataframes: a list of lists of dataframes
    # study_dataframes[[1]]: naiveB
    # study_dataframes[[1]][["S1"]]: S1 dataframe
    # study_dataframes[[1]][["S2"]]: S2 dataframe
    # default_bin: The value at which 0-bin ("off") items turn to 1-bin items ("on")


    binarized_context_dataframes <- list()
    for (context in context_dataframes) {

        # Set the index to append to binarized_dataframes
        index <- 1

        # Get the "naiveB" portion from "naiveB_S1R1"
        # Split the string at the underscore and get the first element
        # Access the first item from the list, and take the first table in the list
        # This will give us the column names of the data frame (ENTRZ_GENE_ID, immNK_S1R1, immNKS1R2, etc.)
        col_names <- colnames(context[1][[1]])
        context_name <- stringr::str_split(col_names[[2]], "_", simplify = TRUE)[[1]]
        context_dataframe <- c()
        for (study in context) {
            # Get the S## from the study name
            # Split the column name at "_", and take the "S##R##" portion
            # From here, use regex to get the S## value
            study_name <- stringr::str_split(colnames(study)[2], "_", simplify = TRUE)[[2]]
            study_name <- stringr::str_extract(string = study_name, pattern = "S\\d{1,2}")

            # Get the number of columns in the dataframe
            n_cols <- ncol(study)

            # Values less than or equal to -3 will be set to 0
            # Values greater than -3 will be set to 1
            # Start on column 2, since column 1 is the ENTREZ_GENE_ID column
            binarized_study <- ifelse(study[,2:n_cols] <= -3, 0, 1)
            context_dataframe[[study_name]] <- binarized_study
        }

        binarized_context_dataframes[[context_name]] <- context_dataframe

    }
    return(binarized_context_dataframes)
}


cluster_dataframe <- function(master_dataframe, cluster_neighbors = 30) {
    if (cluster_neighbors > nrow(master_dataframe)) {
        print("Number of neighbors is greater than number of rows in the dataframe. Setting number of neighbors to number of rows in thedataframe.")
        cluster_neighbors <- nrow(master_dataframe)
    } else if (cluster_neighbors < 2) {
        print("Number of neighbors is less than 2. Setting number of neighbors to 2.")
        cluster_neighbors <- 2
    }

    cluster_umap <- uwot::umap(
      X = master_dataframe,
      n_threads = 12,
      n_neighbors = cluster_neighbors,
      min_dist = 0.01
    )
    cluster_umap_df <- data.frame(cluster_umap)

    # Rename X0 and X1 to x and y
    # format: rename(NEW = OLD)
    cluster_umap_df <- cluster_umap_df %>% dplyr::rename(x = X1, y = X2)

    # Create a new column called "context_study", equal to the row names. Then remove row names
    cluster_umap_df$context_study <- row.names(cluster_umap_df)
    row.names(cluster_umap_df) <- NULL

    # Separate the context_study (naiveB_S1R1) into context (naiveB) and study (S1R1)
    cluster_umap_df <- cluster_umap_df %>% tidyr::separate(data = ., col = context_study, into = c("context", "study"), sep = "_")
    return(cluster_umap_df)
}


plot_cluster_dataframe <- function(clustered_dataframe) {
    file_save_path <- "figures/cluster_umap.pdf"
    pdf(file_save_path)
    p <- ggplot2::ggplot(clustered_dataframe, ggplot2::aes(x=x, y=y, label=study, color=context)) +
      ggplot2::geom_point(alpha=0.7) +  # Plot points
      ggrepel::geom_text_repel(max.overlaps = Inf, show.legend = FALSE, max.time = Inf, max.iter = 300000) +  # Add text labels
      ggplot2::labs(x="Dim 1", y="Dim 2")  # Add axis labels

    print(p)  # Save the plot to the file name
    dev.off()

    return (file_save_path)
}


create_all_replicates_all_studies_all_context_one_source_matrix <- function(binarized_dataframes) {
    # This function is responsible for clustering all replicates across all studies within a single data source for all contexts
    # For example, clustering all naiveB and immNK studies within RNA-seq data

    master_dataframe <- data.frame()
    row_names <- c()

    for (context_num in seq_along(binarized_dataframes)) {
        current_context_dataframe <- binarized_dataframes[[context_num]]
        for (study_num in seq_along(current_context_dataframe)) {

            current_study <- current_context_dataframe[[study_num]]

            context_name <- names(binarized_dataframes)[context_num]  # Cell type (naiveB, immNK, etc.)
            study_name <- names(current_context_dataframe)[study_num]     # S## value (S1, S2, etc.)
            replicate_name <- names(current_study)[1]                     # R## value (R1, R2, etc.)

            # '~' is an anonymous function (i.e., execute max() on the incoming data)
            # From: https://stackoverflow.com/a/54834959
            # Collect all values for each entrez gene ID, and then take the max() of those values
            # Taken in part from: https://stackoverflow.com/a/63960722
            current_study <- current_study %>%
              dplyr::group_by(ENTREZ_GENE_ID) %>%  # Group by entrez gene ID
              dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = list(name = ~max(.))))

            if (context_num == 1) {
                master_dataframe <- current_study
            } else {
                master_dataframe <- master_dataframe %>% dplyr::left_join(current_study, by="ENTREZ_GENE_ID")
            }
        }
    }

    master_dataframe_transpose <- t(master_dataframe[-1])
    colnames(master_dataframe_transpose) <- master_dataframe$ENTREZ_GENE_ID

    return (master_dataframe_transpose)
}

cluster_all_replicates_all_studies_single_context_all_source <- function() {
    # This function will cluster all replicates across all studies within a single context
    # For example, clustering all immNK studies within all data sources provided
}

cluster_all_studies_all_contexts_single_source <- function() {
    # This function will cluster all studies across all contexts within a single data source
    # For example, cluster all immNK and naiveB within all RNAseq data
}

cluster_all_studies_one_context_all_source <- function () {
    # This function will cluster all studies within a single context for all given data sources
    # For example, clustering all immNK studies across all data sources given
}

main <- function(
  results_directory,
  context_names,
  source_type,
  use_trna,
  use_mrna,
  binarize_data,
  default_bin
) {
    
    context_files <- get_replicate_files(results_directory = results_directory, context_names = context_names,  source_type = source_type, use_trna = use_trna, use_mrna = use_mrna)
    context_dataframes <- read_matrix_values(context_files = context_files)

    if (binarize_data == TRUE)
      context_dataframes <- binarize_matrix_values(context_dataframes = context_dataframes, default_bin = default_bin)

    print("Creating matrix")
    all_replicate_all_studies_all_context_matrix <- create_all_replicates_all_studies_all_context_one_source_matrix(binarized_dataframes = context_dataframes)

    print("Clustering")
    all_replicate_all_studies_all_context_clustered <- cluster_dataframe(master_dataframe = all_replicate_all_studies_all_context_matrix)

    print("Plotting")
    all_replicate_all_studies_all_context_clustered_plot_path <- plot_cluster_dataframe(clustered_dataframe = all_replicate_all_studies_all_context_clustered)

    print("DONE")
}


results_directory <- "/Users/joshl/docker/madrid/local_files/results"
context_names <- list("immNK", "naiveB")
source_type <- "zFPKM"
use_trna <- TRUE
use_mrna <- TRUE
binarize_data <- FALSE
default_bin <- -3
main(results_directory = results_directory, context_names = context_names, source_type = source_type, use_trna = use_trna, use_mrna = use_mrna, binarize_data = binarize_data, default_bin = default_bin)
