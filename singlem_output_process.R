# The BiocManager amd pacman packages are R package management tools. 
# They make it more straightforward to install certain packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("pacman"))
  install.packages("pacman")

# Use the p_load function from the pacman package to check if packages are installed. 
# If not they are not installed, the function will download and and install them for you
pacman::p_load("tidyverse", "reshape2", "openxlsx", "data.table", "reshape2")

# SingleM is a tool that finds the abundances of discrete operational taxonomic units (OTUs) 
# directly from shotgun metagenome data. Specifically, SingleM searches the metagenomic reads 
# using a combination of BLAST and HMMs to find sequences encoding a set of conserved, 
# taxonomically annotated, single-copy marker genes (derived from the Genome Taxonomy Database (GTDB); 
# https://gtdb.ecogenomic.org). It then filters these sequences and clusters them into OTUs. 
# The taxonomy of each OTU is then assigned based on a consensus of the taxonomic assignment of all 
# reads that go into that OTU.


# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# Utility functions

#' Convert a dataframe to a matrix. First column becomes row names
#' @param mydataframe input dataframe
df2m <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}

#' Convert a matrix to a dataframe. Rows become first column.
#' @param column_name Name of the new column (defaults to "Row_variable")
m2df <- function(mymatrix, column_name = "Row_variable"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, column_name] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(column_name,cur_names)]
  return(mydf)
}

write_neat_xlsx <- function(dataframe_list.l, output_filename){
  
  # Create workbook
  wb <- createWorkbook()
  
  # Create styles
  header_style <- createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center"
  )
  
  body_style <- createStyle(
    halign = "center",
    valign = "center"
  )
  
  # Add data and styles to each sheet
  for (sheet_name in names(dataframe_list.l)) {
    df <- dataframe_list.l[[sheet_name]]
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, x = df, startRow = 1)
    
    # Apply header style
    addStyle(wb, sheet = sheet_name, style = header_style,
             rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
    
    # Apply body style
    addStyle(wb, sheet = sheet_name, style = body_style,
             rows = 2:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
    
    # Optional: auto-fit columns
    setColWidths(wb, sheet = sheet_name, cols = 1:ncol(df), widths = "auto")
  }
  saveWorkbook(wb, output_filename, overwrite = TRUE)
}


load_and_process_singleM_condensed_data <- function(singleM_output_file){
  # Load singleM data
  singlem_data.df <- read.table(singleM_output_file,sep = "\t", header = T)
  # singlem_data.df$sample <- gsub("_.*", "", singlem_data.df$sample)
  
  # Remove columns
  singlem_data.df$gene <- NULL
  singlem_data.df$sequence <- NULL
  
  # Remove blanks
  singlem_data.df <- singlem_data.df[singlem_data.df$taxonomy != "",]
  
  # Remove spaces around semicolon
  singlem_data.df$taxonomy <- gsub("; ", ";", singlem_data.df$taxonomy)
  
  # Separate taxonomies
  singlem_data.df <- separate(singlem_data.df, "taxonomy", into = c("Root", "Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"), remove =F, sep = ";")
  
  # Splitting taxa strings that are not specified at certain taxa levels will produce NA entries at those levels. 
  # NA entries should be changed to "Unassigned"
  singlem_data.df[is.na(singlem_data.df)] <- "Unassigned"
  
  # Construct taxonomy string
  singlem_data.df <- 
    singlem_data.df %>%
    dplyr::mutate(taxonomy_phylum = paste(Domain, Phylum, sep =";"),
                  taxonomy_class = paste(Domain, Phylum, Class, sep =";"),
                  taxonomy_order = paste(Domain, Phylum, Class, Order, sep =";"),
                  taxonomy_family = paste(Domain, Phylum, Class, Order, Family, sep =";"),
                  taxonomy_genus = paste(Domain, Phylum, Class, Order, Family, Genus, sep =";"),
                  taxonomy_species = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep =";")
    ) %>%
    as.data.frame()
  
  # Generate coverage table for each taxonomy level
  output.l <- list()
  taxonomy_string_levels.v <- c("taxonomy_phylum","taxonomy_class","taxonomy_order","taxonomy_family","taxonomy_genus","taxonomy_species")
  for (tax_string_level in taxonomy_string_levels.v){
    mydata.df <- singlem_data.df[,c("sample",tax_string_level, "coverage")]
    cov.df <- reshape2::dcast(mydata.df, 
                              formula = get(tax_string_level)~sample,
                              value.var = "coverage",
                              fill = 0,sum)
    cov.m <- df2m(cov.df)
    rel.m <- t(t(cov.m)/colSums(cov.m)) * 100
    rel.df <- m2df(rel.m,tax_string_level)
    # output.l[tax_string_level] <- list(tax_string_level = rel.df)
    output.l[[tax_string_level]] <- list("coverage" = cov.m, "abundances" = rel.m)
  }
  output.l
}

load_and_process_singleM_otu_data <- function(singleM_output_file){
  
  singlem_otu_data.df <- read.delim(singleM_output_file, sep = "\t", header = T)
  
  singlem_otu_data.df <- 
    singlem_otu_data.df %>% 
    # dplyr::mutate(sample = gsub("_.*","", sample)) %>%
    dplyr::select(-coverage) %>%
    tidyr::pivot_longer(cols = num_hits, values_to = "number_of_hits") %>%
    dplyr::select(-name) %>%
    dplyr::group_by(sequence) %>%
    dplyr::mutate(OTU_ID = paste0("OTU_", cur_group_id()), gene = gsub("\\.gpkg", "",gene)) %>% 
    dplyr::arrange(OTU_ID) %>%
    dplyr::group_by(sample, gene, taxonomy,sequence,OTU_ID) %>%
    dplyr::summarise_all(list(sum)) %>% 
    dplyr::mutate(taxonomy = gsub("Root;", "", gsub("; ", ";", taxonomy))) %>%
    tidyr::separate(., "taxonomy", into = c("Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"), remove =F, sep = ";") %>% 
    dplyr::mutate(across(c("Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"),~ifelse(is.na(.x),"Unassigned",.x))) %>%
    dplyr::mutate(taxonomy_phylum = paste(Domain, Phylum, sep =";"),
                  taxonomy_class = paste(Domain, Phylum, Class, sep =";"),
                  taxonomy_order = paste(Domain, Phylum, Class, Order, sep =";"),
                  taxonomy_family = paste(Domain, Phylum, Class, Order, Family, sep =";"),
                  taxonomy_genus = paste(Domain, Phylum, Class, Order, Family, Genus, sep =";"),
                  taxonomy_species = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep =";")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-taxonomy) %>%
    dplyr::relocate(sample, gene, OTU_ID,number_of_hits) %>%
    dplyr::relocate(sequence,.after = -1) %>%
    dplyr::rename(Sample_name = sample, Gene = gene, Count = number_of_hits, Sequence = sequence) %>%
    as.data.frame()
  
  singlem_otu_data.df
  
}

taxa_counts_and_abundances_from_singlem_otu_data <- function(singlem.df, taxonomy_map.df, sample_ids.v){
  # Combine the ASV counts back with the taxonomy data
  singlem_taxonomy_merged.df <- dplyr::left_join(singlem.df, taxonomy_map.df, by = "OTU_ID")
  
  # Run the taxa summarising function at each taxonomy level 
  tax_data.l <- generate_tax_level_data(singlem_taxonomy_merged.df, 
                                        tax_string_levels = c("OTU_ID", "taxonomy_species", "taxonomy_genus", "taxonomy_family", "taxonomy_class", "taxonomy_order", "taxonomy_phylum"),
                                        sample_ids = sample_ids.v,
                                        remove_zero_row_entries = T)
  tax_data.l
}

extract_counts_and_abundances_from_singlem_otu_data <- function(singlem_otu_data.df, gene = NULL){
  genes.v <- singlem_otu_data.df %>% dplyr::pull(Gene) %>% unique()
  
  if (is.null(gene)){
    error_message <- paste0("Nominate a marker gene to extract data for. Choose from:\n ", paste(genes.v, collapse = "\n "))
    stop(error_message)
  }
  if (!gene %in% genes.v){
    error_message <- paste0("Gene ", gene, " is not in the SingleM input data. Choose from:\n ", paste(genes.v, collapse = "\n "))
    stop(error_message)
  }
  
  gene_taxonomy_data.df <-
    singlem_otu_data.df %>%
    dplyr::filter(Gene == gene) %>%
    dplyr::select(-Sample_name) %>%
    unique() %>%
    as.data.frame()
  
  singlem_gene.df <- 
    singlem_otu_data.df %>%
    dplyr::filter(Gene == gene) %>% 
    tidyr::pivot_wider(id_cols = OTU_ID,names_from = Sample_name, values_from = Count, values_fill = 0) %>%
    as.data.frame()
  
  sample_ids_shotgun.v <- colnames(df2m(singlem_gene.df))
  singlem_tax_data.l <- taxa_counts_and_abundances_from_singlem_otu_data(singlem_gene.df, gene_taxonomy_data.df,sample_ids_shotgun.v)
  singlem_tax_data.l
}

#' Take a dataframe with the feature counts (first column ID and remaining columns the sample counts and taxonomy labels),
#' a list of the samples and the taxa labels you wish to sum over.
#' @return a list with the counts and relative abundances at each specified taxonomy level
generate_tax_level_data <- function(feature_count_taxonomy_data.df, sample_ids, tax_string_levels, remove_zero_row_entries = F,
                                    normalise_abundances = T){
  
  df2m <- function(mydataframe){
    mymatrix <- mydataframe
    rownames(mymatrix) <- mydataframe[,1]
    mymatrix[,1] <- NULL
    mymatrix <- as.matrix(mymatrix)
    mymatrix
  }
  
  m2df <- function(mymatrix, column_name = "Row_variable"){
    mydf <- as.data.frame(mymatrix)
    cur_names <- names(mydf)
    mydf[, column_name] <- rownames(mydf)
    rownames(mydf) <- NULL
    mydf <- mydf[,c(column_name,cur_names)]
    return(mydf)
  }
  output <- list()
  for (tax_string_level in tax_string_levels){
    counts.df <- reshape2::melt(feature_count_taxonomy_data.df[c(tax_string_level, sample_ids)], measure.vars = c(sample_ids))
    counts.df <- reshape2::dcast(counts.df, get(tax_string_level) ~ variable, sum)
    names(counts.df)[1] <- tax_string_level
    
    # Now create the relative abundance matrix at the current taxa level
    abundances.m <- counts.df
    
    rownames(abundances.m) <- abundances.m[[tax_string_level]]
    abundances.m[tax_string_level] <- NULL
    abundances.m <- as.matrix(abundances.m)
    if (normalise_abundances){
      abundances.m <- t(t(abundances.m[,,drop =F]) / colSums(abundances.m[,,drop =F]))* 100
    }
    abundances.m[is.nan(abundances.m)] <- 0
    
    counts.m <- df2m(counts.df)
    if (remove_zero_row_entries == T){
      abundances.m <- abundances.m[apply(abundances.m, 1, max) > 0,,drop =F]
      counts.m <- counts.m[apply(counts.m, 1, max) > 0,,drop =F]
    }
    output[[tax_string_level]] <- list("counts" = counts.m, "abundances" = abundances.m)
  }
  output
}

# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------

# Run the helper function "load_and_process_singleM_otu_data" to load and process the 'full' OTU output from SingleM.
# E.g.:
# gene	sample	sequence	num_hits	coverage	taxonomy
# S3.1.ribosomal_protein_L2_rplB	my_sample	GAAATGAAGCCGGGCAAGGGCGCGCAGCTCGCGCGAAGCGCCGGCGCCTCGGTGCAGATC	87	147.86	Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Rariloculales; f__Rariloculaceae
# S3.1.ribosomal_protein_L2_rplB	my_sample	GAAATGCGTCCTGGCAAGGGTGGTCAGCTGGCACGCTCTGCCGGTACCTCGGTTCAGGTG	1	1.64	Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Pseudomonadales; f__Alcanivoracaceae; g__Isoalcanivorax
# ...
# 
# This will produce a single table listing, per sample and marker gene, each OTU (ID + sequences), and its read counts, and corresponding lineage information
singlem_otu_data.df <- load_and_process_singleM_otu_data(singleM_output_file = "example_data/metagenome.otu_table.tsv")

# Run the helper function "extract_counts_and_abundances_from_singlem_otu_data" to process the table produced by "load_and_process_singleM_otu_data".
# This will extract the read counts and relative abundances (calculated from the coverage values) for a nominated marker gene, and return
# a list of counts and abundances at each taxonomy level (OTU_ID -> taxonomy_phylum) for each sample
# E.g.
# $taxonomy_phylum
#    -----$counts
#    -----$abundances
#                                   my_sample
# d__Archaea;p__Halobacteriota      0.04789272
# d__Archaea;p__Thermoplasmatota    0.14367816
# d__Bacteria;p__Acidobacteriota    0.19157088
singlem_counts_abundances_specific_marker.l <- extract_counts_and_abundances_from_singlem_otu_data(singlem_otu_data.df, gene = "S3.1.ribosomal_protein_L2_rplB")

# Run the helper function "load_and_process_singleM_condensed_data" to load and process the 'condensed' output from SingleM.
# This output describes the mean coverage across each lineage, calculated across all marker genes.
# E.g.:
# sample	coverage	taxonomy
# my_sample	8.21	Root; d__Archaea
# my_sample	35.44	Root; d__Bacteria
# my_sample	0.23	Root; d__Archaea; p__Halobacteriota
# ...
# 
# Will return a list of coverages and abundances at each taxonomy level 
# (OTU_ID -> taxonomy_phylum) for each sample
# $taxonomy_phylum
#    -----$coverage
#    -----$abundances
singlem_condensed_data.l <- load_and_process_singleM_condensed_data(singleM_output_file = "example_data/metagenome_condensed.tsv")

# Load the microbial fractions output file from SingleM
# E.g.:
# sample	bacterial_archaeal_bases	metagenome_size	read_fraction	average_bacterial_archaeal_genome_size	warning
# my_sample	12450640677	16039119487.0	71.99	3643684	
singlem_microbial_fractions.df <- 
  read.delim("example_data/microbial_fractions.tsv") %>%
  dplyr::rename("Sample_ID" = sample) %>%
  dplyr::select(-warning)

# Extract the coverage (cov) and abundance data from 'singlem_condensed_data.l'
species_cov_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_species$coverage, "taxonomy_species")
genus_cov_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_genus$coverage, "taxonomy_genus")
family_cov_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_family$coverage, "taxonomy_family")
order_cov_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_order$coverage, "taxonomy_order")
class_cov_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_class$coverage, "taxonomy_class")
phylum_cov_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_phylum$coverage, "taxonomy_phylum")

species_rel_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_species$abundances, "taxonomy_species")
genus_rel_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_genus$abundances, "taxonomy_genus")
family_rel_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_family$abundances, "taxonomy_family")
order_rel_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_order$abundances, "taxonomy_order")
class_rel_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_class$abundances, "taxonomy_class")
phylum_rel_singlem.df <- m2df(singlem_condensed_data.l$taxonomy_phylum$abundances, "taxonomy_phylum")

# Scale the relative abundances by the microbial fractions
# Requires the sample names to match across the various singleM results
scale_singlem_dataframe <- function(mydata.df){
  mydata.df %>%
    tidyr::pivot_longer(cols = -1,
                        names_to = "Sample_ID", values_to = "Relative_abundance") %>%
    dplyr::left_join(., singlem_microbial_fractions.df %>% dplyr::select(Sample_ID, read_fraction), by = "Sample_ID") %>%
    dplyr::mutate(Relative_abundance = Relative_abundance * read_fraction/100) %>%
    dplyr::select(-read_fraction) %>%
    tidyr::pivot_wider(names_from = Sample_ID, values_from = Relative_abundance, values_fill = 0) %>%
    dplyr::arrange(1) %>%
    as.data.frame()
}

species_rel_singlem_scaled.df <- scale_singlem_dataframe(species_rel_singlem.df)
genus_rel_singlem_scaled.df <- scale_singlem_dataframe(genus_rel_singlem.df)
family_rel_singlem_scaled.df <- scale_singlem_dataframe(family_rel_singlem.df)
order_rel_singlem_scaled.df <- scale_singlem_dataframe(order_rel_singlem.df)
class_rel_singlem_scaled.df <- scale_singlem_dataframe(class_rel_singlem.df)
phylum_rel_singlem_scaled.df <- scale_singlem_dataframe(phylum_rel_singlem.df)


# Write the coverages and abundances to file
tax_levels.v <- c("taxonomy_species", "taxonomy_genus", "taxonomy_family", "taxonomy_class", "taxonomy_order", "taxonomy_phylum")

coverages.l <- list()
abundances.l <- list()
abundances_scaled.l <- list()
for (tax_string_level in tax_levels.v){
  coverages.l[[tax_string_level]]  <-  m2df(singlem_condensed_data.l[[tax_string_level]][["coverage"]], tax_string_level)
  abundances.l[[tax_string_level]]  <-  m2df(singlem_condensed_data.l[[tax_string_level]][["abundances"]], tax_string_level)
  abundances_scaled.l[[tax_string_level]] <- get(paste0(gsub("taxonomy_", "", tax_string_level),"_rel_singlem_scaled.df"))
  
}
# openxlsx::write.xlsx(c(coverages.l,list("microbial_fractions" = singlem_microbial_fractions.df)), file = "SingleM_coverages.xlsx" ,overwrite = T)
# openxlsx::write.xlsx(c(abundances.l, list("microbial_fractions" = singlem_microbial_fractions.df)), file = "SingleM_abundances.xlsx" ,overwrite = T)
# openxlsx::write.xlsx(c(abundances_scaled.l, list("microbial_fractions" = singlem_microbial_fractions.df)), file = "SingleM_scaled_abundances.xlsx" ,overwrite = T)
write_neat_xlsx(c(coverages.l,list("microbial_fractions" = singlem_microbial_fractions.df)), output_filename = "SingleM_coverages.xlsx")
write_neat_xlsx(c(abundances.l,list("microbial_fractions" = singlem_microbial_fractions.df)), output_filename = "SingleM_abundances.xlsx")
write_neat_xlsx(c(abundances_scaled.l,list("microbial_fractions" = singlem_microbial_fractions.df)), output_filename = "SingleM_scaled_abundances.xlsx")
