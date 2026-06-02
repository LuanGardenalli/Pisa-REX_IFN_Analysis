library(tidyverse)
library(readxl)
data.proteins <- read_excel("/path/data")
data.peptides <- read_csv("/path/data")

#Pre processing data:
pisaRex.preprocess <- function(data.peptides,
                               data.proteins,
                               cell.line,
                               type,
                               facet){
  library(tidyverse)
  
  #Index peptides, filter out contaminants, select columns of interest, and apply normalisations for the specified facet:
  if (facet == "REX") {
    data <- data.peptides
    data <- data %>% 
      mutate(peptideRowID = row_number()) %>% 
      filter(Contaminant=="FALSE") %>% 
      dplyr::select(peptideRowID,`Annotated Sequence` : `Modifications in Master Proteins`, 
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains(facet)) %>% 
      drop_na(starts_with("Abundance:")) %>% 
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x)))
    
    #Non-Cys peptide normalization:
    cys.peptides <- data %>% 
      filter(grepl('C', `Annotated Sequence`)) %>% 
      filter(grepl("Carbamidomethyl", Modifications))
    non.cys.peptides <- data %>% 
      filter(!grepl('C', `Annotated Sequence`))
    non.cys.accessions <- non.cys.peptides %>%
      dplyr::select(`Master Protein Accessions`) %>% 
      pull()
    non.cys.id <- non.cys.peptides %>% 
      mutate(non.cys.id = row_number())
    cys.id <- cys.peptides %>% 
      mutate(cys.id = row_number()) %>% 
      inner_join(non.cys.id %>% dplyr::select(non.cys.id, `Master Protein Accessions`), relationship = "many-to-many")
    denominators <- cys.id %>%
      dplyr::select(cys.id, non.cys.id) %>% 
      left_join(non.cys.id, by = "non.cys.id") %>%
      dplyr::select(cys.id , starts_with("Abundance:")) %>% 
      group_by(cys.id) %>% 
      summarise_all(~ sum(.))
    data <- cys.id %>% 
      dplyr::select(!c(non.cys.id, cys.id)) %>% 
      distinct() %>% 
      mutate(across(starts_with("Abundance:"), ~ .x / denominators$.x))
    data$`Master Protein Accessions` <- sub("-.*", "", data$`Master Protein Accessions`)
    data <- data %>% mutate(`Master Protein Accessions` = word(`Master Protein Accessions`, 1, sep = ";"))
    
    #Normalisations for PISA/EXP:
  } else if (facet == "PISA") {
    denominators <- data.proteins %>% 
      filter(`Protein FDR Confidence: Combined`=="High" & Contaminant=="FALSE") %>% 
      dplyr::select(`Accession`, `Gene Symbol`,
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains("EXP")) %>% 
      drop_na() %>% 
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x))) %>% 
      mutate(mean.control = rowMeans(pick(3:5), na.rm = TRUE),
             mean.sample = rowMeans(pick(6:8), na.rm = TRUE)) %>% 
      dplyr::select(Accession, `Gene Symbol`, mean.control, mean.sample)
    
    data <- data.proteins %>% 
      filter(`Protein FDR Confidence: Combined`=="High" & Contaminant=="FALSE") %>% 
      dplyr::select(`Accession`, `Gene Symbol`,
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains(facet)) %>% 
      drop_na() %>%
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x)))
    
    data <- merge(data, denominators, .by = c("Accession", "Gene Symbol")) %>% 
      mutate(across(contains("Control") & contains("Abundance:"), ~ .x / mean.control),
             across(contains("Sample") & contains("Abundance:"), ~ .x / mean.sample)) %>% 
      dplyr::select(!mean.control & !mean.sample)
    
  } else if ( facet == "EXP") {
    data <- data.proteins
    data <- data %>% 
      filter(`Protein FDR Confidence: Combined`=="High" & Contaminant=="FALSE") %>% 
      dplyr::select(`Accession`, `Gene Symbol`,
                    contains("Abundance:") & 
                      contains(cell.line) & 
                      contains(type) & 
                      contains(facet)) %>% 
      drop_na() %>% 
      mutate(across(starts_with("Abundance:"), ~.x / sum(.x)))
  }
  
  #CV filter and log2 transformation:
  Control <- names(data %>% dplyr::select(contains("Control")))
  Sample <- names(data %>% dplyr::select(contains("Sample")))
  data <- data %>% 
    mutate(`Control CV (%)` = (apply(.[Control], 1, sd, na.rm=TRUE) / rowMeans(.[Control], na.rm=TRUE)) * 100, 
           `Sample CV (%)` = (apply(.[Sample], 1, sd, na.rm=TRUE) / rowMeans(.[Sample], na.rm=TRUE)) * 100) %>% 
    filter(`Control CV (%)` < 30 & `Sample CV (%)` < 30)
  data[c(Control, Sample)][data[c(Control, Sample)] == 0] <- NA
  data[c(Control, Sample)] <- log2(data[c(Control, Sample)])
  return(data)
}
