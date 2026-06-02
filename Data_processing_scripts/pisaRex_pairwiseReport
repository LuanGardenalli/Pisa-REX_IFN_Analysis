#Statistical analysis with limma:
pisaRex.pairwiseReport <- function(data.peptides,
                                   data.proteins,
                                   cell.line,
                                   type,
                                   facet){
  library(limma)
  
  data <- pisaRex.preprocess(data.peptides, data.proteins, cell.line, type, facet)
  
  Control <- names(data %>% 
                     dplyr::select(contains("Control")) %>% 
                     dplyr::select(!`Control CV (%)`))
  
  Sample <- names(data %>% 
                    dplyr::select(contains("Sample")) %>% 
                    dplyr::select(!`Sample CV (%)`)) 
  
  #Differential Expression with Limma
  samples = grep('Abundance:', colnames(data), value=T)
  treatment <- c(rep("Control", length(Control)), rep("Sample", length(Sample)))
  replicates <- c(rep(c(1:length(Sample)), 2))
  sample.table <- tibble("Sample" = samples,
                         "Treatment" = treatment,
                         "Replicates" = replicates)
  sample.table <- column_to_rownames(sample.table, var = "Sample")
  design <- sapply(unique(sample.table$Treatment), function(x) sample.table$Treatment %in% x)
  design = cbind(design)+0
  
  #This accounts for the "inverted" experiment design of the REX facet:
  if (facet == "REX") {
    contrast = makeContrasts(
      Control - Sample,
      levels = design)
  } else {
    contrast = makeContrasts(
      Sample - Control,
      levels = design)
  }
  fit <- lmFit(data[c(Control, Sample)], design)
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)
  if (facet == "REX") {
    DEGs <- topTable(fit, 
                     number = Inf, 
                     genelist = data$`Master Protein Accessions`, 
                     adjust.method = "BH", 
                     sort.by = "none", 
                     p.value = 1, 
                     lfc = 0, 
                     confint = F) %>% 
      mutate(`ID` = word(ID, 1, sep = ";"))
  } else {
    DEGs <- topTable(fit, 
                     number = Inf, 
                     genelist = data$`Gene Symbol`, 
                     adjust.method = "BH", 
                     sort.by = "none", 
                     p.value = 1, 
                     lfc = 0, 
                     confint = F) %>% 
      mutate(`ID` = word(ID, 1, sep = ";"))
  }
  
  #Calculate protein rankings:
  DEGs <- DEGs %>% 
    mutate(FC.rank = rank(-abs(logFC))) %>% 
    mutate(P.value.rank = rank(adj.P.Val)) %>% 
    mutate(Rank_sum = P.value.rank + FC.rank) %>% 
    mutate(Rank = rank(Rank_sum)) %>% 
    dplyr::select(!FC.rank & !P.value.rank & !Rank_sum)
  
  #Translate accessions to gene symbols for REX facet: 
  if (facet == "REX") {
    DEGs <- DEGs %>% mutate(peptideRowID = data$peptideRowID)
    DEGs$ID <- sub("-.*", "", DEGs$ID)
    gene.map <- as_tibble(clusterProfiler::bitr(DEGs$ID, 
                               fromType = "UNIPROT",  
                               toType = c("SYMBOL"),    
                               OrgDb = org.Hs.eg.db::org.Hs.eg.db))
    colnames(gene.map) <- c("ID", "Symbol")
    DEGs <- merge(DEGs, gene.map, by.x = "ID", by.y = "ID", all.x = TRUE)
    
    rowID <- data %>% 
      dplyr::select(peptideRowID, `Positions in Master Proteins`, `Annotated Sequence`)
    
    DEGs <- merge(DEGs, rowID, by.x = "peptideRowID", by.y = "peptideRowID", all.x = F)
  } else {
    DEGs <- DEGs %>% 
      rename(Symbol = ID) %>% 
      mutate(proteinRowID = row_number())
  }
  return(DEGs %>% 
           drop_na())
}

#Uses protein accession as ID instead of gene symbol:
pisaRex.pairwiseReport.accession <- function(data.peptides,
                                   data.proteins,
                                   cell.line,
                                   type,
                                   facet){
  library(limma)
  
  data <- pisaRex.preprocess(data.peptides, data.proteins, cell.line, type, facet)
  
  Control <- names(data %>% 
                     dplyr::select(contains("Control")) %>% 
                     dplyr::select(!`Control CV (%)`))
  
  Sample <- names(data %>% 
                    dplyr::select(contains("Sample")) %>% 
                    dplyr::select(!`Sample CV (%)`)) 
  
  #Differential Expression with Limma
  samples = grep('Abundance:', colnames(data), value=T)
  treatment <- c(rep("Control", length(Control)), rep("Sample", length(Sample)))
  replicates <- c(rep(c(1:length(Sample)), 2))
  sample.table <- tibble("Sample" = samples,
                         "Treatment" = treatment,
                         "Replicates" = replicates)
  sample.table <- column_to_rownames(sample.table, var = "Sample")
  design <- sapply(unique(sample.table$Treatment), function(x) sample.table$Treatment %in% x)
  design = cbind(design)+0
  
  #This accounts for the "inverted" experiment design of the REX facet:
  if (facet == "REX") {
    contrast = makeContrasts(
      Control - Sample,
      levels = design)
  } else {
    contrast = makeContrasts(
      Sample - Control,
      levels = design)
  }
  fit <- lmFit(data[c(Control, Sample)], design)
  fit = contrasts.fit(fit, contrast)
  fit = eBayes(fit)
  if (facet == "REX") {
    DEGs <- topTable(fit, 
                     number = Inf, 
                     genelist = data$`Master Protein Accessions`, 
                     adjust.method = "BH", 
                     sort.by = "none", 
                     p.value = 1, 
                     lfc = 0, 
                     confint = F)
  } else {
    DEGs <- topTable(fit, 
                     number = Inf, 
                     genelist = data$Accession, 
                     adjust.method = "BH", 
                     sort.by = "none", 
                     p.value = 1, 
                     lfc = 0, 
                     confint = F)
  }
  
  #Calculate protein rankings:
  DEGs <- DEGs %>% 
    mutate(FC.rank = rank(-abs(logFC))) %>% 
    mutate(P.value.rank = rank(adj.P.Val)) %>% 
    mutate(Rank_sum = P.value.rank + FC.rank) %>% 
    mutate(Rank = rank(Rank_sum)) %>% 
    dplyr::select(!FC.rank & !P.value.rank & !Rank_sum)
  
  if (facet == "REX") {
    DEGs <- DEGs %>% 
      mutate(peptideRowID = data$peptideRowID) %>% 
      rename(Accession = ID)
  } else {
    DEGs <- DEGs %>% 
      rename(Accession = ID) %>% 
      mutate(proteinRowID = row_number())
  }
    return(DEGs %>% 
             drop_na())
}
