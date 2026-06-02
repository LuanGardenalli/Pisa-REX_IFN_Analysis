library(ComplexHeatmap)
library(circlize)

P <- 0.05
#Extract all accessions of proteins reaching significance in at least one group:
facet <- "REX"
{
  accessions.TRA  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                            cell.line="THP1",
                                            type="Alpha",
                                            facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(ID) %>% pull()
  
  accessions.TRB  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                            cell.line="THP1",
                                            type="Beta",
                                            facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(ID) %>% pull()
  
  accessions.TRG  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                            cell.line="THP1",
                                            type="Gamma",
                                            facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(ID) %>% pull()
  
  accessions.HRA  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                            cell.line="HL60",
                                            type="Alpha",
                                            facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(ID) %>% pull()
  
  accessions.HRB  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                            cell.line="HL60",
                                            type="Beta",
                                            facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(ID) %>% pull()
  
  accessions.HRG  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                            cell.line="HL60",
                                            type="Gamma",
                                            facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(ID) %>% pull()
  
  accessions <- unique(c(accessions.TRA, accessions.TRB, accessions.TRG, 
                         accessions.HRA, accessions.HRB, accessions.HRG))
}

#Colapse peptides into proteins, calculate Log2 Ratios for extracted accessions:
{
  rowID.TRA  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="THP1",
                                       type="Alpha",
                                       facet=facet) %>% 
    arrange(Rank) %>%  
    distinct(ID, .keep_all = T) %>% 
    filter(ID %in% accessions) %>% 
    select(peptideRowID) %>% pull()
  TRA.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="THP1",
                                       type="Alpha",
                                       facet=facet) %>% 
    filter(peptideRowID %in% rowID.TRA)
  TRA.log2ratios <- TRA.replicates %>% 
    select(`Master Protein Accessions`, contains("Abundance")) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    select(`Master Protein Accessions`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~control_mean - .x)) %>% 
    select(!control_mean)
  
  rowID.TRB  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="THP1",
                                       type="Beta",
                                       facet=facet) %>% 
    arrange(Rank) %>%  
    distinct(ID, .keep_all = T) %>% 
    filter(ID %in% accessions) %>% 
    select(peptideRowID) %>% pull()
  TRB.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="THP1",
                                       type="Beta",
                                       facet=facet) %>% 
    filter(peptideRowID %in% rowID.TRB)
  TRB.log2ratios <- TRB.replicates %>% 
    select(`Master Protein Accessions`, contains("Abundance")) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    select(`Master Protein Accessions`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~control_mean - .x)) %>% 
    select(!control_mean)
  
  rowID.TRG  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="THP1",
                                       type="Gamma",
                                       facet=facet) %>% 
    arrange(Rank) %>%  
    distinct(ID, .keep_all = T) %>% 
    filter(ID %in% accessions) %>% 
    select(peptideRowID) %>% pull()
  TRG.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="THP1",
                                       type="Gamma",
                                       facet=facet) %>% 
    filter(peptideRowID %in% rowID.TRG)
  TRG.log2ratios <- TRG.replicates %>% 
    select(`Master Protein Accessions`, contains("Abundance")) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    select(`Master Protein Accessions`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~control_mean - .x)) %>% 
    select(!control_mean)
  
  rowID.HRA  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="HL60",
                                       type="Alpha",
                                       facet=facet) %>% 
    arrange(Rank) %>%  
    distinct(ID, .keep_all = T) %>% 
    filter(ID %in% accessions) %>% 
    select(peptideRowID) %>% pull()
  HRA.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="HL60",
                                       type="Alpha",
                                       facet=facet) %>% 
    filter(peptideRowID %in% rowID.HRA)
  HRA.log2ratios <- HRA.replicates %>% 
    select(`Master Protein Accessions`, contains("Abundance")) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    select(`Master Protein Accessions`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~control_mean - .x)) %>% 
    select(!control_mean)
  
  rowID.HRB  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="HL60",
                                       type="Beta",
                                       facet=facet) %>% 
    arrange(Rank) %>%  
    distinct(ID, .keep_all = T) %>% 
    filter(ID %in% accessions) %>% 
    select(peptideRowID) %>% pull()
  HRB.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="HL60",
                                       type="Beta",
                                       facet=facet) %>% 
    filter(peptideRowID %in% rowID.HRB)
  HRB.log2ratios <- HRB.replicates %>% 
    select(`Master Protein Accessions`, contains("Abundance")) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    select(`Master Protein Accessions`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~control_mean - .x)) %>% 
    select(!control_mean)
  
  rowID.HRG  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="HL60",
                                       type="Gamma",
                                       facet=facet) %>% 
    arrange(Rank) %>%  
    distinct(ID, .keep_all = T) %>% 
    filter(ID %in% accessions) %>% 
    select(peptideRowID) %>% pull()
  HRG.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="HL60",
                                       type="Gamma",
                                       facet=facet) %>% 
    filter(peptideRowID %in% rowID.HRG)
  HRG.log2ratios <- HRG.replicates %>% 
    select(`Master Protein Accessions`, contains("Abundance")) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    select(`Master Protein Accessions`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~control_mean - .x)) %>% 
    select(!control_mean)

  
  heatmap.df.TR <- merge(TRA.log2ratios, 
                         TRB.log2ratios , 
                         by= "Master Protein Accessions", 
                         all.x=T, all.y=T)
  heatmap.df.TR <- merge(heatmap.df.TR, 
                         TRG.log2ratios , 
                         by= "Master Protein Accessions", 
                         all.x=T, all.y=T)
  
  heatmap.df.HR <- merge(HRA.log2ratios , 
                         HRB.log2ratios , 
                         by= "Master Protein Accessions", 
                         all.x=T, all.y=T)
  heatmap.df.HR <- merge(heatmap.df.HR, 
                         HRG.log2ratios , 
                         by= "Master Protein Accessions", 
                         all.x=T, all.y=T)
  
  heatmap.df <- merge(heatmap.df.TR,
                      heatmap.df.HR, 
                      by= "Master Protein Accessions", 
                      all.x=T, all.y=T)
  heatmap.df[is.na(heatmap.df)] <- 0
  
  heatmap.matrix <- heatmap.df %>%
    column_to_rownames("Master Protein Accessions") %>%
    as.matrix()
}

#Assign clusters based on the Leiden clustering file:
heatmap.df$Cluster <- read_csv("Leiden_Clusters_REX.txt", col_names = F) %>% pull()

#Draw Heatmap:
{
  #Terms annotation:
  enrichment_terms <- c(
    "1:\nIFN signalling.\nImmune response.",
    
    "2:\nMitochondrial ribosome.\n",
    
    "3:\nMethylation.",
    
    "4:\nProteasome complex.\nRNA processing.\nAminoacyl-tRNA synthetase.",
    
    "5:\nER-Golgi transport.",
    
    "6:\nIFN-G response.\nGTP binding.",
    
    "7:\nIFN signalling.\nImmune response.",
    
    "8:\nAntiviral defense.")
  
  
  names <- c("", "THP1, Alpha", "",
             "", "THP1, Beta", "",
             "", "THP1, Gamma", "",
             
             "", "HL60, Alpha", "",
             "", "HL60, Beta", "",
             "", "HL60, Gamma", "")
  
  my_colors <- colorRamp2(c(-1, 0, 1), 
                          c("#2166AC", "#F7F7F7","#B2182B"))
  
  ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                      width = max_text_width(enrichment_terms)))
  
  HM.REX <- Heatmap(heatmap.matrix,
                    name = "Log2\nRatio",
                    col = my_colors,
                    cluster_rows = T,
                    cluster_columns = F,
                    
                    row_split = heatmap.df$Cluster,       
                    row_title = NULL,              
                    show_row_names = FALSE,        
                    right_annotation = ha,
                    
                    show_column_names = TRUE,
                    column_labels = names,
                    column_names_rot = 90,
                    column_names_gp = gpar(fontsize = 10),
                    show_heatmap_legend = FALSE)
}
draw(HM.REX)
for(i in 1:8) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(enrichment_terms[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left",
              gp = gpar(
                fontsize = 9, 
                fontfamily = "sans", 
                fontface = "bold"))
  })
}
