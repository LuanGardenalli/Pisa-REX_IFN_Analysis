library(ComplexHeatmap)
library(circlize)

P <- 0.05
facet <- "REX"

#Extract all accessions of proteins reaching significance in at least one group:
extract.sig.proteins <- function(data.peptides,
                                 data.proteins,
                                 cell.line,
                                 type,
                                 facet){
  rowID  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                   cell.line,
                                   type,
                                   facet) %>% 
    filter(adj.P.Val < P) %>%
    select(ID) %>% pull()
}
{
  TA <- extract.sig.proteins(data.peptides,
                             data.proteins,
                             cell.line="THP1",
                             type="Alpha",
                             facet)
  
  TB <- extract.sig.proteins(data.peptides,
                             data.proteins,
                             cell.line="THP1",
                             type="Beta",
                             facet)
  
  TG <- extract.sig.proteins(data.peptides,
                             data.proteins,
                             cell.line="THP1",
                             type="Gamma",
                             facet)
  
  HA <- extract.sig.proteins(data.peptides,
                             data.proteins,
                             cell.line="HL60",
                             type="Alpha",
                             facet)
  
  HB <- extract.sig.proteins(data.peptides,
                             data.proteins,
                             cell.line="HL60",
                             type="Beta",
                             facet)
  
  HG <- extract.sig.proteins(data.peptides,
                             data.proteins,
                             cell.line="HL60",
                             type="Gamma",
                             facet)
  accessions <- unique(c(TA, TB, TG, HA, HB, HG))
}

#Colapse peptides into proteins, calculate Log2 Ratios for extracted accessions:
Log2Ratios <- function(data.peptides,
                       data.proteins,
                       cell.line,
                       type,
                       facet,
                       accessions){
  rowID  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                   cell.line,
                                   type,
                                   facet) %>% 
    arrange(Rank) %>%  
    distinct(ID, .keep_all = T) %>% 
    filter(ID %in% accessions) %>% 
    select(peptideRowID) %>% pull()
  replicates <- pisaRex.preprocess(data.peptides,
                                   data.proteins,
                                   cell.line,
                                   type,
                                   facet) %>% 
    filter(peptideRowID %in% rowID)
  log2ratios <- replicates %>% 
    select(`Master Protein Accessions`, contains("Abundance")) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    select(`Master Protein Accessions`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~ control_mean - .x)) %>% 
    select(!control_mean)
}
{
  TA <- Log2Ratios(data.peptides,
                   data.proteins,
                   cell.line="THP1",
                   type="Alpha",
                   facet,
                   accessions)
  
  TB <- Log2Ratios(data.peptides,
                   data.proteins,
                   cell.line="THP1",
                   type="Beta",
                   facet,
                   accessions)
  
  TG <- Log2Ratios(data.peptides,
                   data.proteins,
                   cell.line="THP1",
                   type="Gamma",
                   facet,
                   accessions)
  
  HA <- Log2Ratios(data.peptides,
                   data.proteins,
                   cell.line="HL60",
                   type="Alpha",
                   facet,
                   accessions)
  
  HB <- Log2Ratios(data.peptides,
                   data.proteins,
                   cell.line="HL60",
                   type="Beta",
                   facet,
                   accessions)
  
  HG <- Log2Ratios(data.peptides,
                   data.proteins,
                   cell.line="HL60",
                   type="Gamma",
                   facet,
                   accessions)
  
  #Combine data into one matrix:
  heatmap.df.T <- merge(TA, 
                        TB, 
                        by= "Master Protein Accessions",
                        all.x=T, all.y=T)
  heatmap.df.T <- merge(heatmap.df.T, 
                        TG, 
                        by= "Master Protein Accessions", 
                        all.x=T, all.y=T)
  
  heatmap.df.H <- merge(HA, 
                        HB, 
                        by= "Master Protein Accessions", 
                        all.x=T, all.y=T)
  heatmap.df.H <- merge(heatmap.df.H, 
                        HG, 
                        by= "Master Protein Accessions", 
                        all.x=T, all.y=T)
  
  heatmap.df <- merge(heatmap.df.T,
                      heatmap.df.H, 
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
