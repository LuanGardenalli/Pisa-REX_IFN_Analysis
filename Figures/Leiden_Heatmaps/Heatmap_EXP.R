library(ComplexHeatmap)
library(circlize)

P <- 0.05
facet <- "EXP"

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
    select(proteinRowID) %>% pull()
  accession <- pisaRex.preprocess(data.peptides, data.proteins, 
                                  cell.line,
                                  type,
                                  facet) %>% 
    filter(row_number() %in% rowID) %>% 
    select(Accession) %>% pull()
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

#Calculate Log2 Ratios for extracted accessions:
Log2Ratios <- function(data.peptides,
                       data.proteins,
                       cell.line,
                       type,
                       facet,
                       accessions){
  replicates <- pisaRex.preprocess(data.peptides,
                                   data.proteins,
                                   cell.line,
                                   type,
                                   facet) %>% 
    filter(Accession %in% accessions)
  log2ratios <- replicates %>% 
    select(Accession, 3:8) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(3:4), na.rm = TRUE)) %>% 
    select(Accession, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
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
                        by= "Accession",
                        all.x=T, all.y=T)
  heatmap.df.T <- merge(heatmap.df.T, 
                        TG, 
                        by= "Accession", 
                        all.x=T, all.y=T)
  
  heatmap.df.H <- merge(HA, 
                        HB, 
                        by= "Accession", 
                        all.x=T, all.y=T)
  heatmap.df.H <- merge(heatmap.df.H, 
                        HG, 
                        by= "Accession", 
                        all.x=T, all.y=T)
  
  heatmap.df <- merge(heatmap.df.T,
                      heatmap.df.H, 
                      by= "Accession", 
                      all.x=T, all.y=T)
  heatmap.df[is.na(heatmap.df)] <- 0
  
  heatmap.matrix <- heatmap.df %>%
    column_to_rownames("Accession") %>%
    as.matrix()
}

#Assign clusters based on the Leiden clustering file:
heatmap.df$Cluster <- read_csv("Leiden_Clusters_EXP.txt", col_names = F) %>% pull()

#Draw Heatmap:
{
  
  #Terms annotation:
  enrichment_terms <- c(
    "1:\nIFN signalling.\nImmune response.",
    
    "2:\nIFN-G response.\nGTP binding.",
    
    "3:\nIFN-G response.\nAntigen processing.\nT-cell regulation.",
    
    "4:\nProteasome activation.",
    
    "5:\nImmune response.\nAnti-viral response.",
    
    "6:\n---",
    
    "7:\nMetallothionein activity.\nResponse to metal ions.",
    
    "8:\nIFN alpha/beta signaling.\nProtein modification.")
  
  names <- c("","THP1, Alpha", "",
             "","THP1, Beta", "",
             "","THP1, Gamma", "",
             
             "","HL60, Alpha", "",
             "","HL60, Beta", "",
             "","HL60, Gamma", "")
  
  my_colors <- colorRamp2(c(-1, 0, 1), 
                          c("#2166AC", "#F7F7F7", "#B2182B"))
  
  ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                      width = max_text_width(enrichment_terms) + unit(4, "mm")))
  
  HM.EXP <- Heatmap(heatmap.matrix,
                    name = "Log2 Ratio:\nEXP",
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
draw(HM.EXP)
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
