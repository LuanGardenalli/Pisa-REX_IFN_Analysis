library(ComplexHeatmap)
library(circlize)

P <- 0.05
#Extract all accessions of proteins reaching significance in at least one group:
facet <- "PISA"
{
  rowID.TPA  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="THP1",
                                       type="Alpha",
                                       facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(proteinRowID) %>% pull()
  accession.TPA <- pisaRex.preprocess(data.peptides, data.proteins, 
                                      cell.line="THP1",
                                      type="Alpha",
                                      facet=facet) %>% 
    filter(row_number() %in% rowID.TPA) %>% 
    select(Accession) %>% pull()
  
  rowID.TPB  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="THP1",
                                       type="Beta",
                                       facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(proteinRowID) %>% pull()
  accession.TPB <- pisaRex.preprocess(data.peptides, data.proteins, 
                                      cell.line="THP1",
                                      type="Beta",
                                      facet=facet) %>% 
    filter(row_number() %in% rowID.TPB) %>% 
    select(Accession) %>% pull()
  
  rowID.TPG  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="THP1",
                                       type="Gamma",
                                       facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(proteinRowID) %>% pull()
  accession.TPG <- pisaRex.preprocess(data.peptides, data.proteins, 
                                      cell.line="THP1",
                                      type="Gamma",
                                      facet=facet) %>% 
    filter(row_number() %in% rowID.TPG) %>% 
    select(Accession) %>% pull()
  
  rowID.HPA  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="HL60",
                                       type="Alpha",
                                       facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(proteinRowID) %>% pull()
  accession.HPA <- pisaRex.preprocess(data.peptides, data.proteins, 
                                      cell.line="HL60",
                                      type="Alpha",
                                      facet=facet) %>% 
    filter(row_number() %in% rowID.HPA) %>% 
    select(Accession) %>% pull()
  
  rowID.HPB  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="HL60",
                                       type="Beta",
                                       facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(proteinRowID) %>% pull()
  accession.HPB <- pisaRex.preprocess(data.peptides, data.proteins, 
                                      cell.line="HL60",
                                      type="Beta",
                                      facet=facet) %>% 
    filter(row_number() %in% rowID.HPB) %>% 
    select(Accession) %>% pull()
  
  rowID.HPG  <- pisaRex.pairwiseReport(data.peptides, data.proteins, 
                                       cell.line="HL60",
                                       type="Gamma",
                                       facet=facet) %>% 
    filter(adj.P.Val < P) %>%
    select(proteinRowID) %>% pull()
  accession.HPG <- pisaRex.preprocess(data.peptides, data.proteins, 
                                      cell.line="HL60",
                                      type="Gamma",
                                      facet=facet) %>% 
    filter(row_number() %in% rowID.HPG) %>% 
    select(Accession) %>% pull()
  
  accessions <- unique(c(accession.TPA, accession.TPB, accession.TPG, accession.HPA, accession.HPB, accession.HPG))
}

#Calculate Log2 Ratios for extracted accessions:
{
  TPA.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="THP1",
                                       type="Alpha",
                                       facet=facet) %>% 
    filter(Accession %in% accessions)
  TPA.log2ratios <- TPA.replicates %>% 
    select(1:6) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(3:4), na.rm = TRUE)) %>% 
    select(Accession, `Gene Symbol`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
    select(!control_mean)
  
  TPB.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="THP1",
                                       type="Beta",
                                       facet=facet) %>% 
    filter(Accession %in% accessions)
  TPB.log2ratios <- TPB.replicates %>% 
    select(1:6) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(3:4), na.rm = TRUE)) %>% 
    select(Accession, `Gene Symbol`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
    select(!control_mean)
  
  TPG.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="THP1",
                                       type="Gamma",
                                       facet=facet) %>% 
    filter(Accession %in% accessions)
  TPG.log2ratios <- TPG.replicates %>% 
    select(1:6) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(3:4), na.rm = TRUE)) %>% 
    select(Accession, `Gene Symbol`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
    select(!control_mean)
  
  HPA.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="HL60",
                                       type="Alpha",
                                       facet=facet) %>% 
    filter(Accession %in% accessions)
  HPA.log2ratios <- HPA.replicates %>% 
    select(1:6) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(3:4), na.rm = TRUE)) %>% 
    select(Accession, `Gene Symbol`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
    select(!control_mean)
  
  HPB.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="HL60",
                                       type="Beta",
                                       facet=facet) %>% 
    filter(Accession %in% accessions)
  HPB.log2ratios <- HPB.replicates %>% 
    select(1:6) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(3:4), na.rm = TRUE)) %>% 
    select(Accession, `Gene Symbol`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
    select(!control_mean)
  
  HPG.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line="HL60",
                                       type="Gamma",
                                       facet=facet) %>% 
    filter(Accession %in% accessions)
  HPG.log2ratios <- HPG.replicates %>% 
    select(1:6) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(3:4), na.rm = TRUE)) %>% 
    select(Accession, `Gene Symbol`, contains("Sample"), control_mean) %>% 
    mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
    select(!control_mean)


  heatmap.df.TP <- merge(TPA.log2ratios %>% select(!`Gene Symbol`), 
                         TPB.log2ratios %>% select(!`Gene Symbol`), 
                         by= "Accession", 
                         all.x=T, all.y=T)
  heatmap.df.TP <- merge(heatmap.df.TP, 
                         TPG.log2ratios %>% select(!`Gene Symbol`), 
                         by= "Accession", 
                         all.x=T, all.y=T)
  
  heatmap.df.HP <- merge(HPA.log2ratios %>% select(!`Gene Symbol`), 
                         HPB.log2ratios %>% select(!`Gene Symbol`), 
                         by= "Accession", 
                         all.x=T, all.y=T)
  heatmap.df.HP <- merge(heatmap.df.HP, 
                         HPG.log2ratios %>% select(!`Gene Symbol`), 
                         by= "Accession", 
                         all.x=T, all.y=T)
  
  heatmap.df <- merge(heatmap.df.TP,
                      heatmap.df.HP, 
                      by= "Accession", 
                      all.x=T, all.y=T)
  heatmap.df[is.na(heatmap.df)] <- 0 
}

#Assign clusters based on the Leiden clustering file:
heatmap.df$Cluster <- read_csv("Leiden_Clusters_PISA.txt", col_names = F) %>% pull()

#Remove clusters with less than 5 proteins:
heatmap.df %>% group_by(Cluster) %>% 
  summarise(n = n())

heatmap.df <- heatmap.df %>% 
  filter(!Cluster == 15)

heatmap.matrix <- heatmap.df %>%
  select(!Cluster) %>% 
  column_to_rownames("Accession") %>%
  as.matrix()

#Draw Heatmap:
{
  #Terms annotation:
  enrichment_terms <- c(
    "1:\nLSU-rRNA.\nsnoRNP complex.\nRibosome biogenesis.",
    
    "2: IFN signalling.\n    Immune response.",
    
    "3: Intracellular organelle.\n    Cell cycle.",
    
    "4: Acetylation.",
    
    "5:\nVesicle tethering complex.\nMitochondrial inner membrane.",
    
    "6:\nGlycosaminoglycan degradation.\nMembrane-bounded organelle.",
    
    "7:\nCollagen.",
    
    "8:\nCornified envelope.\nFocal adhesion.",
    
    "9:\nIsocitrate dehydrogenase activity.\nCajal bodies.",
    
    "10:\nOxidative phosphorylation.",
    
    "11:\nCytoplasmic ribosomes.",
    
    "12:\nPeroxisomal matrix.",
    
    "13:\nMitochondrion organization.",
    
    "14:\nficolin-1-rich granule.")
  
  names <- c("THP1, Alpha", "",
             "THP1, Beta", "",
             "THP1, Gamma", "",
             
             "HL60, Alpha", "",
             "HL60, Beta", "",
             "HL60, Gamma", "")
  
  my_colors <- colorRamp2(c(-1, 0, 1), 
                          c("#2166AC", "#F7F7F7", "#B2182B"))
  
  ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                      width = max_text_width(enrichment_terms) + unit(4, "mm")))
  
  HM.PISA <- Heatmap(heatmap.matrix,
                     name = "Log2 Ratio:\nPISA",
                     col = my_colors,
                     cluster_rows = T,
                     cluster_columns = F,
                     
                     row_split = heatmap.df$Cluster,       
                     row_title = NULL,              
                     show_row_names = FALSE,        
                     right_annotation = ha,
                     
                     show_column_names = T,
                     column_labels = names,
                     column_names_rot = 90,
                     column_names_gp = gpar(fontsize = 10),
                     show_heatmap_legend = FALSE)
  }
draw(HM.PISA)
for(i in 1:14) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(enrichment_terms[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left",
              gp = gpar(
                fontsize = 8, 
                fontfamily = "sans", 
                fontface = "bold"))
  })
}
