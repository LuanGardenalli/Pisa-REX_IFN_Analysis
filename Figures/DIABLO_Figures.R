library(mixOmics)
library(patchwork)
library(uniqtag)
library(circlize)
library(ComplexHeatmap)

#Sample Projection Plots (3x1 format, C2-3):
{
  #Uses patchwork to plot DIABLO sample projection plots in a 3x1 grid:
  p_exp <- plotIndiv(final_multiclass_diablo, blocks = "expression", ind.names = F, 
                     group = Y_condition, comp = c(2, 3), legend = FALSE)$graph
  
  p_sol <- plotIndiv(final_multiclass_diablo, blocks = "solubility", ind.names = F, 
                     group = Y_condition, comp = c(2, 3), legend = FALSE)$graph

  p_red <- plotIndiv(final_multiclass_diablo, blocks = "redox", ind.names = F, 
                     group = Y_condition, comp = c(2, 3), legend = T)$graph
  
  final_1x3_plot <- p_exp | p_sol | p_red
  final_1x3_plot
}

#Heatmaps of discriminant proteins. DIABLO hits are filtered by requiring orthogonal significance from unadjusted P.val from the limma analysis.
cell.line <- "THP1"
c <- 1  #Set component of interest
targets <- selectVar(final_multiclass_diablo, comp = c)

#EXP/PISA:
facet <- "PISA"
P <- 0.05
adj.p <- 1
{
  #Filter DIABLO hits by limma significance:
  {
    rowID.A  <- pisaRex.pairwiseReport.accession(data.peptides, data.proteins, 
                                                 cell.line=cell.line,
                                                 type="Alpha",
                                                 facet=facet) %>% 
      filter(P.Value < P & adj.P.Val < adj.p & Accession %in% 
               if (facet=="EXP")targets$expression$name else 
                 if (facet=="PISA")targets$solubility$name) %>%
      dplyr::select(Accession) %>% pull()
    
    rowID.B  <- pisaRex.pairwiseReport.accession(data.peptides, data.proteins, 
                                                 cell.line=cell.line,
                                                 type="Beta",
                                                 facet=facet) %>% 
      filter(P.Value < P & adj.P.Val < adj.p & Accession %in% 
               if (facet=="EXP")targets$expression$name else 
                 if (facet=="PISA")targets$solubility$name) %>%
      dplyr::select(Accession) %>% pull()
    
    rowID.G  <- pisaRex.pairwiseReport.accession(data.peptides, data.proteins, 
                                                 cell.line=cell.line,
                                                 type="Gamma",
                                                 facet=facet) %>% 
      filter(P.Value < P & adj.P.Val < adj.p & Accession %in% 
               if (facet=="EXP")targets$expression$name else 
                 if (facet=="PISA")targets$solubility$name) %>%
      dplyr::select(Accession) %>% pull()
    
    accessions <- unique(c(rowID.A, rowID.B, rowID.G))
  }
  
  #Calculate Log2 Ratios:
  {
    A.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line=cell.line,
                                       type="Alpha",
                                       facet=facet) %>% 
      filter(Accession %in% accessions)
    A.log2ratios <- A.replicates %>% 
      dplyr::select(`Gene Symbol`, starts_with("Abundance")) %>% 
      rowwise() %>% 
      mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
      dplyr::select(`Gene Symbol`, contains("Sample"), control_mean) %>% 
      mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
      dplyr::select(!control_mean)
    
    #Some symbols can be repeated due to different isoforms. This handles it by numbering repeated names with parentheses:
    A.log2ratios$`Gene Symbol` <- paste0(make_unique(A.log2ratios$`Gene Symbol`, sep = " ("),")")
    A.log2ratios$`Gene Symbol` <- if_else(
      str_detect(A.log2ratios$`Gene Symbol`, "\\(\\d+\\)$"),
      A.log2ratios$`Gene Symbol`,
      str_remove(A.log2ratios$`Gene Symbol`, "\\)$"))
    
    B.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line=cell.line,
                                       type="Beta",
                                       facet=facet) %>% 
      filter(Accession %in% accessions)
    B.log2ratios <- B.replicates %>% 
      dplyr::select(`Gene Symbol`, starts_with("Abundance")) %>% 
      rowwise() %>% 
      mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
      dplyr::select(`Gene Symbol`, contains("Sample"), control_mean) %>% 
      mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
      dplyr::select(!control_mean)
    
    B.log2ratios$`Gene Symbol` <- paste0(make_unique(B.log2ratios$`Gene Symbol`, sep = " ("),")")
    B.log2ratios$`Gene Symbol` <- if_else(
      str_detect(B.log2ratios$`Gene Symbol`, "\\(\\d+\\)$"),
      B.log2ratios$`Gene Symbol`,
      str_remove(B.log2ratios$`Gene Symbol`, "\\)$"))
    
    G.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line=cell.line,
                                       type="Gamma",
                                       facet=facet) %>% 
      filter(Accession %in% accessions)
    G.log2ratios <- G.replicates %>% 
      dplyr::select(`Gene Symbol`, starts_with("Abundance")) %>% 
      rowwise() %>% 
      mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
      dplyr::select(`Gene Symbol`, contains("Sample"), control_mean) %>% 
      mutate(across(starts_with("Abundance:"), ~.x - control_mean)) %>% 
      dplyr::select(!control_mean)
    
    G.log2ratios$`Gene Symbol` <- paste0(make_unique(G.log2ratios$`Gene Symbol`, sep = " ("),")")
    G.log2ratios$`Gene Symbol` <- if_else(
      str_detect(G.log2ratios$`Gene Symbol`, "\\(\\d+\\)$"),
      G.log2ratios$`Gene Symbol`,
      str_remove(G.log2ratios$`Gene Symbol`, "\\)$"))
    
    heatmap.df <- merge(A.log2ratios, 
                        B.log2ratios, 
                        by= "Gene Symbol",
                        all.x=T, all.y=T)
    
    heatmap.df <- merge(heatmap.df, 
                        G.log2ratios, 
                        by= "Gene Symbol", 
                        all.x=T, all.y=T)
    
    heatmap.matrix <- heatmap.df %>%
      column_to_rownames("Gene Symbol") %>%
      as.matrix()
  }
  
  #Draw HM:
  {
    names <- c(rep("Alpha",3),
               rep("Beta",3),
               rep("Gamma",3))
    
    names_pisa <- c(rep("Alpha",2),
                    rep("Beta",2),
                    rep("Gamma",2))
    
    my_colors <- colorRamp2(c(-1, 0, 1), 
                            c("#2166AC", "#F7F7F7", "#B2182B"))
    
    Heatmap(heatmap.matrix,
            name = paste0("Log2 Ratio:\n(s)PLS-DA C",c,"\n",facet),
            col = my_colors,
            cluster_rows = T,
            cluster_columns = F,
            
            row_title = NULL,              
            show_row_names = T,        
            
            show_column_names = TRUE,
            column_labels = if (facet=="EXP") names else
              if (facet=="PISA") names_pisa,
            column_names_rot = 90,
            column_names_gp = gpar(fontsize = 6),
            row_names_gp = gpar(fontsize = 7),
            
            heatmap_legend_param = list(
              legend_height = unit(2, "cm"),   
              legend_width = unit(2, "cm"),    
              labels_gp = gpar(fontsize = 5), 
              title_gp = gpar(fontsize = 7)))
  }
}

#REX:
P <- 0.05
adj.p <- 1
{
  #Filter DIABLO hits by limma significance:
  {
    facet <- "REX"
    rowID.A  <- pisaRex.pairwiseReport.accession(data.peptides, data.proteins, 
                                                 cell.line=cell.line,
                                                 type="Alpha",
                                                 facet=facet) %>% 
      filter(P.Value < P & adj.P.Val < adj.p & peptideRowID %in% targets$redox$name) %>%
      dplyr::select(peptideRowID) %>% pull()
    
    rowID.B  <- pisaRex.pairwiseReport.accession(data.peptides, data.proteins, 
                                                 cell.line=cell.line,
                                                 type="Beta",
                                                 facet=facet) %>% 
      filter(P.Value < P & adj.P.Val < adj.p & peptideRowID %in% targets$redox$name) %>%
      dplyr::select(peptideRowID) %>% pull()
    
    rowID.G  <- pisaRex.pairwiseReport.accession(data.peptides, data.proteins, 
                                                 cell.line=cell.line,
                                                 type="Gamma",
                                                 facet=facet) %>% 
      filter(P.Value < P & adj.P.Val < adj.p & peptideRowID %in% targets$redox$name) %>%
      dplyr::select(peptideRowID) %>% pull()
    
    accessions <- unique(c(rowID.A, rowID.B, rowID.G))
  }
  
  #Calculate Log2 Ratios:
  {
    A.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line=cell.line,
                                       type="Alpha",
                                       facet=facet) %>% 
      filter(peptideRowID %in% accessions)
    
    A.log2ratios <- A.replicates %>% 
      dplyr::select(peptideRowID,`Master Protein Accessions`, starts_with("Abundance")) %>% 
      rowwise() %>% 
      mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
      dplyr::select(peptideRowID, `Master Protein Accessions`,  contains("Sample"), control_mean) %>% 
      mutate(across(starts_with("Abundance:"), ~ control_mean - .x)) %>% 
      dplyr::select(!control_mean)
    
    B.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line=cell.line,
                                       type="Beta",
                                       facet=facet)  %>% 
      filter(peptideRowID %in% accessions)
    
    B.log2ratios <- B.replicates %>% 
      dplyr::select(peptideRowID,`Master Protein Accessions`, starts_with("Abundance")) %>% 
      rowwise() %>% 
      mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
      dplyr::select(peptideRowID, `Master Protein Accessions`,  contains("Sample"), control_mean) %>% 
      mutate(across(starts_with("Abundance:"), ~ control_mean - .x)) %>% 
      dplyr::select(!control_mean)
    
    G.replicates <- pisaRex.preprocess(data.peptides,
                                       data.proteins,
                                       cell.line=cell.line,
                                       type="Gamma",
                                       facet=facet) %>% 
      filter(peptideRowID %in% accessions)
    
    G.log2ratios <- G.replicates %>% 
      dplyr::select(peptideRowID, `Master Protein Accessions`,  starts_with("Abundance")) %>% 
      rowwise() %>% 
      mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
      dplyr::select(peptideRowID, `Master Protein Accessions`,  contains("Sample"), control_mean) %>% 
      mutate(across(starts_with("Abundance:"), ~ control_mean - .x)) %>% 
      dplyr::select(!control_mean)
    
    heatmap.df <- merge(A.log2ratios, 
                        B.log2ratios, 
                        by= c("peptideRowID", "Master Protein Accessions"),
                        all.x=T, all.y=T)
    
    heatmap.df <- merge(heatmap.df, 
                        G.log2ratios, 
                        by= c("peptideRowID", "Master Protein Accessions"),
                        all.x=T, all.y=T)
    heatmap.df[is.na(heatmap.df)] <- 0
    
    
    #Some UniProt IDs fail to map to Gene Symbol by cluster profiler, due to alternative isoform identifiers. This code manually replaces such occurrences with the main protein accession:
    heatmap.df <- heatmap.df %>% 
      mutate(`Master Protein Accessions` = str_replace_all(`Master Protein Accessions`, fixed("A0A140T9T7"), "Q03518"),
             `Master Protein Accessions` = str_replace_all(`Master Protein Accessions`, fixed("A0A7I2V5Y7"), "Q15233"),
             `Master Protein Accessions` = str_replace_all(`Master Protein Accessions`, fixed("A0A0D9SGE8"), "Q8IWS0"),
             `Master Protein Accessions` = str_replace_all(`Master Protein Accessions`, fixed("A0A286YFM8"), "O43175"),
             `Master Protein Accessions` = str_replace_all(`Master Protein Accessions`, fixed("A0A2R8Y793"), "P60709"),
             `Master Protein Accessions` = str_replace_all(`Master Protein Accessions`, fixed("A0A804HJZ5"), "O60907"),
             `Master Protein Accessions` = str_replace_all(`Master Protein Accessions`, fixed("K7EJJ5"), "Q00610"))
    
    #Mapping accessions to symbols:
    heatmap.df$`Master Protein Accessions` <- sub("-.*", "", heatmap.df$`Master Protein Accessions`)
    gene.map <- as_tibble(clusterProfiler::bitr(heatmap.df$`Master Protein Accessions`, 
                               fromType = "UNIPROT",  
                               toType = c("SYMBOL"),    
                               OrgDb = org.Hs.eg.db::org.Hs.eg.db))
    heatmap.df <- merge(heatmap.df, 
                        gene.map, 
                        by.x ="Master Protein Accessions", 
                        by.y = "UNIPROT", 
                        all.x = TRUE)
    heatmap.df[is.na(heatmap.df)] <- 0
    
    heatmap.df <- heatmap.df %>% 
      dplyr::select(!peptideRowID & !`Master Protein Accessions`)
    
    #Adresses potential repeated symbol names in case there are peptide hits belonging to the same protein: 
    heatmap.df$SYMBOL <- paste0(make_unique(heatmap.df$SYMBOL, sep = " ("),")")
    heatmap.df$SYMBOL <- if_else(
      str_detect(heatmap.df$SYMBOL, "\\(\\d+\\)$"),
      heatmap.df$SYMBOL,
      str_remove(heatmap.df$SYMBOL, "\\)$"))
    
    heatmap.matrix <- heatmap.df %>%
      column_to_rownames("SYMBOL") %>%
      as.matrix()
  }
  
  #Draw HM:
  {
    names <- c(rep("Alpha",3),
               rep("Beta",3),
               rep("Gamma",3))
    
    my_colors <- colorRamp2(c(-1, 0, 1), 
                            c( "#2166AC", "grey","#B2182B"))
    
    Heatmap(heatmap.matrix,
            name = paste0("Log2 Ratio:\n(s)PLS-DA C",c,"\n",facet),
            col = my_colors,
            cluster_rows = T,
            cluster_columns = F,
            
            row_title = NULL,              
            show_row_names = T,        
            
            show_column_names = TRUE,
            column_labels = names,
            column_names_rot = 90,
            column_names_gp = gpar(fontsize = 7),
            row_names_gp = gpar(fontsize = 6),
            
            heatmap_legend_param = list(
              legend_height = unit(2, "cm"),
              legend_width = unit(2, "cm"),
              labels_gp = gpar(fontsize = 5), 
              title_gp = gpar(fontsize = 7)))
  }
}
