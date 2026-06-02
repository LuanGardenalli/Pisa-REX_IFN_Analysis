library(ComplexHeatmap)
library(circlize)

facet <- "REX"
proteins <- c("P42224", "P05161", "O14933", "Q5EBM0", "Q9Y3Z3", "P29590", "P19474", "O95786", "P84243")

#calculate Log2 Ratios and Cys positions:
log2ratios <- function(data.peptides,
                       data.proteins,
                       cell.line,
                       type,
                       proteins) {
  limma  <- pisaRex.pairwiseReport.accession(data.peptides,
                                             data.proteins,
                                             cell.line,
                                             type,
                                             facet="REX") %>% 
    filter(Accession %in% proteins)
  
  replicates <- pisaRex.preprocess(data.peptides,
                                   data.proteins,
                                   cell.line,
                                   type,
                                   facet="REX") %>% 
    filter(`Master Protein Accessions` %in% proteins) %>% 
    add_column(FDR = limma$adj.P.Val, P = limma$P.Value) %>% 
    mutate(start_pos = as.numeric(str_extract(`Positions in Master Proteins`, "(?<=\\[)\\d+(?=-)")),
           clean_seq = str_extract(`Annotated Sequence`, "(?<=\\.)[A-Z]+(?=\\.)"),
           `Cys Position` = map2_chr(clean_seq, start_pos, ~ {
             
             c_rel_pos <- str_locate_all(.x, "C")[[1]][, "start"]
             if (length(c_rel_pos) == 0) {
               return(NA_character_)
             }
             
             abs_pos <- .y + c_rel_pos - 1
             paste0("Cys", paste(abs_pos, collapse = "/"))})) %>%
    dplyr::select(-start_pos, -clean_seq) %>%
    mutate(`Cys Position` = paste0(`Master Protein Accessions`,": ", `Cys Position`))
  
  log2ratios <- replicates %>% 
    dplyr::select(`Master Protein Accessions`, `Cys Position`, FDR, P, contains("Abundance:")) %>% 
    arrange(FDR) %>% 
    distinct(`Cys Position`, .keep_all = T) %>% 
    rowwise() %>% 
    mutate(control_mean = rowMeans(pick(contains("Control")), na.rm = TRUE)) %>% 
    mutate(across(contains("Sample"), ~ control_mean - .x)) %>% 
    dplyr::select(!control_mean & !contains("Control")) %>% 
    mutate(`Cys Position` = str_remove(`Cys Position`, "^.*?: "))
  
  return(log2ratios)
}

#Input:
{
  cell.line <- "THP1"
  TA <- log2ratios(data.peptides,
                   data.proteins,
                   cell.line,
                   type="Alpha",
                   proteins)
  
  TB <- log2ratios(data.peptides,
                   data.proteins,
                   cell.line,
                   type="Beta",
                   proteins)
  
  TG <- log2ratios(data.peptides,
                   data.proteins,
                   cell.line,
                   type="Gamma",
                   proteins)
  
  cell.line <- "HL60"
  HA <- log2ratios(data.peptides,
                   data.proteins,
                   cell.line,
                   type="Alpha",
                   proteins)
  
  HB <- log2ratios(data.peptides,
                   data.proteins,
                   cell.line,
                   type="Beta",
                   proteins)
  
  HG <- log2ratios(data.peptides,
                   data.proteins,
                   cell.line,
                   type="Gamma",
                   proteins)
  
  heatmap.df.T <- full_join(TA, 
                            TB, 
                            by = join_by(`Master Protein Accessions`, `Cys Position`), 
                            keep = F,
                            suffix = c(".alpha", ".beta")) %>% 
    full_join(TG, 
              by = join_by(`Master Protein Accessions`, `Cys Position`), 
              keep = F) %>% 
    dplyr::rename(FDR.gamma = FDR,
           P.gamma = P)
  
  heatmap.df.H <- full_join(HA, 
                            HB, 
                            by = join_by(`Master Protein Accessions`, `Cys Position`), 
                            keep = F,
                            suffix = c(".alpha", ".beta")) %>% 
    full_join(HG, 
              by = join_by(`Master Protein Accessions`, `Cys Position`), 
              keep = F) %>% 
    dplyr::rename(FDR.gamma = FDR,
           P.gamma = P)
  
  heatmap.df <- full_join(heatmap.df.T,
                          heatmap.df.H,
                          by = join_by(`Master Protein Accessions`, `Cys Position`),
                          keep = F,
                          suffix = c(".THP1", ".HL"))
  
  #There was one repeated Cys position, which is fixed by setting a blank space in one of them:
  heatmap.df[9,2] <- " Cys189"
  
  #Creates mapping for FDR and P-value signs:
  FDRs <- heatmap.df %>% 
    dplyr::select(contains("FDR"))
  heatmap.df <- heatmap.df %>% dplyr::select(!contains("FDR"))
  
  Ps <- heatmap.df %>% 
    dplyr::select(contains("P."))
  heatmap.df <- heatmap.df %>% dplyr::select(!contains("P."))
  
  replicate.long <- heatmap.df %>%
    dplyr::select(!`Master Protein Accessions`) %>% 
    pivot_longer(-`Cys Position`,names_to = "sample", values_to = "log2fc") %>% 
    mutate(group = str_remove(sample, "^.*?, Sample, "))
  
  colnames(FDRs) <- unique(replicate.long$group)
  FDRs$`Cys Position` <- heatmap.df$`Cys Position`
  
  colnames(Ps) <- unique(replicate.long$group)
  Ps$`Cys Position` <- heatmap.df$`Cys Position`
  
  fdr.long <- FDRs %>% 
    pivot_longer(-`Cys Position`,names_to = "group", values_to = "fdr")
  
  p.long <- Ps %>% 
    pivot_longer(-`Cys Position`,names_to = "group", values_to = "p")
  
  fdr.matrix <- replicate.long %>%
    left_join(fdr.long, by = c("Cys Position", "group")) %>%
    dplyr::select(`Cys Position`, sample, fdr) %>%
    pivot_wider(names_from = sample, values_from = fdr) %>%
    column_to_rownames("Cys Position") %>%
    as.matrix()
  
  p.matrix <- replicate.long %>%
    left_join(p.long, by = c("Cys Position", "group")) %>%
    dplyr::select(`Cys Position`, sample, p) %>%
    pivot_wider(names_from = sample, values_from = p) %>%
    column_to_rownames("Cys Position") %>%
    as.matrix()
  
  
  heatmap.df[is.na(heatmap.df)] <- 0
  
  heatmap.matrix <- heatmap.df %>%
    dplyr::select(!`Master Protein Accessions`) %>% 
    column_to_rownames("Cys Position") %>%
    as.matrix()
  
  #Maps accessions to Symbols:
  gene.map <- as_tibble(clusterProfiler::bitr(heatmap.df$`Master Protein Accessions`, 
                             fromType = "UNIPROT",  
                             toType = c("SYMBOL"),    
                             OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  heatmap.df <- left_join(heatmap.df, gene.map, join_by(`Master Protein Accessions` == UNIPROT)) %>% 
    distinct(`Cys Position`, .keep_all = T)
}

#Draw HM:
{
  names <- c("", "THP1, Alpha", "",
             "", "THP1, Beta",  "",
             "", "THP1, Gamma", "",
             
             "", "HL60, Alpha", "",
             "", "HL60, Beta",  "",
             "", "HL60, Gamma", "")
  
  my_colors <- colorRamp2(c(-1, 0, 2), 
                          c("#2166AC", "#F7F7F7","#B2182B"))
  
  #Sets Cys positions that should be highlighted: 
  rows <- c(4,6,43,33,3,7,23,29,10,12,39,38,5,1,36,31,28,8,9,41,32,40,21,11,19,30,18,46,44)
  highlight <- heatmap.df %>% 
    dplyr::select(SYMBOL, `Cys Position`) %>% 
    rowid_to_column("row") %>% 
    mutate(color = if_else(row %in% rows, "red", "black"))
  
  
  ha = rowAnnotation(foo = anno_empty(border = FALSE, 
                                      width = max_text_width(proteins)))
  HM <- Heatmap(heatmap.matrix,
                name = "Log2\nRatio",
                col = my_colors,
                cluster_rows = T,
                cluster_columns = F,
                
                row_split = heatmap.df$SYMBOL,       
                row_title = "Oxidized Cysteine Positions",              
                show_row_names = T,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 8, col = highlight$color),
                right_annotation = ha,
                show_row_dend = F,
                
                show_column_names = TRUE,
                column_labels = names,
                column_names_rot = 35,
                column_names_gp = gpar(fontsize = 8),
                show_heatmap_legend = F,
                
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if (!is.na(fdr.matrix[i, j]) && fdr.matrix[i, j] < 0.05 && fdr.matrix[i, j] > 0.001) {
                    grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black"))
                  } else if (!is.na(fdr.matrix[i, j]) && fdr.matrix[i, j] < 0.001) {
                    grid.text("**", x, y, gp = gpar(fontsize = 10, col = "black"))
                  } else if (!is.na(p.matrix[i, j]) && p.matrix[i, j] < 0.05) {
                    grid.text("#", x, y, gp = gpar(fontsize = 7, col = "black"))
                  }
                })
  
  draw(HM)
  names <- names(row_order(HM))
  for(i in 1:9) {
    decorate_annotation("foo", slice = i, {
      grid.rect(x = 0, width = unit(1, "mm"), gp = gpar(fill = "black", col = NA), just = "right")
      grid.text(paste(names[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left",
                gp = gpar(
                  fontsize = 8, 
                  fontfamily = "sans", 
                  fontface = "bold"))})
  }
}
