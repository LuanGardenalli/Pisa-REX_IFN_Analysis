library(tidyverse)
library(mixOmics)

#Formatting function:
format_for_diablo <- function(df) {
  df %>%
    pivot_longer(cols = -Accession, names_to = "Sample", values_to = "Intensity") %>%
    pivot_wider(names_from = Accession, values_from = Intensity) %>%
    column_to_rownames("Sample") %>%
    as.matrix()
}

#Extracting first two replicates of processed data:
diablo.input <- function(data.peptides,
                         data.proteins,
                         cell.line,
                         type) {
  
  names <- c("Accession", 
             paste0("Control1_",type), 
             paste0("Control2_",type), 
             paste0("Sample1_",type), 
             paste0("Sample2_",type))
  
  #Accounting for the inverted contrast of REX facet:
  names_rex <- c("Accession", 
                 paste0("Sample1_",type), 
                 paste0("Sample2_",type), 
                 paste0("Control1_",type), 
                 paste0("Control2_",type))
  
  exp <- pisaRex.preprocess(data.peptides,
                            data.proteins,
                            cell.line,
                            type,
                            facet="EXP")
  exp <- exp %>% dplyr::select(Accession, starts_with("Abundance"))
  exp <- exp %>% dplyr::select(1:3, 5:6)
  colnames(exp) <- names
  
  pisa <- pisaRex.preprocess(data.peptides,
                             data.proteins,
                             cell.line,
                             type,
                             facet="PISA")
  pisa <- pisa %>% dplyr::select(Accession, starts_with("Abundance"))
  colnames(pisa) <- names
  
  rex <- pisaRex.preprocess(data.peptides,
                            data.proteins,
                            cell.line,
                            type,
                            facet="REX")
  rex <- rex %>% dplyr::select(peptideRowID, starts_with("Abundance"))
  rex <- rex %>% dplyr::select(1:3, 5:6)
  colnames(rex) <- names_rex
  rex <- rex %>% dplyr::select(Accession,
                        paste0("Control1_",type), 
                        paste0("Control2_",type), 
                        paste0("Sample1_",type), 
                        paste0("Sample2_",type))
  return(list("exp"=exp, "pisa"=pisa, "rex"=rex))
}

#Run DIABLO on all IFN types:
cell.line <-  "THP1"
{
  #Input:
  type <- "Alpha"
  alpha <- diablo.input(data.peptides,
                        data.proteins,
                        cell.line,
                        type)
  type <- "Beta"
  beta <- diablo.input(data.peptides,
                       data.proteins,
                       cell.line,
                       type)
  type <- "Gamma"
  gamma <- diablo.input(data.peptides,
                        data.proteins,
                        cell.line,
                        type)
  
  X_exp_alpha <- format_for_diablo(alpha$exp)
  X_pisa_alpha <- format_for_diablo(alpha$pisa)
  X_rex_alpha <- format_for_diablo(alpha$rex)
  
  X_exp_beta <- format_for_diablo(beta$exp)
  X_pisa_beta <- format_for_diablo(beta$pisa)
  X_rex_beta <- format_for_diablo(beta$rex)
  
  X_exp_gamma <- format_for_diablo(gamma$exp)
  X_pisa_gamma <- format_for_diablo(gamma$pisa)
  X_rex_gamma <- format_for_diablo(gamma$rex)
  
  #Intersect data to only include proteins shared among all types:
  common_exp_features <- Reduce(intersect, list(
    colnames(X_exp_alpha), 
    colnames(X_exp_beta), 
    colnames(X_exp_gamma)))
  X_exp_alpha_aligned <- X_exp_alpha[, common_exp_features]
  X_exp_beta_aligned  <- X_exp_beta[, common_exp_features]
  X_exp_gamma_aligned <- X_exp_gamma[, common_exp_features]
  X_exp_combined <- rbind(X_exp_alpha_aligned, X_exp_beta_aligned, X_exp_gamma_aligned)
  
  common_pisa_features <- Reduce(intersect, list(
    colnames(X_pisa_alpha), 
    colnames(X_pisa_beta), 
    colnames(X_pisa_gamma)))
  X_pisa_alpha_aligned <- X_pisa_alpha[, common_pisa_features]
  X_pisa_beta_aligned  <- X_pisa_beta[, common_pisa_features]
  X_pisa_gamma_aligned <- X_pisa_gamma[, common_pisa_features]
  X_pisa_combined <- rbind(X_pisa_alpha_aligned, X_pisa_beta_aligned, X_pisa_gamma_aligned)
  
  common_rex_features <- Reduce(intersect, list(
    colnames(X_rex_alpha), 
    colnames(X_rex_beta), 
    colnames(X_rex_gamma)))
  X_rex_alpha_aligned <- X_rex_alpha[, common_rex_features]
  X_rex_beta_aligned  <- X_rex_beta[, common_rex_features]
  X_rex_gamma_aligned <- X_rex_gamma[, common_rex_features]
  X_rex_combined <- rbind(X_rex_alpha_aligned, X_rex_beta_aligned, X_rex_gamma_aligned)
  
  #Setting blocks and condition:
  X_list <- list(
    expression = X_exp_combined,
    solubility = X_pisa_combined,
    redox = X_rex_combined)
  
  Y_condition <- factor(c(
    "Control", "Control", "Alpha", "Alpha", 
    "Control", "Control", "Beta", "Beta",    
    "Control", "Control", "Gamma", "Gamma"))
  
  #Accounts for the experiments paired design:
  paired_design <- data.frame(
    sample_pair = factor(c(
      "Pair_A1", "Pair_A2", "Pair_A1", "Pair_A2",  
      "Pair_B1", "Pair_B2", "Pair_B1", "Pair_B2",  
      "Pair_C1", "Pair_C2", "Pair_C1", "Pair_C2" )))
  
  X_list_within <- list(
    expression = withinVariation(X = X_list$expression, design = paired_design),
    solubility = withinVariation(X = X_list$solubility, design = paired_design),
    redox      = withinVariation(X = X_list$redox,      design = paired_design))
  
  #Running DIABLO:
  design_matrix <- matrix(0.1, ncol = length(X_list), nrow = length(X_list), 
                          dimnames = list(names(X_list), names(X_list)))
  diag(design_matrix) <- 0
  
  keepX_grid <- c(10, 20, 30, 50, 75)
  
  tune_parameters <- list(
    solubility = keepX_grid,
    redox = keepX_grid,
    expression = keepX_grid)
  
  library(BiocParallel)
  parallel_setup <- SnowParam(workers = 6)
  
  diablo_tuning <- tune.block.splsda(
    X = X_list_within, 
    Y = Y_condition, 
    ncomp = 3, 
    test.keepX = tune_parameters, 
    design = design_matrix, 
    validation = 'loo', 
    dist = "centroids.dist", 
    BPPARAM = parallel_setup,
    progressBar = TRUE)
  optimal_keepX <- diablo_tuning$choice.keepX
  
  final_multiclass_diablo <- block.splsda(
    X = X_list_within, 
    Y = Y_condition, 
    ncomp = 3, 
    keepX = optimal_keepX, 
    design = design_matrix)
}
