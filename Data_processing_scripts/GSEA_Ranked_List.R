library(tidyverse)

#Creates a dataframe with t-values for each facet, for each pair:
mergeFacets <- function(data.peptides, 
                        data.proteins, 
                        cell.line, 
                        type){
  proteins.exp <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="EXP") %>% 
    arrange(desc(abs(t))) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.EXP=t)
  proteins.pisa <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="PISA") %>% 
    arrange(desc(abs(t))) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.PISA=t)
  proteins.rex <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="REX") %>% 
    arrange(desc(abs(t))) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.REX=t)
  
  input <- merge(proteins.exp, proteins.pisa, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T)
  input <- merge(input, proteins.rex, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T) %>% 
    mutate(Type = type)%>% 
    mutate(Cell.Line = cell.line)
}
{
  cell.line <- "THP1"
  TA <- mergeFacets(data.peptides, 
                    data.proteins, 
                    cell.line, 
                    type="Alpha")
  
  TB <- mergeFacets(data.peptides, 
                    data.proteins, 
                    cell.line, 
                    type="Beta")
  
  TG <- mergeFacets(data.peptides, 
                    data.proteins, 
                    cell.line, 
                    type="Gamma")
  
  cell.line <- "HL60"
  HA <- mergeFacets(data.peptides, 
                    data.proteins, 
                    cell.line, 
                    type="Alpha")
  
  HB <- mergeFacets(data.peptides, 
                    data.proteins, 
                    cell.line, 
                    type="Beta")
  
  HG <- mergeFacets(data.peptides, 
                    data.proteins, 
                    cell.line, 
                    type="Gamma")
}

#Binds everything into an unified dataframe:
input.gsea.combined <- bind_rows(TA, 
                                 TB, 
                                 TG, 
                                 HA, 
                                 HB, 
                                 HG)

#Extract Pair of Interest:
line <- "THP1"
type <- "Alpha"
facet <-"EXP"
{
  GSEA <- input.gsea.combined %>% 
    filter(Cell.Line==line,
           Type==type) %>% 
    select(Symbol, 
           if (facet=="EXP") {"T.EXP"} 
           else if (facet=="PISA") {"T.PISA"} 
           else if (facet=="REX") {"T.REX"} ) %>% 
    arrange(desc(if (facet=="EXP") {T.EXP} 
                 else if (facet=="PISA") {T.PISA} 
                 else if (facet=="REX") {T.REX} )) %>% 
    drop_na()
}
