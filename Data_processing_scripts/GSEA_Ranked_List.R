library(tidyverse)

#Creates a dataframe with t-values for each facet, for each pair:
{
  cell.line <- "THP1"
  type <- "Alpha"
  
  proteins.exp.TA <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="EXP") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.EXP=t)
  proteins.pisa.TA <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="PISA") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.PISA=t)
  proteins.rex.TA <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="REX") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.REX=t)
  
  input.TA <- merge(proteins.exp.TA, proteins.pisa.TA, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T)
  input.TA <- merge(input.TA, proteins.rex.TA, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T) %>% 
    mutate(Type = "Alpha")%>% 
    mutate(Cell.Line = "THP1")
  
} #THP-1, Alpha

{
  cell.line <- "THP1"
  type <- "Beta"
  
  proteins.exp.TB <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="EXP") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.EXP=t)
  proteins.pisa.TB <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="PISA") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.PISA=t)
  proteins.rex.TB <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="REX") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.REX=t)
  
  input.TB <- merge(proteins.exp.TB, proteins.pisa.TB, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T)
  input.TB <- merge(input.TB, proteins.rex.TB, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T) %>% 
    mutate(Type = "Beta")%>% 
    mutate(Cell.Line = "THP1")
} #THP-1,Beta

{
  cell.line <- "THP1"
  type <- "Gamma"
  
  proteins.exp.TG <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="EXP") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.EXP=t)
  proteins.pisa.TG <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="PISA") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.PISA=t)
  proteins.rex.TG <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="REX") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.REX=t)
  
  input.TG <- merge(proteins.exp.TG, proteins.pisa.TG, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T)
  input.TG <- merge(input.TG, proteins.rex.TG, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T) %>% 
    mutate(Type = "Gamma")%>% 
    mutate(Cell.Line = "THP1")
} #THP-1,Gamma

{
  cell.line <- "HL60"
  type <- "Alpha"
  
  proteins.exp.HA <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="EXP") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.EXP=t)
  proteins.pisa.HA <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="PISA") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.PISA=t)
  proteins.rex.HA <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="REX") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.REX=t)
  
  input.HA <- merge(proteins.exp.HA, proteins.pisa.HA, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T)
  input.HA <- merge(input.HA, proteins.rex.HA, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T) %>% 
    mutate(Type = "Alpha") %>% 
    mutate(Cell.Line = "HL60")
} #HL-60, Alpha

{
  cell.line <- "HL60"
  type <- "Beta"
  
  proteins.exp.HB <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="EXP") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.EXP=t)
  proteins.pisa.HB <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="PISA") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.PISA=t)
  proteins.rex.HB <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="REX") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.REX=t)
  
  input.HB <- merge(proteins.exp.HB, proteins.pisa.HB, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T)
  input.HB <- merge(input.HB, proteins.rex.HB, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T) %>% 
    mutate(Type = "Beta") %>% 
    mutate(Cell.Line = "HL60")
} #HL-60, Beta

{
  cell.line <- "HL60"
  type <- "Gamma"
  
  proteins.exp.HG <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="EXP") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.EXP=t)
  proteins.pisa.HG <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="PISA") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.PISA=t)
  proteins.rex.HG <- pisaRex.pairwiseReport(data.peptides, data.proteins, cell.line, type, facet="REX") %>% 
    arrange(abs(t)) %>% 
    distinct(Symbol, .keep_all = T) %>% 
    select(Symbol, t) %>% 
    rename(T.REX=t)
  
  input.HG <- merge(proteins.exp.HG, proteins.pisa.HG, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T)
  input.HG <- merge(input.HG, proteins.rex.HG, by.x="Symbol",by.y="Symbol", all.x=T, all.y=T) %>% 
    mutate(Type = "Gamma") %>% 
    mutate(Cell.Line = "HL60")
} #HL-60, Gamma

#Binds everything into an unified dataframe:
input.gsea <- bind_rows(input.TA, 
                        input.TB, 
                        input.TG, 
                        input.HA, 
                        input.HB, 
                        input.HG)


#Extract Pair of Interest:
line <- "THP1"
type <- "Alpha"
facet <-"EXP"
{
  GSEA <- input.gsea %>% 
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
