#EXP:
# input: heatmap.df[2:19] (Variable from the heatmapping scripts)
# res: 0.3
# K: 20

#PISA:
# input: heatmap.df[2:13]
# res: 0.1
# K: 100

#REX:
# input: heatmap.df[2:19]
# res: 0.2
# K: 25

Leiden_clustering <- function(res,
                              K,
                              input){
  library(dbscan)
  library(igraph)
  
  knn <- kNN(input, k = K)
  
  edges <- data.frame(
    from = rep(1:nrow(input), each = K),
    to = as.vector(t(knn$id)))
  
  knn_graph <- graph_from_data_frame(edges, directed = FALSE)
  Leiden <- cluster_leiden(knn_graph, resolution = res)
  
  clusters <- Leiden$membership
}
