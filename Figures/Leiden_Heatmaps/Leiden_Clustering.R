#kNN graph construction and Leiden clustering code: 

library(dbscan)
library(igraph)

#EXP:
res <- 0.3
K <- 20

knn <- kNN(heatmap.df[2:19], k = K)

edges <- data.frame(
  from = rep(1:nrow(heatmap.df[2:19]), each = K),
  to = as.vector(t(knn$id)))

knn_graph <- graph_from_data_frame(edges, directed = FALSE)
Leiden <- cluster_leiden(knn_graph, resolution = res)

clusters <- Leiden$membership
heatmap.df$Cluster <- as.factor(clusters)

#PISA:
res <- 0.1
K <- 100

knn <- kNN(heatmap.df[2:13], k = K)

edges <- data.frame(
  from = rep(1:nrow(heatmap.df[2:13]), each = K),
  to = as.vector(t(knn$id)))

knn_graph <- graph_from_data_frame(edges, directed = FALSE)
Leiden <- cluster_leiden(knn_graph, resolution = res)

clusters <- Leiden$membership
heatmap.df$Cluster <- as.factor(clusters)

#REX:
res <- 0.2
K <- 25

knn <- kNN(heatmap.df[2:19], k = K)

edges <- data.frame(
  from = rep(1:nrow(heatmap.df[2:19]), each = K),
  to = as.vector(t(knn$id)))

knn_graph <- graph_from_data_frame(edges, directed = FALSE)
Leiden <- cluster_leiden(knn_graph, resolution = res)

clusters <- Leiden$membership
heatmap.df$Cluster <- as.factor(clusters)
