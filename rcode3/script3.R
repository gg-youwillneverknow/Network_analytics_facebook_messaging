#LEGGIAMO IL FILE CON I GLI EDGES
filename="ia-fb-messages.txt"
# Read the text file as a data frame
data <- read.table(filename, header = FALSE, sep = " ")  # Replace with the path to your text file and appropriate settings
colnames(data)=c("Source","Target")
write.csv(data, file = "output.csv", row.names = FALSE)

#INIZIO ANALISI
library(igraph)
library(igraphdata)
library(Matrix)
library(ggplot2)
library(poweRlaw)

#CARICA E VISUALIZZA LA RETE 
net <- graph_from_data_frame(d = data, directed = FALSE)
# Get the number of nodes
N <- vcount(net)
N
# Get the number of edges
L <- ecount(net)
L
plot(net, edge.arrow.size=.2,vertex.label=NA)
is_connected(net)
sample_size=100
# Randomly sample nodes
sampled_nodes <- sample(V(net), size = sample_size)

# Subset the network with the sampled nodes
sampled_network <- induced_subgraph(net, sampled_nodes)

# Visualize the sampled network
plot(sampled_network,edge.arrow.size=.1)

#TROVA LA DIMENSIONE DELLA CLIQUE MASSIMA
dim_clique_max=clique_num(net)
dim_clique_max
#TROVA I NODI  DELLE CLIQUE MASSIMALI DEL GRAFO
cliques <- max_cliques(net)
cliques
# TROVA L'INSIEME DELLE CLIQUES MASSIMALI CHE PARTIZIONANO I NODI DEL GRAFO'
cliqueBP <- matrix(c(rep(paste0("cl", seq_along(cliques)), sapply(cliques, length)), names(unlist(cliques))), ncol=2, )
bp <- graph_from_edgelist(cliqueBP, directed = F)
V(bp)$type <- grepl("cl", V(bp)$name)
length(unique(cliqueBP[,1]))
plot(bp, layout=layout_as_bipartite)
#TROVA LE CLIQUE CON MINIMO 4 NODI 
cliques_dim=cliques(net, min=25)
#TROVA LA CLIQUE MASSIME
largest_cliques <- largest_cliques(net)
num_clique_max=length(largest_cliques)
max_cliques(net, min = NULL, max = NULL, subset = NULL, file = "cliques_network")
# verifico se è triangolato
is_chordal(net)

#DIAMETER
diam=diameter(net, directed=FALSE, weights=NA)
diam
diam_nodes=get_diameter(net, directed=FALSE)
diam_nodes

#DEGREES
deg <- degree(net, mode="all")
hist(deg, breaks=sqrt(length(deg)),main="Istogramma grado dei nodi", xlab = "Grado", ylab = "Frequenza")
deg.dist = degree_distribution(net, cumulative=FALSE, mode="all")
degrees=0:max(deg)
degree_frame = data.frame(degrees,deg.dist)
colnames(degree_frame)=c('degree','pk')
degree_frame=degree_frame[2:nrow(degree_frame),]
degree_frame <- degree_frame[degree_frame$pk != 0.000000e+00, ]
x=degree_frame$degree
y=degree_frame$pk
plot(x=x, y=y, pch=19, cex=1.2, col="orange",
     xlab="Grado k", ylab="P(k)", main="Distribuzione del grado dei nodi")


#SCALE-FREE
plot(x=log10(x), y=log10(y), pch=19, cex=1.2, col="orange",
      xlab="Grado", ylab="P(k)", main="log-log distribuzione del grado dei nodi")
abline(lm(log10(y) ~ log10(x), data = degree_frame), col = "blue")
# fit linear model
linear_model <- lm(log10(degree) ~ log10(pk), data=degree_frame)
summary(linear_model)
st_error=c(0.03)
estimate_alpha=c(-0.60)
p_value=c(2e-16)
R2=c(0.88)
summary_estimate=data.frame(estimate_alpha,st_error,p_value,R2)
write.csv(summary_estimate, file = "summary.csv", row.names = FALSE)
boxplot(linear_model[['residuals']],main='Boxplot: Residuals',ylab='residual value')
ggplot(data=degree_frame,aes(x = log10(x), y = log10(y))) +
  geom_point(colour = "orange") +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("Log(Grado)") + ylab("Log(Pk)") + 
  ggtitle("log-log distribuzion grado e intervallo di confidenza") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold",hjust = 0.5))

#Gilbert Random Network model 
k = 2 * L / N 
k=mean(deg)
p = k / (N-1)
poisson_dist <- function(gamma, lambda) {
  return((lambda^gamma) * exp(-lambda) / factorial(gamma))
}
# Add Poisson distribution function layer
poisson_values <- poisson_dist(x, k)   # Y-axis values for Poisson distribution
ggplot(degree_frame, aes(x = x, y = y)) + 
  geom_point(colour="orange") + 
  xlab("Grado") + ylab("Pk") + 
  ggtitle("Distribuzion grado e poisson") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold",hjust = 0.5)) +
  geom_line(aes(x = x, y = poisson_values), color = "red", linetype = "dashed")   # Customize the Poisson distribution

lnN=log(N, base = exp(1)) # since k is bigger than this value but greater than 1 we are in the connected regime
NG=components(net, mode = c("weak", "strong"))$csize[1]
ncc=components(net, mode = c("weak", "strong"))$no[1]
NG_random=N
avg_distance <- mean_distance(net)
random_avg_distance=log(N, base = exp(1))/log(k, base = exp(1)) # since the random_avg_distance is smaller than the avg_distance the network does not display the small world phenomena
clustering_coefficient <- transitivity(net, type = "local")
data2=data.frame(clustering_coefficient,deg)
colnames(data2)=c("clustering_coefficient","degree")
random_clustering = k/N
clustering_correlations=data.frame(degree=numeric(),clustering_coefficient=numeric())
for (k in unique(data2$degree)) {
  # Subset the dataframe for nodes with degree k
  degree_k_nodes_clustering <- data2[data2$degree == k,"clustering_coefficient"]
  average=mean(degree_k_nodes_clustering,na.rm = TRUE)
  new_row <- data.frame(degree = k, clustering_coefficient = average)
  clustering_correlations=rbind(clustering_correlations,new_row)
}
x=clustering_correlations$degree
y=clustering_correlations$clustering_coefficient
ggplot(data = clustering_correlations, aes(x = x, y = y)) +
  geom_point(color="black") + ggtitle("Coefficiente di clustering vs Grado") +
  xlab("k") + ylab("C(k)") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=0.05, color = "green")

#WALTZ STROGATZ
sw <- sample_smallworld(dim=2, size=N, nei=1, p=p)
plot(sw, vertex.size=6, vertex.label=NA, layout=layout_in_circle)
avg_clus_coef_sw=transitivity(sw, type="global")
avg_distance_sw=mean_distance(sw)

#CLOSENESS
closeness = closeness(net, mode="all", weights=NA)
hist(closeness, breaks=100, main="Histogram of node closeness")
betweenness = betweenness(net, directed=F, weights=NA)
hist(betweenness, breaks=100, main="Histogram of node betweenness")
edge_betweenness = edge_betweenness(net, directed=F, weights=NA)
hist(edge_betweenness, breaks=100, main="Histogram of edge betweenness")

avg_clus_coef=transitivity(net, type="global")
summary=data.frame(L=c(L),N=c(N),diameter=c(diam),k=c(k),avg_clust_coef=c(avg_clus_coef),
                   NG=c(NG),num_connected_components=c(ncc),dim_clique_max=c(dim_clique_max),
                   num_clique_max=c(num_clique_max))
write.csv(summary, file = "summary.csv", row.names = FALSE)

#COMMUNITIES
# proviamo con algoritmo divisivo edge betweenness
eb <- edge.betweenness.community(as.undirected(net))
# suppongo no = 3 comunità
membership <- cut_at(eb, no = 4)
#plot dendogram 
plot_dendrogram(eb)
plot(net, vertex.color= rainbow(3, .8, .8, alpha=.8))
sizes=sizes(eb)
# Calculate the probability of obtaining a community with each size
community_probabilities=prop.table(table(sizes))
# Create a dataframe with "Size" and "Probability" columns
communities_eb <- data.frame(Size = as.numeric(names(community_probabilities)),
                                 Probability = as.numeric(community_probabilities))

# proviamo con algoritmo divisivo random walk
wt <- walktrap.community(net)
# suppongo no = 4 comunità
membership <- cut_at(wt, no = 10)
#plot dendogram 
plot_dendrogram(wt)
plot(net, vertex.color= rainbow(3, .8, .8, alpha=.8)[membership])
sizes=sizes(wt)
# Calculate the probability of obtaining a community with each size
community_probabilities=prop.table(table(sizes))
# Create a dataframe with "Size" and "Probability" columns
communities_wt <- data.frame(Size = as.numeric(names(community_probabilities)),
                                 Probability = as.numeric(community_probabilities))

# proviamo con la modularità con il greedy method (hiearchical, fast method) 
c1 = cluster_fast_greedy(net) # modularity measure
modularity(c1)
B = modularity_matrix(net)
membership=membership(c1)
plot(net, vertex.color= rainbow(3, .8, .8, alpha=.8))
sizes=sizes(c1)
# Calculate the probability of obtaining a community with each size
community_probabilities=prop.table(table(sizes))
# Create a dataframe with "Size" and "Probability" columns
communities_greedy <- data.frame(Size = as.numeric(names(community_probabilities)),
                                 Probability = as.numeric(community_probabilities))
# Subset the network based on the nodes belonging to the specified community
community_nodes <- V(net)$name[membership == 4]
subgraph <- induced_subgraph(net, community_nodes)
layout <- layout.fruchterman.reingold(subgraph)
# Plot the subgraph representing the specific community
plot(subgraph, layout=layout, vertex.label = NA)

deg_sub=degree(subgraph)
# Define a color palette based on degree
color_palette <- colorRampPalette(c("blue", "red"))(max(deg_sub) + 1)
# Set the node color attribute based on degree
V(subgraph)$color <- color_palette[deg_sub + 1]
# Visualize the network with node colors based on degree
plot(subgraph, vertex.color = V(subgraph)$color)
 
#DEGREE CORRELATIONS 
adj_matrix <- get.adjacency(net)

# Compute the sum of degrees for each neighbour
neighbors_sum <- adj_matrix %*% data2$degree

# Compute the count of neighbors for each node
neighbor_counts <- rowSums(adj_matrix)

# Compute the average degree of neighbors for each node
avg_degrees <- neighbors_sum / neighbor_counts
data2$knn = as.vector(avg_degrees)
degree_correlations=data.frame(degree=numeric(),knn=numeric())
for (k in unique(data2$degree)) {
  # Subset the dataframe for nodes with degree k
  degree_k_nodes_knn <- data2[data2$degree == k,"knn"]
  average=mean(degree_k_nodes_knn)
  new_row <- data.frame(degree = k, knn = average)
  degree_correlations=rbind(degree_correlations,new_row)
}
x=degree_correlations$degree
y=degree_correlations$knn
ggplot(data = degree_correlations, aes(x = x, y = y)) +
  geom_point(color="black") + ggtitle("Degree correlations") +
  geom_smooth(method = "lm", se = TRUE) +
  xlab("k") + ylab("knn(k)") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

#ROBUSTNESS 
k2=sum(deg**2)/1266
molly_reed_criterion=k2/k
critical_threshold = 1 - (1 / (molly_reed_criterion - 1))
fER = 1 - (1/k)
robustness=data.frame(fER=c(fER),critical_threshold=c(critical_threshold))
write.csv(robustness, file = "robustness.csv", row.names = FALSE)


