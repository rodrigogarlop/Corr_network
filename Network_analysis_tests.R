# Chord plot
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
# TEST 1 # General adjacency matrix example
	# Library
	library(igraph)
	# Create data
	set.seed(1)
	data <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.8,0.2)), nc=10)
	network <- graph_from_adjacency_matrix(data , mode='undirected', diag=F )
	# Default network
	par(mar=c(0,0,0,0))
	plot(network)
	
# TEST 2 # two-column table example
# create data:
links <- data.frame(
    source=c("A","A", "A", "A", "A","F", "B"),
    target=c("B","B", "C", "D", "F","A","E")
    )
# TEST 3 # literal input
network <- graph_from_literal(A-B-C-D-A,A-E,A-C);plot(network) 
network <- graph_from_literal(A-B-C-D,E-A-E-A,D-C-A,D-A-D-C,D-C, simplify=FALSE);plot(network) # In this other version, edges can be repeated
#this is not the same, as edges are weighted differently
# create the network object
network <- graph_from_data_frame(d=links, directed=F) 

# plot it
plot(network)

library(igraph)
setwd("/home/rod/Documents/01_Projects/16S")
name="3yno2017-I-lvlASVgt00_Avg"
rho_Gr=0.6 # set group rho
rho_Ft=0.8 # set feature rho
min_freq=0.001 # set corresponding abs or rel min frequency
freq=0.01
df <- read.table(paste0("Filter_core_I/grp_collated_3yno2017/",name,"_grp_stats.tsv"), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, check.names = FALSE, row.names=1) # load the average for class taxonomy level of the 3 year core with no 2017 samples (this includes the average for each grouping, such as whole, by year, by pond and (year/pond). IMPORTANT: this is normally not a square adjeacency matrix but rather a contingency (incidence) matrix (e.g. ASVs vs group matrix)
# dim(df)
# [1] 36 16
df[df<min_freq]=0 # First, ignore those observations not accounting for at least 0.1% (0.001 prob). This is done so that very small values do not influence the result
df <- df[,colSums(df)>0]# Then, purge empty columns and rows
df <- df[rowSums(df)>0,]
# dim(df)
# [1] 23 16
# With this incidence (contingency) matrix we can create a network relating rows and columns:
pdf(paste0("net-cross_",name,"Freq",freq,".pdf"))
network <- graph_from_incidence_matrix(df, multiple=FALSE,weighted=TRUE)
plot(network, vertex.label.family="Helvetica", vertex.label.color="black", edge.color="lightgray", edge.lty=1, vertex.color = c(rep("chartreuse3",nrow(df)),rep("coral1",ncol(df))), vertex.frame.color = c(rep("forestgreen",nrow(df)),rep("orangered",ncol(df))),vertex.label.cex=0.5, vertex.size=5, main=paste0("Feats vs groups - Set: ",name," Freq cutoff ",min_freq))
df2=df
df2[df2<freq]=0 # trim those with freq lower than this other set
df2 <- df2[,colSums(df2)>0]# Then, purge empty columns and rows
df2 <- df2[rowSums(df2)>0,]
network <- graph_from_incidence_matrix(df2, multiple=FALSE,weighted=TRUE)
plot(network, vertex.label.family="Helvetica", vertex.label.color="black", edge.color="lightgray", edge.lty=1, vertex.color = c(rep("chartreuse3",nrow(df2)),rep("coral1",ncol(df2))), vertex.frame.color = c(rep("forestgreen",nrow(df2)),rep("orangered",ncol(df2))),vertex.label.cex=0.5, vertex.size=5,main=paste0("Feats vs groups - Set: ",name," Freq cutoff ",freq))
dev.off()
cor_mat_c <-cor(df, method="spearman") # get a spearman correlation matrixby group, (non-parametric, which accounts to the lack of normal data in the sets) which is expected as they follow a non-normal law distro. IMPORTANT: This is the square adjacency matrix that will be used for the actual network (and should be symmetrical)
# net <- graph_from_adjacency_matrix(cor_mat_c, weighted=TRUE, mode='upper', diag=F)
# plot(net) # All non-zero connections were made since it is highly connected (groups produce a very similar "means" vector), we need to reduce total edges. We will keep only those between 0.8 and 1 (no negative correlations are present
cor_mat_c[cor_mat_c<rho_Gr]=0
net <- graph_from_adjacency_matrix(cor_mat_c, weighted=TRUE, mode='upper', diag=F)
pdf(paste0("net-groups_",name,"-rho-",rho_Gr,".pdf"))
par(mar=c(1,1,1,1))
plot(net, main=paste0("Groups - Set: ",name," Rho: ",rho_Gr), vertex.label.family="Helvetica", vertex.label.color="black", edge.color="lightgray", edge.lty=2, vertex.color = "cornflowerblue", vertex.frame.color = "blue3", vertex.label.cex=0.5, layout=layout.fruchterman.reingold)
dev.off()
# Now create the same type of network for the actual taxa:
cor_mat_r <-cor(t(df), method="spearman") # this has positive and negative correlations
# net <- graph_from_adjacency_matrix(cor_mat_r, weighted=TRUE, mode='upper', diag=F)
# plot(net) # too many correlations are non-zero
cor_mat_r[(cor_mat_r>(rho_Ft*(-1)))&(cor_mat_r<rho_Ft)]=0 # set threshold for usefull rho_Ft values
net <- graph_from_adjacency_matrix(cor_mat_r, weighted=TRUE, mode='upper', diag=F)
pdf(paste0("net-feats_",name,"-rho-",rho_Ft,".pdf"))
par(mar=c(1,1,1,1))
plot(net, main=paste0("Features - Set: ",name," Rho: ",rho_Ft), vertex.label.family="Helvetica", vertex.label.color="black", edge.color="lightgray", edge.lty=1, vertex.color = "firebrick", vertex.frame.color = "red3",vertex.label.cex=0.5, vertex.size=3, layout=layout.circle)
dev.off()

library(reshape2)
long <- melt(cor_mat_r) # use long format for edges #This creates double items, one for A<-B and one for B->A
# Maybe include a step to clean the upper or lower triangles prior to melting
# dim(long)
# [1] 26896     3
long <- long[long[,3]!=0,] # Remove empty edges (0s in the table)
# dim(long)
# [1] 590   3
range <- seq(-1,1,0.01) # Create a range for expected rho values (201 items, from -1 to 1 by default)
dict <- colorRampPalette(c("red", "white", "blue"))(length(range)) # and create a dictionary for colors
names(dict) <- range

# library(RColorBrewer)
# brewer.pal(3,"Reds")

# Now, for customizing the plot, we create a metadata dataframe:
bin_df <- df; bin_df[bin_df>0] = 1 #Copy and create a binary table from the original data
range_0_2 <- function(x){(x-min(x))/(max(x)-min(x))*2}
dist <- range01(rowSums(bin_df))+5 # Get a range between 0.5 and 1.5 to represent how many samples have each feature
pos <- 1:nrow(cor_mat_r)
names <- rownames(cor_mat_r)
names <- as.numeric(gsub("ASV_","",rownames(cor_mat_r))) # This is for ASVs only since sorting will result in a bad order if done alfabetically
ord <- order(names)
meta <- data.frame(pos, names, ord, dist) # Create a metadata dataframe
meta <- meta[ord,] # now, sort it alphabetically
col <- rep(topo.colors(trunc(nrow(meta)/10)+1),each=10)[1:nrow(meta)]
meta <- cbind(meta,col) #Add colors to the dataframe
meta <- meta[meta$pos,] # finally, revert to the original order
pdf(paste0("net-feats_",name,"-rho-",rho_Ft,"color.pdf"))
plot(net, main=paste0("Features - Set: ",name," Rho: ",rho_Ft), vertex.label.family="Helvetica", vertex.label.color="black", edge.color=dict[as.character(round(long[,3],2))], edge.lty=1, vertex.color = meta$col, vertex.frame.color = "white",vertex.label.cex=0.5, layout=layout.circle, edge.width=abs(long[,3]*1.5), vertex.size=meta$dist)
dev.off()

# Measurements
# Centrality
deg <- degree(net,mode="ALL") # in or out degree (either in this case, as it's undirected)
eig <- evcent(net)$vector #eigenvector (based on neighboring nodes)
# bet <- betweenness(net, directed=FALSE) #This won't accept negative values for edge weight
net2 <- net # instead, we copy the network
E(net2)$weight <- abs(E(net2)$weight) # use absolute values intestead (correlations range from -1 to 1) so this should not affect the overal relationships
bet <- betweenness(net2, directed=FALSE) # now we can measure betweenness
# Network Structure
edge_density(net) # whole network edge_density
# [1] 0.01593596
target1 <- met[met[,6]=="#0019FFFF",1] # Name of vertices from Group 2 as defined by color in met 
sub1 <- induced_subgraph(net,target1, impl="auto") #This creates a subset. The implementation parameter allows for create_from_scratch or copy_and_delete (the last is to split selected nodes then delete those edges with missing nodes
edge_density(sub1)
# [1] 0.02222222 # higher density but because there are only two nodes
target2 <- met[met[,6]=="#4C00FFFF",1]
sub2 <- induced_subgraph(net,target2, impl="auto")
edge_density(sub2)
# [1] 0.04444444
# for assortativity we require a factor that matches the nodes, we'll use met[,6]
assortativity_nominal(net, types=met[,6])
# [1] 0.07760395 # We do not know if it is high or low so we should compare to randomized items
rand_data <- sapply(1:1000, function(x){assortativity_nominal(net, types=sample(met[,6]))}) # create a randomized data by sampling our existing groups and shuffling them 1k times to calculate assortativity
hist(rand_data,xlim=c(min(rand_data),0.08))
abline(v=assortativity_nominal(net, types=met[,6]))
# as_long_data_frame(net) # save as long format
# Other usefull plotting parameters (use with plot())
# vertex.frame.color="transparent"
# Alternative version (load as edge list + metadata)
# long # This object has the long format for edges and weights as a df
# meta # This object has the metadata for nodes as df
met <- cbind(rownames(meta),meta) # meta was missing the names in the first column, we added those. The first row in both must match with node names (all in edges included in node list)
graph_from_data_frame(long, directed=FALSE, vertices = met)
n <- graph_from_data_frame(long, directed=FALSE, vertices = met, weighted=TRUE)

# checkpoint
save.image("chkp1.RData")
load("chkp1.RData")

# clusters
clusters <- cluster_louvain(net2) # This requires only positive values for edge weights
pdf("test.pdf");plot(clusters, net2,vertex.size=sqrt(deg),vertex.label=NA,edge.width=E(net2)$weight);dev.off()

# OTHER IMPORTANT FUNCTIONS (examples)
plot(sample_gnm(n=100, m=40), vertex.size=6, vertex.label=NA) # Random graph
rewire(rn, each_edge(prob=0.1)) # randomly rewire net rn
plot(rn %du% tr, vertex.size=10, vertex.label=NA)  # merge two graph plots
simplify(net, remove.multiple = F, remove.loops = T) # remove multiple and loops


# TEST cluster_edge_betweenness
par(mfrow = c(2,3))

lo <- layout_with_fr(in_net) 
clu_edg_btw <- cluster_edge_betweenness(in_net, directed=F)
len <- length(clu_edg_btw1)
membership(clu_edg_btw1)
modularity(clu_edg_btw1)

# clu_edg_btw0.1 <- cluster_edge_betweenness(in_net, directed=F, weights=temp)
# len <- length(clu_edg_btw0.1)
# membership(clu_edg_btw0.1)
# modularity(clu_edg_btw0.1)


temp = E(in_net)$scaled_01+0.001
clu_edg_btw1 <- cluster_edge_betweenness(in_net, directed=F, weights=temp)
len <- length(clu_edg_btw1)
membership(clu_edg_btw1)
modularity(clu_edg_btw1)

rev = 1.01-E(in_net)$scaled_01
clu_edg_btw2 <- cluster_edge_betweenness(in_net, directed=F, weights=rev)
len <- length(clu_edg_btw2)
membership(clu_edg_btw2)
modularity(clu_edg_btw2)

dendPlot(clu_edg_btw, mode="hclust")
dendPlot(clu_edg_btw1, mode="hclust")
dendPlot(clu_edg_btw2, mode="hclust")

plot(clu_edg_btw, in_net, layout=lo)
plot(clu_edg_btw1, in_net, layout=lo)
plot(clu_edg_btw2, in_net, layout=lo)

par(mfrow = c(1,1))
