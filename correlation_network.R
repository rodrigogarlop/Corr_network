# 2019-09-18 By Rodrigo García-López
# This script was tested with R 3.6.0 - "Planting of a Tree"
# It uses a contingency table (no initial commented lines are present so edit the input if required. Best if a reduced sparsity/transformed table is used.
# The input table should have samples as columns and features as rows
# The script depends on libraries corr and ggplot2

# Run as follows:
# cat table.tsv|~/bin/R-3.6.0/bin/Rscript correlation_network.R <output_prefix>
# Tested as: # cat 12_diversity_analyses/02_using_sparsity_reduced_tables/07_transformed_tables/01_Multiregion/06_Final_mix_gg_6/06_Final_mix_gg-lvl6-rar-10000-perm.tsv|~/bin/R-3.6.0/bin/Rscript correlation_network.R <output_prefix>
# df <- read.table("Filter_core_I/rarefied_6.6k-4y-I-lvl6_grp_I_.tsv", sep="\t",header=T, skip=0, fileEncoding = "UTF-8", comment.char='', check.names=FALSE)
# prefix="test"

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) { # at least, one arguments should be included: <out_name_prefix>
  stop("A minimum of 1 argument is mandatory: cat table.tsv|Rscript multi_distance_metrics.R <output_file>", call.=FALSE)
}
prefix <- as.character(args[1]) # Get a string handle to create output names
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, fileEncoding = "UTF-8", comment.char='', check.names=FALSE)
df <- rowsum(df[2:length(df)], group=df[[1]]) # Sum repeated rows
df <- as.matrix(df[order(rowSums(df),decreasing=T),])

### Correlation networks ###
# For a last analyses, we want to see if samples cluster by correlation (spearman)
library(corrr) # This requires the correlation package for r
library(ggplot2) # And one for the graphical output
library(corrplot)
subjects <- correlate(df, method = "spearman", use = "pairwise.complete.obs",diagonal = 1) #load a subject correlation object
write.table(fashion(stretch(subjects),leading_zeros = TRUE,decimals=10),paste(prefix,"subject_all_pair_cor.tsv", sep="-"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) #Export the whole subject pairwise comparison list
out <- fashion(subjects,leading_zeros = TRUE,decimals=4)
names(out) <- c("Columns",out[,1])
write.table(out,paste(prefix,"sq_col_cor_mat.tsv", sep="-"), sep="\t", quote=FALSE, row.names=FALSE) #Export the square matrix
pdf(paste(prefix,"sample_net_all_gt_0.3.pdf", sep="-"))
network_plot(subjects,min_cor = 0.3) #Only correlations over 0.70 will be drawn
dev.off()
pdf(paste(prefix,"sample_net_all_gt_0.7.pdf", sep="-"))
network_plot(subjects,min_cor = 0.7) #Only correlations over 0.70 will be drawn
dev.off()
pdf(paste(prefix,"sample_cor_heatmap_all.pdf", sep="-"))
p <-rplot(subjects) #Also create a correlation heatmap for subject/timepoints
p+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=6))
dev.off()
col1 <- colorRampPalette(c("red", "white", "blue"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")) # This is more similar to the default
pdf(paste(prefix,"corrplot_col.pdf", sep="-"))
# corrplot(cor(df,method="spearman"),method="color",type="upper",order="hclust",tl.cex=0.4,tl.col=1,addrect=1,tl.offset=0.2,hclust.method="complete",diag=F,col=col2(20))
corrplot(cor(df,method="spearman"),method="color",type="full",tl.cex=0.4,tl.col=1,addrect=1,tl.offset=0.2,diag=F,col=col2(20))
dev.off()
# "hclust.method" should be one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
if(ncol(df)>50){
	top50 <- df[1:50,]
	subjects <- correlate(top50, method = "spearman", use = "pairwise.complete.obs",diagonal = 1) #load a subject correlation object
	write.table(fashion(stretch(subjects),leading_zeros = TRUE,decimals=4),paste(prefix,"subject_top50_pair_cor.tsv", sep="-"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) #Export the subject matrix
	pdf(paste(prefix,"sample_net_top50_gt_0.3.pdf", sep="-"))
	network_plot(subjects,min_cor = 0.3) #Only correlations over 0.70 will be drawn
	dev.off()
	pdf(paste(prefix,"sample_net_top50_gt_0.7.pdf", sep="-"))
	network_plot(subjects,min_cor = 0.7) #Only correlations over 0.70 will be drawn
	dev.off()
	pdf(paste(prefix,"sample_cor_heatmap_top50.pdf", sep="-"))
	p <-rplot(subjects) #Also create a correlation heatmap for subject/timepoints
	p+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=6))
	dev.off()
}

#There are way too many otus to print them all, so we'll start by changing the labels to numbers
write.table(cbind("OTU-ID"=c(1:nrow(df)),"OTU Taxonomy"=row.names(df)),paste(prefix,"feature_list.tsv", sep="-"), sep="\t", quote=FALSE, row.names=F) #Export the list of otus
mat_otu <- t(df)
otus <- correlate(mat_otu, method = "spearman", use = "pairwise.complete.obs",diagonal = 1) # get correlation object for otus
write.table(fashion(stretch(otus),leading_zeros = TRUE,decimals=4),paste(prefix,"features_all_pair_cor.tsv", sep="-"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA) #Export the whole subject pairwise comparison list
out <- fashion(otus,leading_zeros = TRUE,decimals=4)
names(out) <- c("Rows",out[,1])
write.table(out,paste(prefix,"sq_row_cor_mat.tsv", sep="-"), sep="\t", quote=FALSE, row.names=FALSE) #Export the square matrix
# The last step produces a long list of all item pairwise correlations (all permutations) but cannot be used for plotting networks nor correlation heatmaps, thus we must reduce the total items
colnames(mat_otu) <- 1:nrow(df) #replace names with numbers as the complete taxonomies may be too long for plots
if(nrow(df)>100){
	mat_otu <- mat_otu[,1:100] # and keep only the top 100 most abundant
	otus <- correlate(mat_otu, method = "spearman", use = "pairwise.complete.obs",diagonal = 1) # get correlation object for otus
	pdf(paste(prefix,"feature_net_top100_gt_0.7.pdf", sep="-"))
	network_plot(otus,min_cor = 0.7) #Only correlations over 0.70 will be drawn
	dev.off()
}
if(nrow(df)>50){
	mat_otu <- mat_otu[,1:50] # now keep only the top 50 most abundant
	otus <- correlate(mat_otu, method = "spearman", use = "pairwise.complete.obs",diagonal = 1) # recalculate correlation object for otus
	pdf(paste(prefix,"feature_net_top50_gt_0.7.pdf", sep="-"))
	network_plot(otus,min_cor = 0.7) #Only correlations over 0.70 will be drawn
	dev.off()
	pdf(paste(prefix,"feature_cor_heatmap_top50.pdf", sep="-"))
	p <-rplot(otus) #Also create a correlation heatmap for subject/timepoints
	p+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),axis.text.y = element_text(size=6))
	dev.off()
}

pdf(paste(prefix,"corrplot_row.pdf", sep="-"))
# corrplot(cor(t(df),method="spearman"), method="color",type="upper",order="hclust",tl.cex=0.3,tl.col=1,addrect=1,tl.offset=0.2,hclust.method="complete",diag=F,col=col2(20))
corrplot(cor(t(df),method="spearman"),method="color",type="upper",tl.cex=0.3,tl.col=1,addrect=1,tl.offset=0.2,diag=F,col=col2(20))
dev.off()
