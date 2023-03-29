# Started on 2020-10-08
# by Rodrigo García-López for Adrian Ochoa-Leyva's Metagenomics and Metatranscriptomics Lab at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# Disclamer: This script was written for use with specific data and are not therefore warrantied to be usable with different sets. They were created to carry out particular analyses and were originally intended for the lab's requirements, not commercial nor distribution.
# This script was tested with R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
# The script is intended to calculate between and within differences of the groups defined as parameters, and differences by each of the groups defined.

# The input is a square simmetrical correlation matrix with values ranging from -1 to 1

# Run as follows:
# cat dm.tsv|Rscript Compare_correlations_per_grp.R <prefix_output>  <metric_title> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>

# Tested with command:
# cat /3y-I-lvl6-sq_col_cor_mat.tsv|Rscript PCoA_from_dm.R //3y-I-lvl6-sample 3y-I-lvl6 2015_Z1 2016_Z2 2017_Z3 2018_Z4 D W IH
# Test in R:
# df <- read.table("/home/rod/Documents/01_Projects/16S/Health_Disease/Filter_core_I/correlations/3y-I-lvl6-sq_col_cor_mat.tsv", sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F)
# prefix="test"
# metric="3y_core_lvl6" # The name of the metric in use, for the graph titles
# # groups <- c("Z1","Z2","Z3","Z4","Z5")
# # groups <- c("2015","2016","2017","2018")
# # groups <- c("H","I","S")
# # groups <- c("HE","IE","HL","IL","HM","IM")
# # groups <- c("E","L","M")
# groups <- c("2015_Z1_H","2015_Z1_I","2015_Z2_H","2015_Z2_I","2016_Z2_H","2016_Z2_I","2015_Z3_I","2015_Z3_S","2016_Z3_S","2017_Z3_H","2017_Z3_I","2017_Z3_S","2018_Z4_I","2018_Z5_H")
# # groups <-("2015_Z1_I", "2015_Z2_H", "2016_Z2_H", "2016_Z2_I", "2017_Z3_H", "2017_Z3_I", "2018_Z4_I", "2018_Z5_H")
# groups <- c("2015_Z1", "2016_Z2", "2017_Z3", "2018_Z4", "D", "W", "IH")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) { # at least, two arguments should be included: <min_nonNAs> <prefix_output>  <name_of_alpha_metric>
  stop("A minimum of 4 arguments are mandatory: cat table.tsv|Rscript PCoA_from_dm.R <prefix_output>  <metric_title> <#_group_name_1> <#_group_name_2> ... <#_group_name_n>", call.=FALSE)
}
prefix <- as.character(args[1]) # Get a string handle to create output names
metric <- as.character(args[2])  # and the name of the metric (free strings)
groups <- as.character(args[3:length(args)]) # Create a vector with the name of the groups (this should be included in the sample names and be exclusive to each group for this to work)
df <- read.table(file('stdin'), sep="\t",header=T, skip=0, comment.char='',quote="",fill=F, row.names=1, check.names=F) # 
print("Loaded parameters:")
print(paste("Prefix:",prefix))
print(paste("Metric name:",metric))
print("Groups that should be present:")
print(groups)

####################################### Pre-load #######################################
 ### Define colours ### 
# Colors will only be used if the total number of groups is the same or less than the number of colors
# pair or color sets used for different comparisons
# group_cols <- c("coral1","red","cornflowerblue","blue") # I've used this for random 5 group items
# group_cols <- c("black","red","green","blue","magenta","orange","coral1","yellow") # I've used this for runs (5 grps)
# group_cols <- c("coral1","cornflowerblue","turquoise1") # I've used these for ponds (3 grps)
# group_cols <- c("springgreen4","darkorchid3","tan1") # I've used these for organs (3 grps)
# group_cols <- c("navyblue","hotpink","forestgreen","orangered") # I've used these for years (4 grps)
# group_cols <- c("dodgerblue4","plum","greenyellow","tan1") # Alternative for years (DEPRECATED)
# group_cols <- c("chartreuse2","darkorchid3","darkgreen","magenta","seagreen3","violetred") # I've used these for 2organs/3ponds (6 gprs)
# group_cols <- c("chartreuse2","orchid","limegreen","orchid","darkseagreen1","pink") #alternative for organ/ponds (NOW DEPRECATED)
group_cols <- c("black","red","green","blue","magenta","orange","coral1","yellow") # I've used this for runs (5 grps)
# group_cols <- c("navyblue","cornflowerblue","royalblue4","darkorchid","orange","coral1","yellow") # I've used this for years/runs/organs
# group_cols <- c("red4","navyblue","orangered4","darkorchid4","red1","cornflowerblue","darkslateblue","black","gray30","tan4","royalblue4","darkslategray","darkorchid","gold2") # alternative version for years/runs/organs (DEPRECATED)
# group_cols <- rep("gray",length(group_cols))
 ### Define Functions ###
# Predefine some important global objects
prepare_data <- function(df) { # create presets (sample totals and groups # MOD: ignored empty groups
	samples <- names(df) # extract names
	grp_sam <- sapply(groups,function(x) grep(x,samples)) # create index of samples in each group
	grp_less <- sapply(grp_sam,length) # and count the totals # This will now consider empty groups
	grp_all <- grp_less # Since this was modified, we still required the old version
	grp_sam <- grp_sam[(grp_less>=1)] # Remove groups that ended up with 0 samples
	grp_less <- sapply(grp_sam,length) # This will now consider only non-empty groups
	len <- length(unlist(grp_sam)) # Get total items
	grp <- rep("Other",len) # create a template for considering non-group items
	for(i in 1:length(grp_sam)){ # This dual cycle gets the actual sample distribution for all groups
		for(j in 1:length(grp_sam[[i]])){
		grp[grp_sam[[i]][j]] <- names(grp_less)[i];
		}
	}
	out <-list(grp_sam,grp_less,as.factor(grp),grp_all)
	names(out) <- c("grp_sam","grp_less","grp", "grp_all")
	return(out)
}
set_cols <- function(grp_sam, grp_less){
	# Assign colors (if no larger than the input color vectors)
	group_cols <- group_cols[1:length(groups)][grp_less>0] # Trim to grp total, then mask those not present
	len <- length(unlist(grp_sam))
	cols <- rep("gray",len) # set default colors
	for(i in 1:length(grp_sam)){ # This dual cycle gets the actual color for that sample
		for(j in 1:length(grp_sam[[i]])){
		cols[grp_sam[[i]][j]] <- group_cols[i];
		}
	}
	rep <- sapply(grp_sam,length) # ALTERNATIVE VERSION: this uses a block color (for ordered samples)
	cols_ord <- rep(group_cols,times=rep)
	out <-list(group_cols,cols,cols_ord)
	names(out) <- c("group_cols","cols","cols_ord")
	return(out)
}

cluster_consistency <- function(in_matrix,presets) {
	# UPDATE 2020-09-04: Fixed error when group items total 1, then skip those
	mask <- unlist(lapply(presets$grp_sam, length))>1
	all_labels <- names(presets$grp_less[mask]) # This is imported from outside
	beta_dist <- data.frame(matrix(NA,nrow=1,ncol=length(all_labels)*4)) # Create an empty matrix
	group_list <- vector(mode = "list", length = length(all_labels)*2) # as well as an empty list
	names(beta_dist) <- c(paste0(all_labels,"_w"),paste0(all_labels,"_b"),paste0(all_labels,"_pval"),paste0(all_labels,"_qval"))
	for(i in 1:length(all_labels)){
	# 	print(all_labels[i])
		list <- unlist(presets$grp_sam[mask][i]) # get the list of the current group
		sub <- in_matrix[list,list] # subset the matrix to retain only those in the group
		within <- sub[lower.tri(sub)] # keep only the lower triangle (of a symmetric matrix)
		between <- unlist(in_matrix[list,setdiff((1:nrow(in_matrix)),list)], use.names = FALSE) # get the rows of the current group but the complement in the columns (all distances between the group and other samples)
		min_set <- min(length(between),length(within)) # group and non-group comparisons are normally uneven in numbers, keep track of the smallest
		group_list[[seq(1,length(group_list),2)[i]]] <- within # save the corresponding collections
		group_list[[seq(2,length(group_list),2)[i]]] <- between
		beta_dist[1,i] <- median(within) # and also the medians, which may vary a lot
		beta_dist[1,i+length(all_labels)] <- median(between)
		MW_test <- function() {wilcox.test(sample(within, min_set, replace=FALSE),sample(between, min_set, replace=FALSE))$p.value} # create a function to carry out a mann whitney test for pair, with a sample for n repetitions
		MW_all <- replicate(1000, MW_test()) # to save some time, calculate with 1000 items
		beta_dist[1,i+(length(all_labels)*2)] <- median(MW_all) # Save the medians for the pvalue
		beta_dist[1,i+(length(all_labels)*3)] <- median(p.adjust(MW_all, method="fdr")) # and its adjusted pvalue
		beta_dist
	}
	y <- trunc(seq(min(in_matrix),1,(1-min(in_matrix))/20)*1000)/1000
	pdf(paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"between_within_groups.pdf", sep="_"))
	boxplot(group_list, border=c("cornflowerblue","coral1"),col=c("darkslategray2","bisque"),main=paste("Between (b) and within (w) dissimilarities per group.",metric),frame=FALSE,outline=FALSE,yaxt='n',xaxt='n')
	mtext(prefix)
	axis(1,las=2,at=seq(1:(length(all_labels)*2)),labels=as.character(rbind(paste0(all_labels,"_w"), paste0(all_labels,"_b"))))
	axis(2,las=1,at=y)
	dev.off()
	write.table(cbind(prefix,beta_dist),paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"between_within_groups.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

cluster_groups <- function(in_matrix,presets) {
	# UPDATE 2020-09-04: Fixed error when group items total 1, then skip those
	mask <- unlist(lapply(presets$grp_sam, length))>1
	all_labels <- names(presets$grp_less[mask]) # This is imported from outside
	perm <- cbind(rbind(all_labels, all_labels), combn(all_labels,2,simplify=T))
	A <- perm[1,1]
	B <- perm[2,1]
	feats <- combn(unique(c(unlist(presets$grp_sam[A]),unlist(presets$grp_sam[B]))),2)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	beta_dist <- data.frame(matrix(NA,nrow=1,ncol=length(all_labels)*4)) # Create an empty matrix
	group_list <- vector(mode = "list", length = length(all_labels)*2) # as well as an empty list
	names(beta_dist) <- c(paste0(all_labels,"_w"),paste0(all_labels,"_b"),paste0(all_labels,"_pval"),paste0(all_labels,"_qval"))
	for(i in 1:length(all_labels)){
	# 	print(all_labels[i])
		list <- unlist(presets$grp_sam[mask][i]) # get the list of the current group
		sub <- in_matrix[list,list] # subset the matrix to retain only those in the group
		within <- sub[lower.tri(sub)] # keep only the lower triangle (of a symmetric matrix)
		between <- unlist(in_matrix[list,setdiff((1:nrow(in_matrix)),list)], use.names = FALSE) # get the rows of the current group but the complement in the columns (all distances between the group and other samples)
		min_set <- min(length(between),length(within)) # group and non-group comparisons are normally uneven in numbers, keep track of the smallest
		group_list[[seq(1,length(group_list),2)[i]]] <- within # save the corresponding collections
		group_list[[seq(2,length(group_list),2)[i]]] <- between
		beta_dist[1,i] <- median(within) # and also the medians, which may vary a lot
		beta_dist[1,i+length(all_labels)] <- median(between)
		MW_test <- function() {wilcox.test(sample(within, min_set, replace=FALSE),sample(between, min_set, replace=FALSE))$p.value} # create a function to carry out a mann whitney test for pair, with a sample for n repetitions
		MW_all <- replicate(1000, MW_test()) # to save some time, calculate with 1000 items
		beta_dist[1,i+(length(all_labels)*2)] <- median(MW_all) # Save the medians for the pvalue
		beta_dist[1,i+(length(all_labels)*3)] <- median(p.adjust(MW_all, method="fdr")) # and its adjusted pvalue
		beta_dist
	}
	y <- trunc(seq(min(in_matrix),1,(1-min(in_matrix))/20)*1000)/1000
	pdf(paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"between_within_groups.pdf", sep="_"))
	boxplot(group_list, border=c("cornflowerblue","coral1"),col=c("darkslategray2","bisque"),main=paste("Between (b) and within (w) dissimilarities per group.",metric),frame=FALSE,outline=FALSE,yaxt='n',xaxt='n')
	mtext(prefix)
	axis(1,las=2,at=seq(1:(length(all_labels)*2)),labels=as.character(rbind(paste0(all_labels,"_w"), paste0(all_labels,"_b"))))
	axis(2,las=1,at=y)
	dev.off()
	write.table(cbind(prefix,beta_dist),paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"between_within_groups.tsv", sep="_"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
}

####################################### MAIN #######################################
############ Part 1: Assess the group distribution, calculate stats ############
 ### Prefilter samples ###
# library("vegan")
presets <- prepare_data(df) # Define how samples from each group are distributed
# Remove samples not in groups requested by the user:
df <- df[unlist(presets$grp_sam),] # and subset the matrix rowwise
df <- df[,unlist(presets$grp_sam)] # as well as columnwise
presets <- prepare_data(df) # Now recalculate presets
# Create color scheme (this will work if the total colors, defined above, is the correct one)
if(length(groups)<=length(group_cols)){ # Reassign colors to existing groups (to make them the same between graphs). This is only carried out if the number of defined colors (solid and border) is the same than the number of groups.
	cols <- set_cols(presets$grp_sam,presets$grp_all)
	s_col <- cols$cols
	g_col <- cols$group_cols
} else {
	s_col <- g_col <- "gray"
}

############ Part 2: Calculate between and within group correlations ############
# Calculate the within and between median for each group
cluster_consistency(df,presets)
# write.table(cbind(prefix,metric,compare),paste(prefix,metric,"grps",gsub(", ","_",toString(groups)),"stats.tsv", sep="_"), sep="\t", quote=FALSE, row.names=F) 

############ Part 3: Plot correlations between groups ############

print("Done")
