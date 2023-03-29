# Started on 2020-06-03 by Rodrigo García-López
# df <- read.table("test_repeated_samples/raw_Z1_Z2-HE1-rar-1000-depth-12469-perm.tsv",comment.char = "", sep = "\t", header = TRUE, row.names = 1,skip = 0,quote="",fill = FALSE)

library("corrr")
# library("ggplot2")
args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) { # five arguments are required
  stop("A minimum of 2 arguments are required", call.=FALSE)
}
out_name <- as.character(args[1])
sample <- as.character(args[2])
df <- read.table(file('stdin'),comment.char = "", sep ="\t", header = TRUE, row.names = 1,skip = 0,quote="",fill = FALSE)
subjects <- correlate(df, method = "spearman", use = "pairwise.complete.obs",diagonal = 1) #load a subject correlation object
out <- focus(subjects,paste0("X",sample))
write.table(out,out_name, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
