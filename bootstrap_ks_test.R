#Danielle Rubin
#
#Purpose: R script to obtain the K-S test p-value for paired single-cell X chromosome
#somatic mutation count data, bootstrap over K-S test p-values for paired single-cell 
#autosome somatic mutation count data as statistical controls, and compare these
#bootstrapped autosome p-values to the X chromosome p-value to determine whether there
#is a statistically significant mutation imbalance between the paired X chromosomes.
#
#Inputs: csv containing the number of somatic mutations on each chromosome in a pair
#for all 22 chromosomes within each single cell sample. Each column corresponds 
#to one chromosome in a pair (e.g. chrX_a, chrX_b, chr1_a, chr1_b, ...),each row 
#corresponds to a single-cell sample ID, and each table cell contains the somatic 
#mutation count for a given chromosome (column) in a given single-cell sample (row).
#
#Outputs: final p-value calculated by computing the proportion of bootstrapped 
#         autosomal control p-values smaller (i.e. more significant) than the 
#         X chromosome p-value;
#         density plot marking this chrX p-value (red line + labelled) on the
#         distribution of the bootstrapped autosomal control p-values


#Read input csv
paired_data <- read.csv("Paired Mutations per Chromosome.csv")
#Create vector containing all chromosome numbers
chrom_nums <- c("X")
for (i in 1:22) {
  chrom_nums[i+1] <- i
}
#Create vector containing all chromosome names
chroms <- sapply(chrom_nums, function(x) paste("chr", x, sep=""))
p_values <- c()
#Returns vectors containing the mean differences and K-S test p-values for
#each chromosome pair
for (i in 1:length(chroms)) {
  index_a <- which(colnames(paired_data) == paste(chroms[i], "_1", sep=""))
  index_b <- which(colnames(paired_data) == paste(chroms[i], "_2", sep=""))
  chr_a_counts <- paired_data[[index_a]]
  if (length(which(chr_a_counts == "not used")) > 0) {
    chr_a_counts <- chr_a_counts[-which(chr_a_counts == "not used")]
  }
  chr_a_counts <- as.numeric(as.character(unlist(chr_a_counts)))
  chr_b_counts <- paired_data[[index_b]]
  if (length(which(chr_b_counts == "not used")) > 0) {
    chr_b_counts <- chr_b_counts[-which(chr_b_counts == "not used")]
  }
  chr_b_counts <- as.numeric(as.character(unlist(chr_b_counts)))
  p_values <- c(p_values, ks.test(chr_a_counts, chr_b_counts)$p.value)
}
#Chromosome X p-value
X_p_value <- p_values[1]

#Generating bootstrapped autosomal p-values

#Parses autosomal paired mutation count data into a list of lists where
#each entry is a list corresponding to a single-cell sample, each entry
#of which is a vector with the two mutation counts for each autosome pair
all_cells_counts <- list()
for (i in 1:nrow(paired_data)) {
  bootstrap_vec <- paired_data[i,4:ncol(paired_data)]
  paired_counts <- list()
  for (j in 2:length(chroms)) {
    index_a <- which(names(bootstrap_vec) == paste(chroms[j], "_1", sep=""))
    index_b <- which(names(bootstrap_vec) == paste(chroms[j], "_2", sep=""))
    paired_counts[[j]] <- c((as.character(bootstrap_vec[1,index_a])), (
                            as.character(bootstrap_vec[1,index_b])))
  }
  paired_counts <- paired_counts[2:length(paired_counts)]
  if (length(paired_counts[-grep("not used", paired_counts)]) > 0) {
    paired_counts <- paired_counts[-grep("not used", paired_counts)]
  }
  all_cells_counts[[i]] <- paired_counts
}
#Function that bootstraps by pooling together all autosomal p-values across all
#samples, preserving pairings, into two conglomerated lists and randomly "flips"
#paired values between the two lists to simulate autosomes with randomly distributed
#mutation counts as controls, returning a K-S test p-value for the flipped lists
bootstrapper <- function() {
  bootstrapped_list <- c()
  for (i in 1:length(all_cells_counts)) {
    cell_counts <- all_cells_counts[[i]]
    bootstrapped_list[[i]] <- as.numeric(unlist(sample(cell_counts, 1, 
                                                       replace=FALSE)))
  }
  chr_a_counts <- as.numeric(sapply(bootstrapped_list, "[[", 1))
  chr_b_counts <- as.numeric(sapply(bootstrapped_list, "[[", 2))
  p_value <- ks.test(chr_a_counts, chr_b_counts)$p.value
  return(p_value)
}

#Performs 10000 calls to bootstrapper(), generating 10000 control p-values
p_value_array <- replicate(10000, bootstrapper())
#Calculates final p-value used for determing significance by measuring the proportion
#of bootstrapped autosomal p-values smaller than the X chromosome p-value
bootstrapped_p_value <- sum(p_value_array <= X_p_value)/length(p_value_array)
print(bootstrapped_p_value)

#Generates density plot of autosomal bootstrapped p-values with the X chromosome 
#p-value marked with a labelled red line on the distribution
p_value_label <- textGrob("chrX p-value", gp=gpar(fontsize=8))
enrichment_graph <- ggplot(data.frame("proportions"=(-log(p_value_array))), 
                    aes(x=proportions)) + geom_density(color="lightblue",
                    fill="lightblue", aes(x=proportions, y=..density..)) + 
                    theme_classic() + labs(title = 
                    "Autosomal Bootstrapped K-S Test P-Values", x = "-log(P-Value)",
                    y = "Density") + theme(plot.title = element_text(hjust = 0.5)) +
                    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=
                    c(0,0)) + theme(text=element_text(size=13)) + geom_vline(
                    xintercept=(-log(X_p_value)), color="red") + annotation_custom(
                    p_value_label,xmin=-log(X_p_value),xmax=-log(X_p_value),
                    ymin=-0.007,ymax=-0.007)
gt <- ggplot_gtable(ggplot_build(enrichment_graph))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)