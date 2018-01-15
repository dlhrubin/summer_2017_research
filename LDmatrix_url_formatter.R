#Danielle Rubin
#
#Purpose: R script for selecting sites based on a chosen bin size (e.g. 100kb)
#and formatting them to input into the NIH LDmatrix tool 
#(https://analysistools.nci.nih.gov/LDlink/)
#
#Inputs: csv with 3 columns- chromosome (CHROM), position (POS), and rs number (ID) 
#corresponding to the mutation in that position
#
#Outputs: csv of formatted url strings corresponding to the LD matrices representing
#each bin of sites


#Reads input csv
data <- read.csv("rsid.csv")
#Vector containing the rs numbers from the csv
rsIDs <- data$ID

#Creating the 100kb bins

#Vector containing the position numbers from the csv
positions <- data$POS
#Calculates the range of the position numbers for splitting them up into bins
chrom_length <- max(positions)-min(positions)
#Sets the size of the bin in bp (here 100,000 b, i.e. 100 kb)
bin_size <- 100000
#Initializes vector of the sites that correspond to the breakpoints between 
#bins, starting with the first site (smallest position) in the list
bin_sites <- c(min(positions))
#Initializes vector of the list indices of the above positions, starting with
#the first site
bin_indices <- c(which(positions == min(positions)))
#Variable corresponding to the first site in a given bin, starting with the 
#first site in the list
start_site <- min(positions)

#Loop to fill the vectors of bin breakpoint sites and their indices
#While the position of the bin start site is smaller than the largest position
#in the list
while (start_site < max(positions)) {
  #Variable corresponding to first site of the next bin. Assigned by selecting
  #the largest position value in the list that is no more than 100,000 bases 
  #larger than the start site
  next_site <- max(positions[which(positions < (start_site + bin_size))])
  #Variable corresponding to the index of the next site
  site_index <- which(positions == next_site)
  #If the site immediately following the start site in the list is more than 
  #100,000 bases larger than the start site, the next_site variable will equal
  #the original start_site, so it must be manually changed
  if (next_site == start_site) {
    next_site <- positions[site_index + 1]
    site_index <- site_index + 1
  }
  #Appends bin_sites vector to add the first site of the next bin
  bin_sites <- c(bin_sites, next_site)
  #Appends the bin_indices vector to add the list index of next_site
  bin_indices <- c(bin_indices, site_index)
  #Assigns next_site as the new start site for creating the next bin
  start_site <- next_site
}

#Checking bin size

#Initializing vector to contain the sizes of all the bins 
#(checking to make sure all bins are ~100kb in size)
bin_sizes <- c()
#Iterating through the vector of positions separating each bin
for (i in 1:(length(bin_sites)-1)) {
  #Appending bin_sizes vector with the size of the ith bin (= the difference
  #of the first site of the next bin and the first site of the current bin)
  bin_sizes <- c(bin_sizes, bin_sites[i+1]-bin_sites[i])
}
#Prints the sizes of the bins that are larger than 100kb
bin_sizes[which(bin_sizes > 100000)]

#Splitting bin breakpoint sites into groups of 100 sites because LDmatrix can
#receive a maximum of 100 sites as input

#Vector containing the indices at which the bin indices will be split into groups of 100
matrix_indices <- seq(1, length(bin_indices), by=99)
#Appending the last bin index in the list of bin indices to the matrix indices
#vector to ensure all sites are included
matrix_indices <- c(matrix_indices, length(bin_indices))
#Initializing list of the rs numbers corresponding to the 100 bin breakpoint 
#sites within each 100-site matrix "chunk" 
rs_matrix_chunks <- list()
#Initializing list of the positions corresponding to the 100 bin breakpoint 
#sites within each 100-site matrix "chunk" 
pos_matrix_chunks <- list()
#Initializing vector of the sizes of each matrix chunk
matrix_chunks_sizes <- c()
#Iterates through the list of matrix breakpoint indices
for (i in 1:(length(matrix_indices)-1)) {
  #Appends pos_matrix_chunks list with the 100 positions corresponding to
  #bin sites contained within the ith matrix chunk
  pos_matrix_chunks[[i]] <- positions[bin_indices[matrix_indices[i]:matrix_indices[i+1]]]
  #Appends matrix_chunks_sizes vector with the size of each matrix chunk  
  matrix_chunks_sizes <- c(matrix_chunks_sizes, length(pos_matrix_chunks[[i]]))
  #Appends rs_matrix_chunks list with the 100 rs numbers corresponding to bin
  #sites contained within the ith matrix chunk
  rs_matrix_chunks[[i]] <- rsIDs[bin_indices[matrix_indices[i]:matrix_indices[i+1]]]
}

#Formatting the url strings needed to programmatically access LDmatrix from 
#bash with the selected rs numbers for each 100-site matrix chunk

#Initializing empty vector containing the url strings for each matrix chunk
urls <- c()
#Iterates through list of rs numbers separated by matrix chunk
for (i in 1:length(rs_matrix_chunks)) {
  #Vector of the rs numbers in the ith matrix chunk
  rsIDs_list <- as.character(unlist(rs_matrix_chunks[[i]]))
  #Variable corresponding to the rs number to be added to the url string, 
  #starting with the first rs number in the matrix chunk
  url_string <- rsIDs_list[1]
  #Iterates through the vector of rs numbers in the ith matrix chunk
  for (j in 2:length(rsIDs_list)) {
    #Appends the url letters and the next rs number to the growing url string
    url_string <- paste(url_string, "%0A", rsIDs_list[j], sep="")
  }
  #Appends the formatted url for the ith matrix chunk to the vector of urls
  #for all matrix chunks
  urls <- c(urls,url_string)
}
#Variable corresponding to data frame of the urls (one url per row)
rsID_urls <- data.frame(urls)
#Writes the url dataframe to a csv to be accessed in bash
write.table(rsID_urls, file = "rsID_urls.csv", col.names = FALSE, row.names = FALSE)