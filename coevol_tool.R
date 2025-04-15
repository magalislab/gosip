#!/usr/bin/env Rscript

library(optparse)

# Menu to select tree + metadata type
cat("Choose tree type and metadata source:\n")
main_choice <- menu(c(
  "1: BEAST tree with states and rescaled in time (1 nexus file)",
  "2: BEAST tree (1 nexus file) + Nextstrain metadata (1 .json file)",
  "3: Treetime tree (rescaled in time) + mugration metadata (2 nexus files)"
), title = "Select input type")

option_list <- switch(main_choice,
                      # Option 1
                      list(
                        make_option(c("-ta", "--treeancestral"), type = "character", default = NULL,
                                    help = "Load phylogenetic tree rescaled in time from BEAST (nexus format)"),
                        make_option(c("-d", "--dpi"), type = "integer", default = NULL,
                                    help = "Number of days post infection at necropsy", metavar = "number"),
                        make_option(c("-w", "--window"), type = "integer", default = 30,
                                    help = "Length of time windows in days [default 30 days]"),
                        make_option(c("-o", "--output"), type = "character", default = "Coevol_output",
                                    help = "Name for the output file")
                      ),
                      
                      # Option 2
                      list(
                        make_option(c("-t", "--tree"), type = "character", default = NULL,
                                    help = "Load phylogenetic tree rescaled in time (nexus format)"),
                        make_option(c("-a", "--ancestral"), type = "character", default = NULL,
                                    help = "Load ancestral state reconstruction (.json file)"),
                        make_option(c("-d", "--dpi"), type = "integer", default = NULL,
                                    help = "Number of days post infection at necropsy", metavar = "number"),
                        make_option(c("-w", "--window"), type = "integer", default = 30,
                                    help = "Length of time windows in days [default 30 days]"),
                        make_option(c("-o", "--output"), type = "character", default = "Coevol_output",
                                    help = "Name for the output file")
                      ),
                      
                      # Option 3
                      list(
                        make_option(c("-t", "--tree"), type = "character", default = NULL,
                                    help = "Load phylogenetic tree rescaled in time (nexus format)"),
                        make_option(c("-a", "--ancestral"), type = "character", default = NULL,
                                    help = "Load ancestral state reconstruction (nexus format)"),
                        make_option(c("-d", "--dpi"), type = "integer", default = NULL,
                                    help = "Number of days post infection at necropsy", metavar = "number"),
                        make_option(c("-w", "--window"), type = "integer", default = 30,
                                    help = "Length of time windows in days [default 30 days]"),
                        make_option(c("-o", "--output"), type = "character", default = "Coevol_output",
                                    help = "Name for the output file")
                      ),
                      
                      stop("No valid option selected. Exiting.")
)

# Parse arguments using selected option list
parser <- parse_args(OptionParser(option_list = option_list))


libs <- c("tidyverse", "treeio", "tidytree", "jsonlite",
  "parallel", "ape", "data.table", "Rgraphviz")

lapply(libs, library, character.only = TRUE)

numCores=detectCores()
`%notin%` = Negate(`%in%`)
options(warn=-1)

# === ANALYSIS by input type ===

if (main_choice == 1) {
  message("Option 1 selected: BEAST tree with states and rescaled in time")
  # Load tree
  b <- read.beast(parser$treeancestral)
  b.tree <- as_tibble(b)
  tree_data <- b.tree %>%
    select(node, parent, height, label, states) %>%
    mutate(height = as.numeric(height),
           end_dpi = parser$dpi - height) %>%
    rename(end_tissue = states) %>%
    mutate(origin_dpi = end_dpi[parent],
           origin_tissue = end_tissue[parent])
  
} else if (main_choice == 2) {
  message("Option 2 selected: BEAST tree + Nextstrain metadata")
  
  b <- read.beast(parser$tree)
  b.tree <- as_tibble(b)
  tree_data <- b.tree %>%
    select(node, parent, height, label, states) %>%
    mutate(height = as.numeric(height),
           end_dpi = parser$dpi - height) %>%
    rename(end_tissue = states)
  
  # Load mutation JSON
  mut_info <- read_json(parser$ancestral)$nodes
  mut_info <- lapply(mut_info, `[[`, "muts")
  mut_info <- lapply(mut_info, function(x) {
    x <- do.call(rbind, x)
    data.frame(muts = x)
  })
  mut_info <- rbindlist(mut_info, idcol = "label") %>%
    mutate(site = as.numeric(gsub("[A-Z-](\\d+)[A-Z-]", "\\1", muts))) %>%
    group_by(site) %>%
    mutate(num_branches = n_distinct(label)) %>%
    filter(num_branches >= 2) %>%
    arrange(label, site)
  
  tree_data <- tree_data %>%
    mutate(origin_dpi = end_dpi[parent],
           origin_tissue = end_tissue[parent])
  
  alldata <- merge(tree_data, mut_info)
  
} else if (main_choice == 3) {
  message("Option 3 selected: Treetime tree + mugration metadata")
  
  mutsbranches <- as_tibble(read.beast(parser$tree))
  mut_info <- data.frame()
  for (i in 1:nrow(mutsbranches)) {
    for (j in seq_along(mutsbranches$mutations[[i]])) {
      mut <- mutsbranches$mutations[[i]][j]
      addcol <- cbind(mutsbranches[i, 1:5], mut)
      mut_info <- rbind(mut_info, addcol)
    }
  }
  
  mut_info <- mut_info %>%
    mutate(site = as.numeric(gsub("[A-Z-](\\d+)[A-Z-]", "\\1", mut))) %>%
    group_by(site) %>%
    mutate(num_branches = n_distinct(label)) %>%
    filter(num_branches >= 2) %>%
    arrange(label, site)
  
  mut_states <- mut_info %>%
    rename(end_dpi = date) %>%
    mutate(end_dpi = as.numeric(end_dpi),
           end_dpi = (end_dpi - 2017) * 365)
  
  # Get origin dpi from parent
  mut_states$origin_dpi <- sapply(mut_states$parent, function(p) {
    a <- which(mut_states$node == p)
    if (length(a) == 1) mut_states$end_dpi[a] else NA
  })
  
  mut_states <- mut_states[!is.na(mut_states$site), ]
}

#sliding window approach
origin_dpi <- 0
size <- parser$window
slide <- parser$window
nec <- parser$dpi

Window <- as.data.frame(cbind(Start=seq(origin_dpi,nec-parser$window,slide), End=seq(size,nec,slide)))

mylist.names <- as.character(seq(size,nec,slide))
slid.win < -vector("list", length(seq(origin_dpi,nec-parser$window,slide)))
names(slid.win) < -mylist.names

slid.win <- mclapply(1:nrow(Window), function(win) {
  final <- mclapply(seq_along(alldata$origin_dpi), function(branch) {
    if(alldata$origin_dpi[branch]<=Window$End[win] & 
       alldata$end_dpi[branch]>=Window$Start[win] ) {
      result = alldata[branch,]
    } else {
      result=NULL
    }
    return(result)
  }, mc.cores=numCores)
  x=do.call(rbind, compact(final))
  return(data.frame("X"=x))
}, mc.cores=numCores)

#changing names of dataframes to dpi
names(slid.win) <- mylist.names
#changing colnames back to original in dataframes in the list of windows
slid.win <- lapply(slid.win, setNames, colnames(alldata))

#colnames(slid.win$`120`) = colnames(alldata)  

require(data.table)

#this will put the sequence only in the last window they appear
slid.win.max <- rbindlist(slid.win, idcol="window") %>%
  dplyr::mutate(window=as.numeric(window)) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(max_win=max(window)) %>%
  dplyr::filter(window==max_win)

#this will put the sequence only in the first window they appear
slid.win.min <- rbindlist(slid.win, idcol="window") %>%
  dplyr::mutate(window=as.numeric(window)) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(min_win=min(window)) %>%
  dplyr::filter(window==min_win)

# Gaps often are not isolated; they are found in stretches, so want to consolidate our gapped sites into discrete stretches by
# first identifying mutations that result in gaps (still split up by branch):
deletions.max <- filter(slid.win.max, grepl("-$", muts)) %>%
  group_by(label) %>%
  group_split()

deletions.min <- filter(slid.win.min, grepl("-$", muts)) %>%
  group_by(label) %>%
  group_split()

for (i in 1) {
  if (is_empty(deletions.max) == TRUE){
    df.max <- slid.win.max
    next
  } else {
    deletions.max <- lapply(deletions.max, function(x) {
      site_batch <- lapply(2:nrow(x), function(y) {
        if(isTRUE(x$site[y] == x$site[y-1]+1)) {
          return(NA)
          } else {
            return(x$site[y])
            }
        })
      site_batch <- c(x$site[1], do.call(rbind, site_batch))
      result <- cbind(x, site_batch)
      return(result)
      })
# Now categorize gap-related mutations as insertions or deletions
    deletions.max <- lapply(deletions.max, function(x) {
      fill(x, site_batch) %>%
        group_by(label, site_batch) %>%
        mutate(batch_length = length(site),
                  batch_end = site_batch + batch_length-1) %>%
      distinct()
    })
# Now take stretches of gaps and concatenate them into one discrete deleted region 
# This is done so that when we look for co-evolving pairs (normally done for individual nucleotide sites),
# We can treat a deletion region as an individual "site" 
    deletions.max <- do.call(rbind, lapply(deletions.max, function(x) {
      x <- x %>%
        ungroup() %>%
        mutate(binned = if_else(batch_length>1, 
                              as.character(paste0(site_batch, "-", batch_end)), 
                              as.character(site_batch))) #%>%
    })) %>%
      group_by(binned) %>%
      mutate(num_branches = length(label))
#Now remove indels from the larger mutation datatable and merge it with the deletions table
    require(dplyr)
    df.max <- merge(slid.win.max, deletions.max, by = c('label','muts','site','window','max_win','node','parent',
                            'origin_dpi','end_dpi'), all = T)#,'origin_tissue','end_tissue'
write.csv(df.max, paste0(parser$output,"branch-site_mutations_max.win.csv"), quote=F, row.names=F)
  }
}



for (j in 1) {
  if (is_empty(deletions.min) == TRUE){
    df.min <- slid.win.min
    next
  } else {
    deletions.min <- lapply(deletions.min, function(x) {
      site_batch <- lapply(2:nrow(x), function(y) {
        if(isTRUE(x$site[y] == x$site[y-1]+1)) {
          return(NA)
        } else {
          return(x$site[y])
        }
      })
      site_batch <- c(x$site[1], do.call(rbind, site_batch))
      result <- cbind(x, site_batch)
      return(result)
    })
    # Now categorize gap-related mutations as insertions or deletions
    deletions.min <- lapply(deletions.min, function(x) {
      fill(x, site_batch) %>%
        group_by(label, site_batch) %>%
        mutate(batch_length = length(site),
               batch_end = site_batch + batch_length-1) %>%
        distinct()
    })
    # Now take stretches of gaps and concatenate them into one discrete deleted region 
    # This is done so that when we look for co-evolving pairs (normally done for individual nucleotide sites),
    # We can treat a deletion region as an individual "site" 
    deletions.min <- do.call(rbind, lapply(deletions.min, function(x) {
      x <- x %>%
        ungroup() %>%
        mutate(binned = if_else(batch_length>1, 
                                as.character(paste0(site_batch, "-", batch_end)), 
                                as.character(site_batch)))
    })) %>%
      group_by(binned) %>%
      mutate(num_branches = length(label))
    #Now remove indels from the larger mutation datatable and merge it with the deletions table
    require(dplyr)
    df.min <- merge(slid.win.min, deletions.min, by = c('label','muts','site','window','min_win','node','parent',
                                                       'origin_dpi','end_dpi'), all = T)#,'origin_tissue','end_tissue'
    
    write.csv(df.min, paste0(parser$output,"branch-site_mutations_min.win.csv"), quote=F, row.names=F)
  }
}


# Now, for each branch, lump changing sites into site pairs by using the following function to create all possible pairwise combinations
expand.grid.unique <- function(x, y, include.equals=FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

##applying BGM
require(dplyr)
require(bnlearn)

max.per_window <- map(mylist.names, ~ subset(df.max, window == .x))
min.per_window <- map(mylist.names, ~ subset(df.min, window == .x))

max.namescol <- colnames(max.per_window[[1]])
max.per_window <- Reduce(function(...) merge(..., by=max.namescol, all=TRUE), max.per_window)

min.namescol <- colnames(min.per_window[[1]])
min.per_window <- Reduce(function(...) merge(..., by=min.namescol, all=TRUE), min.per_window)


#counting how many times a site mutation happen in each time window
max.ctable <- table(dplyr::select(df.max, window, site))
#making into dataframe for visual analysis
df.max.ctable <- as.data.frame.matrix(max.ctable)

#Calculating the value of each cell in a table as a proportion of all values
#to use as input for bnlearn
max.ptable <- prop.table(max.ctable)
df.max.ptable <- as.data.frame.matrix(max.ptable)

#Extracting sites from contingency table that have mutations in more than 1 time window
#setting a list of the columns with their count of zeros
max.zerocount <- as.data.frame(t(colSums(df.max.ctable == 0)))
max.rownum <- ceiling(nrow(df.max.ctable)/2)
#Now if a column has more than half the rows with 0, then will not bt taken
#into account for whitelist
compare_fun <- if (nrow(df.max.ctable) %% 2 == 0) `<=` else `<`
max.forwhitelist <- names(max.zerocount)[compare_fun(max.zerocount, max.rownum)]

#creating all possible combinations of pairs of sites
max.whitelist <- expand.grid(from = max.forwhitelist,
  to = max.forwhitelist,
  stringsAsFactors = FALSE) %>%
  filter(from != to)

#applying bayesian graphical model with whitelist of sites
suppressWarnings({ max.res = iamb(df.max.ptable, whitelist = max.whitelist) })
max.result_pairs = max.res$arcs

library(ggplot2)
library(igraph)
library(RColorBrewer)


for (pair in 1) {
  if(is_empty(max.result_pairs)==TRUE){
    print("No pairs of sites initially found in maximum window")
    next
  }else if (nrow(max.ctable)<=1){
    print("All sites have been assigned to the same window. No relevant pairs of sites can be determined")
    next
  } else {
    max.neg_corr=do.call(rbind,mclapply(1:nrow(max.ctable), function(win) {
      do.call(rbind,lapply(1:nrow(max.result_pairs),function(site_pair) {
        val1=max.ctable[win, max.result_pairs[site_pair,1]]
        val2=max.ctable[win, max.result_pairs[site_pair,2]]
        if(val1==0 & val2==0) {
          neg_corr = data.frame(from=max.result_pairs[site_pair,1],
                                to=max.result_pairs[site_pair,2], score=0, row.names = NULL)
        } else if (val1 >0 & val2 >0) {
          neg_corr = data.frame(from=max.result_pairs[site_pair,1],
                                to=max.result_pairs[site_pair,2], score=1, row.names=NULL)
        } else {
          neg_corr = data.frame(from=max.result_pairs[site_pair,1],
                                to=max.result_pairs[site_pair,2], score=0, row.names=NULL)
        }
      }))
    }, mc.cores=numCores))
    
    
    #Counting the pairs with a score of 0 and the score of 1
    require(plyr)
    max.zeros = ddply(max.neg_corr,.(from,to,score),nrow)
    detach("package:plyr", unload = TRUE)
    #We can have co-evolving sites where sometimes one site happens without the other
    #Grouping into blacklist according to pairs and find pairs with majority having score of zero
    #This way, it will blacklist the pairs in which the times the score of 0 is more than the same pair that it was scored as 1
    #Counting which pairs have a majority score of 0
    if (nrow(df.max.ctable) %% 2 == 0){
      torem.max = which(max.zeros$score == 0 & max.zeros$V1 > max.rownum)
    } else {
      torem.max = which(max.zeros$score == 0 & max.zeros$V1 >= max.rownum)
    }
    
    removing.max = data.frame()  
    for(i in 1:length(torem.max)){
      removing.max = rbind(removing.max, max.zeros[torem.max[i],])
    }
    #df with pairs only from original df
    max.zeros1 = unique(max.zeros[,1:2])
    #df that contains pairs only now to remove
    max.removing = unique(removing.max[,1:2])
    
    #Pairs with majority of score 0
    #blacklist now will be what is in the df that needs to be removed that is also
    #in original df counting 0 and 1
    max.blacklist = inner_join(max.removing, max.zeros1)
    
    #if a pair of sites appears in both blacklist and whitelist, bnlearn
    #whitelists it. anti_join is looking for what is in whitelist that is
    #not in blacklist to use as new whitelist
    max.whitelist1 = anti_join(max.whitelist, max.blacklist, by=c('from','to'))
    
    # Remove these blacklisted pairs from the original probability table as well
    # After blacklisting, filter pairs based on bootstrap strength:
    suppressWarnings({ max.res=iamb(df.max.ptable, whitelist = max.whitelist1,
                                    blacklist=max.blacklist) })
    max.result_pairs=as.data.frame(max.res$arcs)
    #otaining the arc direction and strength
    suppressWarnings({ bs.max=boot.strength(df.max.ptable, algorithm="iamb") })
    
    #assigning arc direction to pairs of sites in max window
    final.max=do.call(rbind, compact(lapply(split(bs.max, seq(nrow(bs.max))), function(x) {
      x2=paste(x$from,x$to, sep=",")
      bl=do.call(paste, c(max.blacklist, sep=","))
      if(x2 %in% bl) {
        return(NULL)
      } else {
        return(x)
      }
    })) )
    
    #getting those sites that are only in max.whitelist1
    final.max1 = inner_join(final.max, max.whitelist1, by=c('from','to'))
    #filtering to only those sites that significant with a strength of >0
    final.max1 = final.max1 %>% filter(strength > 0)
    fromto = final.max1[,1:2]
    
    if (nrow(fromto)==0){
      fromto = data.frame(from = 0, to = 0)
    }
    
    write.csv(fromto, file = paste0(treename,"maxwin_from_to_sites.csv"), eol = '\n', row.names = T)
    
    #cross referencing to keep sequences that have significant co-evolution 
    #but also not on the same branch
    max.keep.seq = do.call(rbind, mclapply(1:nrow(slid.win.max), function(newid) {
      do.call(rbind, lapply(1:nrow(final.max1), function(sigsite){
        s1=final.max1[sigsite,1]
        s2=final.max1[sigsite,2]
        if(slid.win.max$site[newid] %in% final.max1[sigsite,]){
          keep.seq=data.frame(slid.win.max[newid,])
        }
      }))
    }, mc.cores=numCores))
    
    max.keep.seq = unique(max.keep.seq)
    
    write.csv(max.keep.seq, file = paste0(treename,"maxwin_info.csv"), eol = '\n', row.names = T)
    
    max.sitesonly = as.data.frame(unique(max.keep.seq$site))
    #print(max.sitesonly)
    if (nrow(max.sitesonly)==0){
      max.sitesonly = append(max.sitesonly, 0)
    } else {
      colnames(max.sitesonly)=treename
    }
    
    write.table(max.sitesonly,file = paste0(treename,"_maxwin.sites_only.csv"), eol = '\n', row.names = F, append = T)
    
    # Plots distribution of strengths
    strength_distribution.max = ggplot(final.max1, aes(x=strength)) +
      geom_density() + xlab("Strength (posterior)") + ylab("% total pairs of sites") +
      theme_classic()
    ggsave(paste0(parser$output,"_max.win_strength_distribution.png"), strength_distribution.max, width = 8, height = 8, units = "in", dpi = 600)
    
    #plot network
    network.max = strength.plot(max.res, strength = bs.max, main = "Max_window",
                                threshold = 0.001,
                                shape = "circle",
                                layout = "fdp")
    
    require(igraph)
    network.max1 = graph_from_graphnel(network.max, name = T) 
    network.max1 = delete.edges(network.max1, which(E(network.max1)$weight <0.09))
    png(paste0(parser$output,"_max.win_sites.png"), 1024,1024)
    plot(delete.vertices(simplify(network.max1), degree(network.max1)==0), edge.color="palevioletred4", vertex.label.color="black", vertex.label.font=1,
         vertex.color = "lightpink2", edge.arrow.size = 1, vertex.size = 25,
         vertex.label.cex=2, edge.lty = 1)
    dev.off()
    
  }
}

min.ctable = table(dplyr::select(df.max, window, site))
#making into dataframe for visual analysis
df.min.ctable = as.data.frame.matrix(min.ctable)

#Calculating the value of each cell in a table as a proportion of all values
#to use as input for bnlearn
min.ptable = prop.table(min.ctable)
df.min.ptable = as.data.frame.matrix(min.ptable)

#Extracting sites from contingency table that have mutations in more than 1 time window
#setting a list of the columns with their count of zeros
min.zerocount <- as.data.frame(t(colSums(df.min.ctable == 0)))
min.rownum <- ceiling(nrow(df.min.ctable)/2)
#Now if a column has more than half the rows with 0, then will not bt taken
#into account for whitelist
compare_fun <- if (nrow(df.min.ctable) %% 2 == 0) `<=` else `<`
min.forwhitelist <- names(min.zerocount)[compare_fun(min.zerocount, min.rownum)]

#creating all possible combinations of pairs of sites
min.whitelist <- expand.grid(from = min.forwhitelist,
                             to = min.forwhitelist,
                             stringsAsFactors = FALSE) %>%
  filter(from != to)

#applying bayesian graphical model with whitelist of sites
suppressWarnings({ min.res = iamb(df.min.ptable, whitelist = min.whitelist) })
min.result_pairs = min.res$arcs


for (pair in 1) {
  if(is_empty(min.result_pairs)==TRUE){
    print("No pairs of sites initially found in maximum window")
    next
  }else if (nrow(min.ctable)<=1){
    print("All sites have been assigned to the same window. No relevant pairs of sites can be determined")
    next
  } else {
    min.neg_corr=do.call(rbind,mclapply(1:nrow(min.ctable), function(win) {
      do.call(rbind,lapply(1:nrow(min.result_pairs),function(site_pair) {
        val1=min.ctable[win, min.result_pairs[site_pair,1]]
        val2=min.ctable[win, min.result_pairs[site_pair,2]]
        if(val1==0 & val2==0) {
          neg_corr = data.frame(from=min.result_pairs[site_pair,1],
                                to=min.result_pairs[site_pair,2], score=0, row.names = NULL)
        } else if (val1 >0 & val2 >0) {
          neg_corr = data.frame(from=min.result_pairs[site_pair,1],
                                to=min.result_pairs[site_pair,2], score=1, row.names=NULL)
        } else {
          neg_corr = data.frame(from=min.result_pairs[site_pair,1],
                                to=min.result_pairs[site_pair,2], score=0, row.names=NULL)
        }
      }))
    }, mc.cores=numCores))
    
    require(dplyr)
    #Counting the pairs with a score of 0 and the score of 1
    min.zeros = ddply(min.neg_corr,.(from,to,score),nrow)
    
    #Counting which pairs have a majority score of 0
    if (nrow(df.min.ctable) %% 2 == 0){
      torem.min = which(min.zeros$score == 0 & min.zeros$V1 > min.rownum)
    } else {
      torem.min = which(min.zeros$score == 0 & min.zeros$V1 >= min.rownum)
    }
    
    removing.min = data.frame()  
    for(i in 1:length(torem.min)){
      removing.min = rbind(removing.min, min.zeros[torem.min[i],])
    }
    #df with pairs only from original df
    min.zeros1 = unique(min.zeros[,1:2])
    #df that contains pairs only now to remove
    min.removing = unique(removing.min[,1:2])
    #Pairs with majority of score 0
    #blacklist now will be what is in the df that needs to be removed that is also
    #in original df counting 0 and 1
    min.blacklist = inner_join(min.removing, min.zeros1)
    
    #if a pair of sites appears in both blacklist and whitelist, bnlearn
    #whitelists it. anti_join is looking for what is in whitelist that is
    #not in blacklist to use as new whitelist
    min.whitelist1 = anti_join(min.whitelist, min.blacklist, by=c('from','to'))
    
    # Remove these blacklisted pairs from the original probability table as well
    # After blacklisting, filter pairs based on bootstrap strength:
    suppressWarnings({ min.res=iamb(df.min.ptable, whitelist = min.whitelist1,blacklist=min.blacklist) })
    min.result_pairs=as.data.frame(min.res$arcs)
    #otaining the arc direction and strength
    suppressWarnings({ bs.min=boot.strength(df.min.ptable, algorithm="iamb") })
    
    #assigning arc direction to pairs of sites in max window
    final.min=do.call(rbind, compact(lapply(split(bs.min, seq(nrow(bs.min))), function(x) {
      x2=paste(x$from,x$to, sep=",")
      bl=do.call(paste, c(min.blacklist, sep=","))
      if(x2 %in% bl) {
        return(NULL)
      } else {
        return(x)
      }
    })) )
    
    #getting those sites that are only in min.whitelist1
    final.min1 = inner_join(final.min, min.whitelist1, by=c('from','to'))
    #filtering to only those sites that significant with a strength of >0
    final.min1 = final.min1 %>% filter(strength > 0)
    fromto = final.min1[,1:2]
    
    if (nrow(fromto)==0){
      fromto = data.frame(from = 0, to = 0)
    }
    
    write.csv(fromto, file = paste0(treename,"maxwin_from_to_sites.csv"), eol = '\n', row.names = T)
    
    #cross referencing to keep sequences that have significant co-evolution 
    #but also not on the same branch
    min.keep.seq = do.call(rbind, mclapply(1:nrow(slid.win.min), function(newid) {
      do.call(rbind, lapply(1:nrow(final.min1), function(sigsite){
        s1=final.min1[sigsite,1]
        s2=final.min1[sigsite,2]
        if(slid.win.min$site[newid] %in% final.min1[sigsite,]){
          keep.seq=data.frame(slid.win.min[newid,])
        }
      }))
    }, mc.cores=numCores))
    
    min.keep.seq = unique(min.keep.seq)
    
    write.csv(min.keep.seq, file = paste0(treename,"maxwin_info.csv"), eol = '\n', row.names = T)
    
    min.sitesonly = as.data.frame(unique(min.keep.seq$site))
    #print(min.sitesonly)
    if (nrow(min.sitesonly)==0){
      min.sitesonly = append(min.sitesonly, 0)
    } else {
      colnames(min.sitesonly)=treename
    }
    
    write.table(min.sitesonly,file = paste0(treename,"_maxwin.sites_only.csv"), eol = '\n', row.names = F, append = T)
    
    # Plots distribution of strengths
    strength_distribution.min = ggplot(final.min1, aes(x=strength)) +
      geom_density() + xlab("Strength (posterior)") + ylab("% total pairs of sites") +
      theme_classic()
    ggsave(paste0(parser$output,"_min.win_strength_distribution.png"), strength_distribution.min, width = 8, height = 8, units = "in", dpi = 600)
    
    #plot network
    network.min = strength.plot(min.res, strength = bs.min, main = "Max_window",
                                threshold = 0.001,
                                shape = "circle",
                                layout = "fdp")
    
    require(igraph)
    network.min1 = graph_from_graphnel(network.min, name = T) 
    network.min1 = delete.edges(network.min1, which(E(network.min1)$weight <0.09))
    png(paste0(parser$output,"_min.win_sites.png"), 1024,1024)
    plot(delete.vertices(simplify(network.min1), degree(network.min1)==0), edge.color="darkblue", vertex.label.color="black", vertex.label.font=1,
         vertex.color = "lightblue", edge.arrow.size = 1, vertex.size = 25,
         vertex.label.cex=2, edge.lty = 1)
    dev.off()
    
  }
}

