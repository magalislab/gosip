#!/bin/bash

require(rlist)
require(tidyverse)
require(phytools)
require(parallel)

numCores=detectCores()
set.seed(123)

`%notin%` <- Negate(`%in%`)

## Parameters #####################################################
n=c("A","C","G","T")

#alignment from which random sequences are taken
# Make sure in proper reading frame
VS=as.list(as.data.frame(t(toupper(read.dna("VS.fas", format="fasta", as.character=T)))))
#picking random seq from alignment to evolve
starting_seq=sample(VS, 1)[[1]]
#taking only the first position site for all codons
first_codon=seq(1,length(starting_seq),3)

s=sample(2:10,1) # Number of sites in network
network=sample(first_codon,s)

#matreix for coevolving sites
suppressWarnings({ ## We know it fills diagonals with NAs
  network_rates <- diag(unique(network), nrow = length(unique(network))) 
})
diag(network_rates) = 0
non.diag <- 1/(nrow(network_rates)-1)
network_rates[lower.tri(network_rates)] <- non.diag
network_rates[upper.tri(network_rates)] <- non.diag
colnames(network_rates) <- rownames(network_rates) <- unique(network)

#number of gens to evolve
ngen=180

r=3E-05 # Evolutionary rate of HIV in subs/site/gen

###################################################################
sim_fct = function(starting_seq, ngen, n, r, network, network_rates) {
  
  sites_not_in_net = 1:length(starting_seq)
  sites_not_in_net = sites_not_in_net[sites_not_in_net %notin% network]
  
  r_0=1-(r*3)
  
  seqs = list(list(seq=starting_seq, gen=1))
  names(seqs)=c("1")
  
  i=2
  while (i<=ngen & length(seqs)<1E05) {
    seqs_prev = seqs[as.numeric(names(seqs))==i-1]
    l_prev=length(seqs_prev)

    if(i>7) {
      n_sample = round(0.50*l_prev)
    } else {
      n_sample = l_prev
    }
    seqs_prev = sample(seqs_prev, n_sample, replace=T) # Replacement allows for some sequences to be more fit than others
    # Although this fitness is random, which is not very realistic
    new_seq1 = seqs_prev
    new_seq2 = seqs_prev
    
    for (j in seq_along(seqs_prev)){
      for (k in sites_not_in_net) {
        new_seq1[[j]]$seq[k] = sample(
          c(n[n==new_seq1[[j]]$seq[k]],n[n!=new_seq1[[j]]$seq[k]]), 1, replace=F, c(r_0,rep(r, 3)))
          
          n_s=sample(network, 1)

          n1=new_seq1[[j]]$seq[n_s]
          new_seq1[[j]]$seq[n_s] = sample(
            c(n[n==new_seq1[[j]]$seq[n_s]],n[n!=new_seq1[[j]]$seq[n_s]]), 1, replace=F, c(r_0,rep(r, 3)))
          n2=new_seq1[[j]]$seq[n_s]
          new_seq1[[j]]$gen=i
          
          if(n1 != n2) {
          changed=c(j,n_s)
        } else {
          changed=NA
        } # End if-else statement
        new_seq2[[j]]$seq[k] = sample(
          c(n[n==new_seq2[[j]]$seq[k]],n[n!=new_seq2[[j]]$seq[k]]), 1, replace=F, c(r_0,rep(r, 3)))
        new_seq2[[j]]$gen=i
      } # End loop along non-network sites (k)
      if(length(changed)>1) {
        
        nr_row = network_rates[as.character(n_s),]
        sub_net = sample(network, 1, replace=F, prob=nr_row)
        
        p=sample(0:1, 1, replace=T, prob=c(0.05, 0.95))  # 95% probability that co-variation event occurs   
        if(p==1) {
          j2 = sample(1:length(new_seq2), 1) #Pick a sequence at random and co-evolve
          new_seq2[[j2]]$seq[sub_net] = sample(n[n!=new_seq2[[j2]]$seq[sub_net]], 1)
          new_seq2[[j2]]$gen=i
        }
        } # End if-else statement
    } # End loop along sequences (j)
    #assign new name to comp mut seq
    names(new_seq1)=rep(i, length(new_seq1))
    names(new_seq2)=rep(i, length(new_seq2))
    
    seqs = append(seqs,new_seq1)
    seqs = append(seqs,new_seq2)
    
    
    i=i+1
  } # End while loop
  
  if(length(unique(names(seqs))) < ngen) {
    print("Achieved too large a population before number of specified generations. Try changing some of the parameters.")
    stop()
  }
  # Sample at regular intervals
  sampling_times = seq(14, ngen, trunc(ngen/4))
  if(last(sampling_times) != ngen) {
    sampling_times = c(sampling_times, ngen)
  }
  
  seqs2 = seqs[names(seqs) %in% sampling_times]
  
    seqs2 = lapply(seqs2, "[[", 1)
    seqs2 = lapply(seqs2, unlist, use.names=FALSE)
    for (i in seq_along(seqs2)) {
      .gen=names(seqs2)[i]
      names(seqs2)[i] = paste("isolate", .gen, i, sep="_")
    }
  
  return(seqs2)
}

#applying function to make sequences
seqs=sim_fct(starting_seq, ngen, n, r, network, network_rates)
#writing alignment
write.dna(seqs, paste0('coevol_sims.fasta'), format="fasta")#name for each sim
#creating .csv file with expected coevolving sites
allsites = data.frame(Value = network)
write.table(allsites,paste0('coevol_sites.csv'),row.names = F, col.names = F)


