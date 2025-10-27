rm(list=ls()) #clean, clc, close all
# ***************
# R version 4.4.2 / win
# author: Nikola Jordánová
# *************

# Path
setwd('V:/MPA-PRG/exercise_07') # set working directory

library(Biostrings)

# Task 1
Score <- function(start_idxs, sequences, motif_length){
  
  count_matrix <- matrix(0, nrow = 4, ncol = motif_length)
  rownames(count_matrix) <- c("A", "C", "G", "T")
  
  for (i in 1:length(sequences)) {
    subseq_i <- as.character(subseq(sequences[[i]], start = start_idxs[i], width = motif_length))
    
    chars <- unlist(strsplit(subseq_i, ""))
    for (j in 1:motif_length) {
      base <- chars[j]
      if (base %in% c("A", "C", "G", "T")) {
        count_matrix[base, j] <- count_matrix[base, j] + 1
      }
    }
  }
  
  score <- sum(apply(count_matrix, 2, max))
  
  return(score)
  
}

seqs <- readDNAStringSet("seq_score.fasta")
starts <- c(1, 1, 1, 1, 1)
motif_length <- 5

Score(starts, seqs, motif_length)


# Task 2
NextLeaf <- function(s, t, k){
  # s - array of starting indexes
  # t - number of sequences
  # k - k = n - l + 1, where n is length of sequences and l is motif length
  for (i in t:1){
    if (s[i] < k){
      s[i] <- s[i] + 1
      if (i < t) {
        s[(i+1):t] <- 1
      }
      return(s)
    }
    s[i] <- 1
  }
  return(NULL)
}


# Task 3
BFMotifSearch <- function(DNA, t, n, l){
  # DNA - DNAStringSet object of sequences
  # t - number of sequences
  # n - length of each sequence
  # l - motif length
  
  s <- rep(1, t)
  bestScore <- Score(s, DNA, l)
  bestMotif <- s
  
  while (TRUE){
    s <- NextLeaf(s, t, n - l + 1)
    
    if (is.null(s)) {
      return(bestMotif)
    }
    
    currentScore <- Score(s, DNA, l)
    if (currentScore > bestScore){
      bestScore <- currentScore
      bestMotif <- s
    }
  }
}


DNA <- readDNAStringSet("seq_motif.fasta")

t <- length(DNA)      
n <- width(DNA)[1] 
l <- 3


bestMotif <- BFMotifSearch(DNA, t, n, l)
print(bestMotif)

# Task 4
NextVertex <- function(s, i, t, k){
  # s - array of starting indexes s = (s1 s2 … st), where t is the number of sequences
  # i - level of vertex
  # t - number of sequences
  # k - k = n - l + 1, where n is length of sequences and l is motif length
  
  if (i < t){
    s[i+1] <- 1
    i <- i + 1
    return (list(s=s, i=i))
  }
  else{
    for (j in t:1){
      if (s[j] < k){
        s[j] <- s[j] + 1
        return (list(s=s, i=j))
      }
    }
  }
  return (list(s=s, i=0))
  
}


  
s <- c(1, 1, 1)
i <- 1 
t <- 3 
k <- 5


result <- NextVertex(s, i, t, k)
print(result)


# Task 5
ByPass <- function(s, i, t, k){
  # s - s = (s1 s2 … st); an array of starting indexes, where t is the number of sequences
  # i - level of vertex
  # t - number of DNA sequences
  # k - k = n - l + 1, where n is length of DNA sequences and l is motif length
  
  for (j in i:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return (list(s=s,i=j))
    }
  }
  return (list(s=s, i=0))
}


s <- c(1, 1, 1) 
i <- 3 
t <- 3       
k <- 5  


result <- ByPass(s, i, t, k)
print(result)

# Task 6
BBMotifSearch <- function(DNA, t, n, l){
  # DNA - DNAStringSet object of sequences
  # t - number of sequences
  # n - length of each sequence
  # l - motif length
  
  s <- rep(1, t)
  bestScore <- 0
  i <- 1
  while (i>0){
    if (i<t){
      optimisticScore <- Score(s[1:i], DNA[1:i], l) + (t - i) * l
      if (optimisticScore < bestScore){
        bypass <- ByPass(s, i, t, n - l + 1)
        s <- bypass$s
        i <- bypass$i
      }
      else{
        res <- NextVertex(s, i, t, n - l + 1)
        s <- res$s
        i <- res$i
      }
    }
    else{
      if (Score(s[1:t], DNA[1:t], l) > bestScore){
        bestScore <- Score(s[1:t], DNA[1:t], l)
        bestMotif <- s
      }
      res <- NextVertex(s, i, t, n - l + 1)
      s <- res$s
      i <- res$i
    }
  }
  return (bestMotif)
  
}

DNA <- readDNAStringSet("seq_motif.fasta")

t <- length(DNA)
n <- width(DNA)[1]   # délka sekvence
l <- 3               # délka motivu

bestMotif <- BBMotifSearch(DNA, t, n, l)
print(bestMotif)






