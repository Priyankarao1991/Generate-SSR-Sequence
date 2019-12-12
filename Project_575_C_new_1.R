library(dplyr)
library(stringdist)
library(rowr)
library(seqinr)
library(msa)
library(Biostrings)
library(BiocGenerics)
library(ggplot2)
library(plot.matrix)
data <- read.csv("C:\\PhD\\Fall 2019\\Data_Science_575\\Project_Material\\ssr_data_final.csv",header = FALSE,sep = ",")

#### Below is the table of characters for each position
w1 = table(data$V1)
w2 = table(data$V2)
w3 = table(data$V3)
w4 = table(data$V4)
w5 = table(data$V5)
w6 = table(data$V6)
w7 = table(data$V7)
w8 = table(data$V8)
w9 = table(data$V9)
w10 = table(data$V10)
w11 = table(data$V11)
w12 = table(data$V12)
w13 = table(data$V13)
w14 = table(data$V14)
w15 = table(data$V15)
w16 = table(data$V16)
w17 = table(data$V17)
w18 = table(data$V18)
w19 = table(data$V19)
w20 = table(data$V20)
w21 = table(data$V21)
w22 = table(data$V22)
w23 = table(data$V23)
w24 = table(data$V24)
w25 = table(data$V25)
w26 = table(data$V26)
w27 = table(data$V27)
w28 = table(data$V28)
w29 = table(data$V29)
w30 = table(data$V30)
w31 = table(data$V31)

#### Storing the the frequencies colum of the table in even_indexes
even_indexes<-seq(2,62,2) # change later
n = nrow(data)
ProbMatrix = rowr::cbind.fill(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,w23,w24,w25,w26,w27,w28,w29,w30,w31, fill = NA)
ProbMatrix[even_indexes] <- lapply(ProbMatrix[ even_indexes], as.character)
ProbMatrix[is.na(ProbMatrix)] <- 0
ProbMatrix[even_indexes] <- lapply(ProbMatrix[even_indexes], function(x) as.numeric(as.character(x)))
ProbMatrix[even_indexes] = (ProbMatrix[even_indexes])/n  ### Calculating the probability for each character

#### Applying the inverse transform method

  pos = function(p,u) {

    prob=dplyr::pull(p,Freq)
    letter=dplyr::pull(p,Var1)
    letter=as.character(letter)
    if(u < prob[1]){
      return(letter[1])
    }
    for(idx in 2:length(prob)) {
      if(sum(prob[1:(idx-1)]) < u && u < sum(prob[1:idx]) ) {
        return(letter[idx])
      }
    }
    
}

Nsim=1000;  # desired number of sequences to be generated
num_row = Nsim;
num_col = ncol(data);
result_seq = character();
length_list = c(25:31); # choosing the length of each sequence
y=matrix(NA, nrow = num_row, ncol = num_col) # storing generated sequence as a matrix
dist_matrix = matrix(0, nrow = Nsim, ncol = n) # storing the edit distance as a matrix
result_seq_1 = character()  # storing the generated characters of sequence as a vector
z = character();  # storing the concatenated the charaters as a single sequence
for( i in 1:Nsim){
  k=1;
  length_pick = sample(length_list,1)
  for (j in 1: length_pick){
    u = runif(1)
    p=ProbMatrix[,k:(k+1)]
    p <-p[order(-p$Freq),]
    result= pos(p,u)
    k=k+2;
    result_seq = c(result_seq,result)
    result_seq_1 = result_seq
  }
  
  y[i,1:length(result_seq)]=result_seq
  result_seq_1 = paste(result_seq,collapse = "")
  z=c(z,result_seq_1)
 dist_matrix[i,] = stringdist(result_seq_1,data$V32,method = "lv") # jaccard
  result_seq =character();
  result_seq_1 = character()
}

z

######################################################  Statistical Plots

# get min value for each row in distance matrix

min_dist=apply(dist_matrix, 1, FUN=min)

mean(min_dist)

########## Percent Identity Sequence

########## PID matrix for Given vs Artificial SSRs
z1 = c()
pid_matrix = matrix(0, nrow = 234, ncol = 10)
for (j in 1:10) {
for(i in 1:234){
  palign6 <- pairwiseAlignment(z[j], data[i,32],substitutionMatrix = NULL)
  #palign6
  z1[i] = pid(palign6, type = "PID1")
}
  pid_matrix[,j] = z1
}


######## PID for given Vs Given
z1 = c()
pid_matrix = matrix(0, nrow = 10, ncol = 10)
for (j in 1:10) {
  for(i in 1:10){
    palign6 <- pairwiseAlignment(data[j,32], data[i,32],substitutionMatrix = NULL)
    #palign6
    z1[i] = pid(palign6, type = "PID1")
  }
  pid_matrix[,j] = z1
}

##### Plots for Given vs Given

plot(pid_matrix[1:10,], xlab = "Given SSRs", ylab = "Given SSRs",col = topo.colors,breaks = c(50, 60, 70, 80, 90, 100), digits = 4, cex = 0.5)

##### Plots for Artifical vs Given

plot(pid_matrix[1:10,], xlab = "Artifical SSRs", ylab = "Given SSRs",col = topo.colors,breaks = c(50, 60, 70, 80, 90, 100), digits = 4, cex = 0.5)

plot(pid_matrix, xlab = "Artifical SSRs", ylab = "Given SSRs",col = topo.colors,breaks = c(50, 60, 70, 80, 90, 100))


plot(dist_matrix[1:10,1:10], xlab = "Given SSRs",ylab = "Artifical SSRs")
############################# PID Plot

