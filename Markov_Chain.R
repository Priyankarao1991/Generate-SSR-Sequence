
##Markov Chain
library(stringdist)
data <- read.csv("/Users/chaitrahande/Downloads/ssr_data_final.csv",header = FALSE,sep = ",")
n = nrow(data)
data_k_matrix = data[,1:31]
# transition probability matrix
trans.matrix <- function(X)
{
  mat <- table( c(X[,-ncol(X)]), c(X[,-1]) )
  #mat <- mat / rowSums(mat) 
  mat = sweep(mat,MARGIN=1,FUN="/",STATS=rowSums(mat))
  return(mat)
}
data_NA_matrix= as.matrix(data_k_matrix)
is.na(data_NA_matrix) <-data_NA_matrix==''

k_matrix=trans.matrix(data_NA_matrix)

letter_list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

# initial probability distribution( for the first position of the sequence) 
initial_Prob <- c(0.529914530,0,0.034188034,0,0,0.004273504,0,0,0,0.004273504,0,0.004273504,0,0,0,0.008547009,0.414529915,0,0,0); 
names(initial_Prob) <-letter_list

# number of sequences to generate
Nsim = 1000;
num_row = Nsim;
num_col=ncol(data)
length_list=c(25:31);

y=matrix(NA, nrow = num_row, ncol = num_col)
dist_matrix=matrix(0, nrow = Nsim, ncol = n)
result_seq = character()
z =character();
# Create a function to generate the sequence.

SSRSeq <- function(lengthSeq,alphabet,initial_Prob,ProbMatrix){
  
  # output sequence
  outputSeq <- rep(NA,lengthSeq)
  # generate first letter
  outputSeq[1]  <- sample(alphabet,1,prob=initial_Prob) 
  #generate rest
  for(i in 2:length(outputSeq)){
    prevSeq <- outputSeq[i-1]    
    currentProb <- ProbMatrix[prevSeq,]
    outputSeq[i] <- sample(alphabet,1,prob=currentProb)
  } 
  
  return(outputSeq)
}

j=1;
while(j<(Nsim+1)){
  
  length_pick = sample(length_list,1)
  df1<- SSRSeq(length_pick,letter_list,initial_Prob,k_matrix)
  y[j,1:length(df1)]=df1
  result_seq = paste(df1,collapse = "")
  z=c(z,result_seq)
  dist_matrix[j,] = stringdist(result_seq,data$V32,method = "lv")
  result_seq =character();
  j=j+1;
}







######Pattern frequecy
library(GrpString)
z1=CommonPatt(z, low = 30)

# Convert the last column in data to chanracter
data[32] <- lapply(data[ 32], as.character)
z2=CommonPatt(data$V32, low = 30)

#  HISTOGRAM OF THE PATTERN for generated and given data
par(mfrow=c(1,2))
#layout(1,2)
hist(z1$Freq_grp,breaks=5,main="Given SSR:Distribution of Patterns",xlab="Frequency",freq=FALSE)
lines(density(z1$Freq_grp),col="blue",lwd=2)

hist(z2$Freq_grp,breaks=5,main="Artificial SSR:Distribution of Patterns",xlab="Frequency",freq=FALSE)
lines(density(z2$Freq_grp),col="blue",lwd=2)












