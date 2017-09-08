
#########################################################
# This program is part of the SNN-Cliq method           #
# Contact Chen Xu at UNC-Charlotte for more information.# 
#########################################################
#----- example of use------#
#data<-read.table(infile, header=TRUE, sep="\t", row.names=1);
#data<-log2(data+1)
#source('SNN.R')
#SNN(data, edge_file, k=3, distance='euclidean')
#--------------------------#

SNN<-function(data, outfile, k, distance){
  
  if(missing(data)){
    stop(paste("Input data missing.",help,sep="\n"))
  }
  if(missing(outfile)){
    stop(paste("Output file name missing.",help,sep="\n"))
  }
  if(missing(k)){
    k=3
  }
  if(missing(distance)){
    distance<-"euclidean"  # other distance options refer to dist() in R
  }
  m<-as.data.frame(data)
  numSpl<-dim(data)[1]
  m<-dist(data, distance, diag=TRUE, upper=TRUE)
  x<-as.matrix(m)
  IDX<-t(apply(x,1,order)[1:k,]) # knn list
  
  edges<-list()              # SNN graph
  for (i in 1:numSpl){
    j<-i
    while (j<numSpl){
      j<-j+1
      shared<-intersect(IDX[i,], IDX[j,])
      if(length(shared)>0){			
        s<-k-0.5*(match(shared, IDX[i,])+match(shared, IDX[j,]))
        strength<-max(s)
        if (strength>0)
          edges<-rbind(edges, c(i,j,strength))
      }				
    }
  }
  write.table(edges, outfile, quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
}

