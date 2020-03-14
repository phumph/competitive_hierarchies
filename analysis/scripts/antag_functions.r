#### Functions associated wtih Humphrey et al.
#### last updated: 21 June 2016

#### HEADER ####
# function definitions
c.calc <- function(cmat)
{
  cres <- data.frame(strain.id = row.names(cmat), nO = NA, nD = NA, cO = NA, cD = NA, cW = NA)
  
  # calculate the number of interactions in the defensive and offense positions:
  for(i in 1:length(cmat[,1]))
  {
    # offense
    cres[i,'nO'] <- length(cmat[,i]) - length(is.na(cmat[,i])[is.na(cmat[,i]) == TRUE]) - 1
    # defense
    cres[i,'nD'] <- length(cmat[i,]) - length(is.na(cmat[i,])[is.na(cmat[i,]) == TRUE]) - 1
  }
  
  cres[,'cD.raw'] <- rowSums(cmat, na.rm=TRUE)
  cres[,'cO.raw'] <- colSums(cmat,na.rm=TRUE)
  
  cres[,'cO'] <- round(with(cres,cO.raw/nO),4)
  cres[,'cD'] <- round(with(cres,(nD-cD.raw)/nD),4)
  
  cres[,'cW'] <- round(with(cres, (cO - 1 + cD)),4)
  
  return(cres)
}

# takes matrix same dimensions as cmat as n c.calc() above
i.calc <- function(imat)
{
  require(dplyr)
  
  ires <- data.frame(strain.id = row.names(imat), nI = NA, nR = NA, cR = NA,cT = NA,iW = NA)
  
  for(i in 1:length(imat[,1]))
  {
    # interactions as invader
    ires[i,'nI'] <- length(imat[,i]) - length(is.na(imat[,i])[is.na(imat[,i]) == TRUE]) - 1
    # interactions as defender (this is irrelevant to this section)
    ires[i,'nR'] <- length(imat[i,]) - length(is.na(imat[i,])[is.na(imat[i,]) == TRUE]) - 1
  }
  
  # sum number of strains each kills
  ires[,'cT.raw'] <- colSums(imat,na.rm=TRUE)
  
  # sum number of strains each is killed by (subtract from number of cT>0 strains to get cR):
  # define submatrix with strains with cT>0 as cols:
  imat2 <- imat[,names(imat) %in% filter(ires,cT.raw > 0)[,'strain.id']]
  ires[,'cR.raw'] <- length(imat2[1,]) - rowSums(imat2, na.rm=TRUE) # this needs to be bounded by the number of strains that have nT > 0
  
  ires[,'cR'] <- round(with(ires,cR.raw/length(imat2[1,])),4)
  ires[,'cT'] <- round(with(ires,cT.raw/nI),4)
  
  ires[,'iW'] <- round(with(ires, (cT - 1 + cR)),4)
  
  return(ires)
}


#### stacker() #### takes 96-well plate input and vertically stacks everything either row- or col-wise
stacker <- function(z, row = TRUE)
{
  if (row == TRUE)
  {
    zz <- t(z[1,])
    for(i in 2:length(z[,1]))
    {
      zz <- rbind(zz,t(z[i,]))
    }  
  }
  else
  {
    zz <- paste(z[,1])
    for(i in 2:length(z[1,]))
    {
      zz <- c(zz,paste(z[,i]))
    }  
    t(t(zz))
  }
  return(zz)
}



# TO DO:
# v 1. Find final matrixes for exploitative and interference competition
# √ 2. Re-calculate all components of competitiveness to cross-check approach.
# 3. Re-export Heatmap figure and PCA


#### DATA EXPLORATION ####
c1 <- read.table("c_matrix.txt",row.names=1,header=T,"\t")
row.names(c1) <- names(c1)
diag(c1) <- 0
c1 <- c1/2

i1 <- read.table("i_matrix.txt",row.names=1,header=T,"\t")
row.names(i1) <- names(i1)
diag(i1) <- 0

cres1 <- c.calc(c1)
cres1[,'phylo'] <- c(1:40) # add phylo position to this matrix and sort final combined data by it
cres1[,'clade'] <- c(rep('Psyr',26),rep('Pflu',14))
ires1 <- i.calc(i1)
res1 <- arrange(merge(cres1,ires1),phylo)



