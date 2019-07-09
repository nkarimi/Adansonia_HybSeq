# inspired from evobiR::CalcD but requires package "seqinr" only
calcPatternCounts = function(alignment = "alignment.fasta",
                 ambig="D",           # options are D R I
                 align.format='fasta'){
  # read alignment and make it a matrix
  alg <- read.alignment(alignment, format = align.format, forceToLower=T) # read.alignment from seqinr package

  # check that there are 4 taxa, and all have same alignment length
  # sometimes, 1 or more taxa might be missing a gene of interest
  badalign = FALSE
  if (length(alg$seq) != 4){
    badalign = TRUE
  } else {
    if (nchar(alg$seq[[1]]) != nchar(alg$seq[[2]]) ||
        nchar(alg$seq[[1]]) != nchar(alg$seq[[3]]) ||
        nchar(alg$seq[[1]]) != nchar(alg$seq[[4]])){
      badalign = TRUE
    }
  }
  if (badalign) {
    return(c(NA,NA,NA))
  }
  # convert the alignment to a matrix
  alignment <- matrix(, length(alg$nam), nchar(alg$seq[[1]]))             # empty matrix
  for (i in 1:length(alg$nam)){
    alignment[i, ] <- unlist(strsplit(alg$seq[[i]], ""))                  #  fill in the matrix
  }
  #### This section is being added to deal with reccurent
  #### Requests to deal with ambiguity in sequence data
  # R A or G
  # Y C or T
  # S G or C
  # W A or T
  # K G or T
  # M A or C

  ## First we deal with the situation where the user
  ## wishes to simply drop ambig sites
  if(ambig == "D"){
    target <- c("a","c","g","t")
    keep <- vector(, ncol(alignment))
    for(i in 1:ncol(alignment)){
      keep[i] <- all(alignment[,i] %in% target)
    }
    alignment <- alignment[,keep, drop=F]
  }

  ## Next we deal with the situation where users want to
  ## randomly resolve ambigous sites
  if(ambig == "R"){
    # I still want to limit sites so we first drop
    # those sites that look like 3 or 4 possibilities
    target <- c("a", "c", "g", "t", "r",
                "y", "s", "w", "k", "m")
    keep <- vector(, ncol(alignment))
    for(i in 1:ncol(alignment)){
      keep[i] <- all(alignment[,i] %in% target)
    }
    alignment <- alignment[,keep, drop=F]

    # this function will be applied to each site in our data
    # it resolves ambiguities randomly
    resolver <- function(x){
      if(x=="r") z <- sample(c("a", "g"), 1)
      if(x=="y") z <- sample(c("c", "t"), 1)
      if(x=="s") z <- sample(c("g", "c"), 1)
      if(x=="w") z <- sample(c("a", "t"), 1)
      if(x=="k") z <- sample(c("g", "t"), 1)
      if(x=="m") z <- sample(c("a", "c"), 1)
      if(x %in% c("a", "c", "g", "t")) z <- x
      return(z)
    }
    alignment <- apply(alignment, c(1,2), resolver)
  }

  # now calculate the number of sites with pattern ABBA or BABA
  abba = 0
  baba = 0
  nsites = ncol(alignment)
  if (nsites>0){ for (i in 1:nsites){
    if (length(unique(alignment[, i])) == 2){                        # unique(c(p1,p2,p3,o))==2 aka biallelic
      if (alignment[1, i] != alignment[2, i]){                       # p1 != p2   aka different resolutions in p1 and p2
        if (alignment[4, i] != alignment[3, i]){                     # o  != p3   Durand says "less likely pattern due to seq. errors
          if (alignment[3, i] == alignment[1, i]) {baba <- baba + 1} # add to the count of baba sites
          if (alignment[2, i] == alignment[3, i]) {abba <- abba + 1} # add to the count of abba sites
        }
      }
    }
  }}

  #cat("\nSites in alignment =", ncol(alignment))
  #cat("\nNumber of sites with ABBA pattern =", abba)
  #cat("\nNumber of sites with BABA pattern =", baba,"\n\n")
  return(c(nsites,abba,baba))
}

# calculate the D statistic given a list of
# x=ABBA and y=BABA counts across several genes
calcD = function(x, y){
  sx = sum(x); sy = sum(y)
  return((sx - sy)/(sx + sy))
}

# calculate the standard deviation of the D statistics
# when re-sampling loci (bootstrap)
calcDsd = function(x, y, Nbootstrap=500){
  nloci = length(x)
  bloci = sample.int(nloci, nloci*Nbootstrap, replace=T)
  # next: bootstrap ABBA & BABA counts --all of them: 1 column = 1 bootstrap sample
  bx = x[bloci] # matrix(x[bloci], nloci, Nbootstrap)
  by = y[bloci]
  bsx = .colSums(bx, nloci, Nbootstrap) # summing per column
  bsy = .colSums(by, nloci, Nbootstrap)
  bD = (bsx - bsy)/(bsx + bsy)
  return(sd(bD))
}
