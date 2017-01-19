options(stringsAsFactors=FALSE)

library(GenSA)
library(deconstructSigs)
library(quadprog)

ReadMutationalProfile <- function(file) {
  #Read data as a matrix
  return(data.matrix(read.table(file, sep='\t', header=T, row.names=1, check.names=F)))
}

QPsig <- function(tumour.ref, signatures.ref) {
  # This code was addapted from Lynch, F1000 Research, 2016.
  # The QPsig function is designed to take input in similar format to the
  # whichSignatures function in the deconstructSigs library for easy interchange
  # between the two. The output is limited to the equivalent of the `weights'
  # slot in the output of whichSignatures.

  # we normalize the observations so that they sum to 1
  obs<-as.numeric(tumour.ref/sum(tumour.ref))

  # to allow use of the deconstructSigs objects we convert to matrices
  signatures.ref<-as.matrix(signatures.ref)

  # we use the factorized version of solve.QP -
  # see the helpfile of that function for details of the required form
  # otherwise we would set Dmat = signatures.ref %*% t(signatures.ref) as indicated
  # in the article
  Rinverse <- backsolve(chol(signatures.ref %*% t(signatures.ref)),
                        diag(dim(signatures.ref)[1]))

  # we also need to define the linear part of the function that we are minimizing
  dvec <- (obs) %*% t(signatures.ref)

  # we have one constraint that the sum of weights is equal to 1
  # we have more constraints that each weight is positive
  Amat <- cbind(rep(1,dim(Rinverse)[1]), diag(dim(Rinverse)[1]))
  bvec <- c(1,rep(0,dim(Rinverse)[1]))

  # we now call the solve.QP function from the quadprog library
  myQP<-quadprog::solve.QP(Dmat = Rinverse, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1,
                           factorized = TRUE)
  return(myQP$solution)
}

DecompositionError <- function(observation, weights, signatures) {
  #Squared root of squared sum of errors
  #Remove NA by replacing it with 0
  weights.noNA = weights
  weights.noNA[is.na(weights.noNA)] = 0
  return(sqrt(sum((observation - signatures %*% weights.noNA)^2)))
}

RunGenSA <- function(spectrum, signatures, control) {
  #Wrapper for GenSA function to run simulated annealing
  K = ncol(signatures)
  DecompositionError.local <- function(weights) {
    #objective function to compute estmation error
    estimate = signatures %*% weights
    return(sqrt(sum((spectrum - (estimate / sum(estimate)))^2)))
  }

  #run GenSA
  sa = GenSA(par=NULL, lower=rep(0.0,K), upper=rep(1.0,K), fn=DecompositionError.local, control=control)
  return(sa)
}


RunGenSATrialsWithDesiredError <- function(spectrum, signatures, N, desired.error) {
  #run N trials until reaching expected error to sample landscape of almost optimal solutions
  error = rep(0.0, N)
  weights = matrix(0.0, ncol=N, nrow=ncol(signatures))
  rownames(weights) = colnames(signatures)
  colnames(weights) = paste0('Trial_', seq(N))
  #repeat the simulation N times
  for(i in seq(N)) {
    #reduce maximum number of iterations of the algorithm (the algorithm restarts anyway after 1282 iterations)
    #reduce temperature; we want to explore params around optimal solution (it is a convex function so there is no other local minimum)
    #stop when there is no improvement in 1000 steps as this is simple function
    #stop when reaching desired error (higher then in optimal solution)
    sa = RunGenSA(spectrum, signatures, list(maxit=1000, temperature=100,
      threshold.stop=desired.error, nb.stop.improvement=1000, simple.function=T))
    #record results
    error[i] = sa$value
    weights[,i] = sa$par/sum(sa$par) #normalize weights to 1
  }
  return(list(error=error, weights=weights))
}

PlotCombinedWeights <- function(gsat, gsa.weights, qp.weights, deSig.weights, nmf.weights, spectrum, cosmic, title) {
  #plots results of different methods on one plot
  RSSE <- function(weights) {
    DecompositionError(spectrum,weights,cosmic) #squared root of squared sum of errors
  }
  #find max weight to set y-axis limits
  weight.max = max(gsat$weights, gsa.weights, qp.weights, deSig.weights, nmf.weights, na.rm=T)
  #generate subtitle listing error of each method
  subtitle = sprintf("GenSA+error(median) %.5f, GenSA %.5f, QP %.5f, deconstructSigs %.5f, NMF %.5f",
    median(gsat$error), RSSE(gsa.weights), RSSE(qp.weights), RSSE(deSig.weights), RSSE(nmf.weights))
  #boxplot of GenSA trials with increased error
  names_sigs = sapply(strsplit(rownames(gsat$weights), ' '), '[', 2)
  boxplot(t(gsat$weights), names=names_sigs, col='gray90', xlab='Signatures', ylab='Weights', main=title, sub=subtitle, las=2, ylim=c(0,1.05*weight.max))
  #show optimal results of different methods
  points(nmf.weights, pch=6, col='darkorange')
  points(deSig.weights, pch=2, col='limegreen')
  points(qp.weights, pch=13, col='deepskyblue')
  points(gsa.weights, pch=1, col='firebrick1')
  legend('topright', c('GenSA w/ error','GenSA','QP','deconstructSigs','NMF'),
    col=c('gray90','firebrick1','deepskyblue','limegreen','darkorange'), pch=c(15,1,13,2,6))
}

################################################################################

sink('sa_GenSA_BRCA.output')

#BRCA signatures
selected_sigs = c(1,2,3,5,6,8,13,17,18,20,26,30)

### Read data ###
#samples are in columns, mutation types in rows
#cosmic matrix was downloaded from COSMIC website (96x30)
cosmic = ReadMutationalProfile('data/cosmic_signatures_probabilities.txt')
cosmic = cosmic[,selected_sigs]
#spectra are precomputed mutational profiles for all 560 patiens and for a combined dataset (96x561)
#other data can be provided
spectra = ReadMutationalProfile('data/mutation_spectra_BRCA-EU.txt')
#match mutation type order to cosmic order as in file
spectra = spectra[match(rownames(cosmic),rownames(spectra)), ,drop=F]

### Other methods first ###

### Quadratic programming
#run QPsig for each column of spectra (it expects signatures in deconstructSigs form)
#gets a 30x561 matrix of signature weights for each sample
#some weights are negative, but very very close to 0 (e.g. -1e-20)
print('Quadratic programming')
print(date())
QP = apply(spectra, 2, QPsig, t(cosmic))
print(date())
rownames(QP) = colnames(cosmic)
write.csv(QP, file='weights_BRCA/QP.csv')

### deconstructSigs
#Cosmic and signatures.cosmic matrices contains exactly the same data, but are differently organized.
#They are just differently ordered, so watch out when passing parameters.
#Match mutations order to that in package signatures!
#set 'signature.cutoff' 0.0, so we do not miss any signature
#ncol(spectra)
print("deconstructSigs")
print(date())
deSig = do.call(cbind, lapply(seq(ncol(spectra)), function(i) {
    ds = whichSignatures(tumor.ref=as.data.frame(t(spectra[match(colnames(signatures.cosmic),rownames(spectra)),i,drop=F])),
      signatures.ref=signatures.cosmic, contexts.needed=T, signature.cutoff=0.0, associated=row.names(signatures.cosmic)[selected_sigs])
    return(t(ds$weights))
  }))
deSig = deSig[selected_sigs,] #deconstructSigs still output all signatures, so we have to limit it again
print(date())
rownames(deSig) = colnames(cosmic)
write.csv(deSig, file='weights_BRCA/deSig.csv')

### NMF
print("NMF - read and process data")
print(date())
#read the supplementary table 21 by S. Nik-Zainal (Nature, 2016)
#containing weights of signatures identified by NMF
nmf.data = read.csv('data/Table_SignatureContribution__SupplTab21.csv', check.names=F)
nmf.data_samples = nmf.data[,'Sample Name']
nmf.data = as.matrix(nmf.data[,2:13]) #remove first and last column (not needed) and transform to matrix
rownames(nmf.data) = nmf.data_samples
nmf.data = rbind(nmf.data, All=colSums(nmf.data)) #add summary row and name it 'All' (the same as in spectra)
nmf.data = nmf.data[match(colnames(spectra),rownames(nmf.data)),] #match sample names to that in file in spectra
nmf.data = apply(nmf.data, 2, "/", rowSums(nmf.data))
#create a new matrix with weights for all signatures and update based on nmf.data
NMF = matrix(NA, nrow=ncol(cosmic), ncol=nrow(nmf.data),
  dimnames=list(colnames(cosmic), rownames(nmf.data)))
NMF[match(colnames(nmf.data), colnames(cosmic)),] = t(nmf.data)
write.csv(NMF, file='weights_BRCA/NMF.csv')
print(date())


### Simulated Annealing ###
### Generalized Simulated Annealing by Y. Xiang (The R Journal, 2013) implemented in R package 'GenSA'

#run long SA to determine best solution using default params
set.seed(1234)
print("Simulated annealing - long run")
print(date())
i = 0 #just to print status
GSA = apply(spectra, 2, function(spectrum) {
    i <<- i+1 #just to print status
    print(c(colnames(spectra)[i],date())) #just to print status
    gsa = RunGenSA(spectrum, cosmic, control=list(simple.function=T))
    gsa$par/sum(gsa$par) #normalize to 1
  })
print(date())
rownames(GSA) = colnames(cosmic)
write.csv(GSA, file='weights_BRCA/GenSA.csv')


### Test GenSA  different levels of allowed error (only on combined data)
print("Simulated annealing - testing different errors")
print(date())
set.seed(1234)
pdf('plots/boxplot_BRCA_all_diff_errors.pdf', height=11, width=8.5)
par(mfrow=c(3,1))
error.GSA.all = DecompositionError(spectra[,'All'], GSA[,'All'], cosmic)
for(error.increase in c(0.001,0.01,0.05,0.1)) {
  print(c(error.increase, date()))
  gsat = RunGenSATrialsWithDesiredError(spectra[,'All'], cosmic, 1000, (1+error.increase)*error.GSA.all)

  PlotCombinedWeights(gsat, GSA[,'All'], QP[,'All'], deSig[,'All'], NMF[,'All'], spectra[,'All'], cosmic,
    paste0('Boxplots based on SA for optimal GenSA error * ',1+error.increase,' (',100*error.increase,'%)','\n','All samples (12 sigs: 1,2,3,5,6,8,13,17,18,20,26,30)'))
}
par(mfrow=c(1,1))
dev.off()
print(date())


### Simulated annealing with increased error by 1% and 5% for all samples
for(error.increase in c(0.01,0.05)) {
  print(sprintf("Simulated annealing with error increase by %.2f fold", 1+error.increase))
  dir_weight = paste0('weights_BRCA/GenSA_trials_error_',error.increase)
  dir.create(dir_weight, showWarnings = FALSE)
  print(date())
  set.seed(1234)
  pdf(paste0('plots/boxplot_BRCA_samples_error_',error.increase,'.pdf'), height=11, width=8.5)
  par(mfrow=c(3,1))
  for(i in seq(ncol(spectra))) {
    sample = colnames(spectra)[i]
    print(c(sample,date()))
    error.GSA = DecompositionError(spectra[,i], GSA[,i], cosmic)
    gsat = RunGenSATrialsWithDesiredError(spectra[,i], cosmic, 100, (1+error.increase)*error.GSA)
    PlotCombinedWeights(gsat, GSA[,i], QP[,i], deSig[,i], NMF[,i], spectra[,i], cosmic,
      paste0(sample,'(optimal GSA error * ',1+error.increase,')\n','(12 sigs: 1,2,3,5,6,8,13,17,18,20,26,30)'))
    write.csv(gsat$weights, paste0(dir_weight,'/',sample,'.csv'))
  }
  par(mfrow=c(1,1))
  dev.off()
  print(date())
}


sink()
