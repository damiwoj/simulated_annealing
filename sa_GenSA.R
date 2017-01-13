sink('sa_GenSA.output')
source('~/supercoiling/Software/R_mylib/preloadsettings.R')

library(GenSA)
source('QPsig.R')
library(deconstructSigs)
library(quadprog)

ReadMutationalProfile <- function(file) {
  return(data.matrix(read.table(file, sep='\t', header=T, row.names=1, check.names=F)))
}

cosmic = ReadMutationalProfile('../data/cosmic_signatures_probabilities.txt')
spectra = ReadMutationalProfile('../data/mutation_spectra_BRCA-EU.txt')
spectra = spectra[match(rownames(cosmic),rownames(spectra)), ,drop=F] #match mutations to cosmic order as in file

RunGenSA <- function(spectrum, signatures, control) {
  K = ncol(signatures)

  #objective function to compute estmation error
  DecompositionError <- function(weights) {
    #fn.call <<- fn.call + 1
    estimate = signatures %*% weights
    sse = sqrt(sum((spectrum - (estimate / sum(estimate)))^2))
    return(sse)
  }

  #run GenSA
  sa_out = GenSA(par=NULL, lower=rep(0.0,K), upper=rep(1.0,K),
    fn=DecompositionError,
    control=control)

  return(sa_out)
}

RunGenSATrialsWithIncreasedError <- function(spectrum, signatures, N, error.perc) {
  #run long SA to determine best solution using default params
  gsa_best = RunGenSA(spectrum, signatures, control=list(simple.function=T))

  #run N trials allowing higher solution error to sample landscape of almost optiomal solutions
  error = rep(0.0, N)
  weights = matrix(0.0, ncol=ncol(signatures), nrow=N)
  colnames(weights) = colnames(signatures)
  for(i in seq(N)) {
    #reduce maximum number of iterations of the algorithm (the algorithm restarts anyway after 1282 iterations)
    #reduce temperature; we want to explore params around optimal solution (it is a convex function so there is no other local minimum)
    #increase optimal error (found above) by 'error.perc' percent
    gsa = RunGenSA(spectrum, signatures, list(maxit=1282, temperature=100, threshold.stop=gsa_best$value*(1+error.perc), simple.function=T))
    error[i] = gsa$value
    weights[i,] = gsa$par/sum(gsa$par)
  }

  return(list(error=error, weights=weights, long=list(error=gsa_best$value, weights=gsa_best$par/sum(gsa_best$par))))
}

CombineBoxplot <- function(sat, qp.weights, deSig.weights, nmf.weights, spectrum, title) {
  SSE <- function(weights) {
    sum((spectrum - (cosmic %*% weights))^2)
  }
  weight.max = max(sat$weights, qp.weights, deSig.weights, nmf.weights)
  subtitle = sprintf("GenSA %.5f, GenSA+error(median) %.5f, QP %.5f, deconstructSigs %.5f, NMF %.5f", sat$long$error, median(sat$error), sqrt(SSE(qp.weights)), sqrt(SSE(deSig.weights)), sqrt(SSE(nmf.weights)) )
  boxplot(sat$weights, names=1:30, col='gray90', xlab='Signatures', ylab='Weights', main=title, sub=subtitle, las=2, ylim=c(0,1.05*weight.max))
  points(nmf.weights, pch=6, col='darkorange')
  points(deSig.weights, pch=5, col='limegreen')
  points(qp.weights, pch=2, col='deepskyblue')
  points(sat$long$weights, pch=1, col='firebrick1')
  legend('topright', c('GenSA w/ error','GenSA','QP','deconstructSigs','NMF (12 sigs)'), col=c('gray90','firebrick1','deepskyblue','limegreen','darkorange'), pch=c(15,1,2,5,6))
}

#Quadratic programming
QuProg = QPsig(tumour.ref=spectra[,1], signatures.ref=t(cosmic))

#run deconstructSigs
#match mutations order to that in package signatures!
deSig = whichSignatures(tumor.ref=as.data.frame(t(spectra[match(colnames(signatures.cosmic),rownames(spectra)),1,drop=F])),
  signatures.ref=signatures.cosmic, contexts.needed=T, signature.cutoff=0.0)

#NMF
nmf.data = read.csv('../data/Table_SignatureContribution__SupplTab21.csv', check.names=F)
nmf.data = rbind(nmf.data,c('All',colSums(nmf.data[,2:13]),NA)) #add summary row and name it 'All' (the same as in spectra)
rownames(nmf.data) = nmf.data[,1]
nmf.data = nmf.data[match(colnames(spectra),rownames(nmf.data)),2:13] #match sample names to that in file in spectra
#nmf.data_sum = rowSums(nmf.data) # DO NOT WORK!!!
nmf.data_sum = as.numeric(apply(nmf.data, 1, function(r) sum(as.numeric(r))))
nmf.data = apply(nmf.data, 2, function(col) as.numeric(col)/nmf.data_sum)
rownames(nmf.data) = colnames(spectra) #recover colnames
nmf.weights = matrix(0, nrow=nrow(nmf.data), ncol=ncol(cosmic))
colnames(nmf.weights) = colnames(cosmic)
rownames(nmf.weights) = rownames(nmf.data)
nmf.weights[,match(colnames(nmf.data), colnames(cosmic))] = nmf.data #match order of signature subset to that in cosmic

set.seed(1234)

print("Simulated annealing - testing different errors")
date()
sat0 = RunGenSATrialsWithIncreasedError(spectra[,1], cosmic, 1000, 0.001)
sat1 = RunGenSATrialsWithIncreasedError(spectra[,1], cosmic, 1000, 0.01)
sat2 = RunGenSATrialsWithIncreasedError(spectra[,1], cosmic, 1000, 0.05)
sat3 = RunGenSATrialsWithIncreasedError(spectra[,1], cosmic, 1000, 0.1)
date()
pdf('boxplot_all_diff_errors.pdf', height=11, width=8.5)
par(mfrow=c(2,1))
CombineBoxplot(sat0, QuProg, as.numeric(deSig$weights), nmf.weights['All',], spectra[,1], paste0('Boxplots based on SA for best error + ',100*0.001,'%\n','All samples'))
CombineBoxplot(sat1, QuProg, as.numeric(deSig$weights), nmf.weights['All',], spectra[,1], paste0('Boxplots based on SA for best error + ',100*0.01,'%\n','All samples'))
CombineBoxplot(sat2, QuProg, as.numeric(deSig$weights), nmf.weights['All',], spectra[,1], paste0('Boxplots based on SA for best error + ',100*0.05,'%\n','All samples'))
CombineBoxplot(sat3, QuProg, as.numeric(deSig$weights), nmf.weights['All',], spectra[,1], paste0('Boxplots based on SA for best error + ',100*0.1,'%\n','All samples'))
par(mfrow=c(1,1))
dev.off()


print("Simulated annealing with increased error by 1%")
#weights = list(GenSA = data.frame(), QP = data,frame(), decSigs = data.frame(), NMF = data.frame())
date()
set.seed(1234)
pdf('boxplot_samples_error_0.01.pdf', height=11, width=8.5)
par(mfrow=c(3,1))
for(i in seq(ncol(spectra))) {
  sample = colnames(spectra)[i]
  print(sample)
  qp = QPsig(tumour.ref=spectra[,i], signatures.ref=t(cosmic))
  ds = whichSignatures(tumor.ref=as.data.frame(t(spectra[match(colnames(signatures.cosmic),rownames(spectra)),i,drop=F])), #match mutations order to that in package signatures!
    signatures.ref=signatures.cosmic, contexts.needed=T, signature.cutoff=0.0)
  sat = RunGenSATrialsWithIncreasedError(spectra[,i], cosmic, 100, 0.01)
  CombineBoxplot(sat, qp, as.numeric(ds$weights), as.numeric(nmf.weights[sample,]), spectra[,i], colnames(spectra)[i])
  # #remember weights
  # weights$GenSA = rbind(weights$GenSA, sat$long$weights)
  # weights$QP = rbind(weights$QP, qp)
  # weights$decSigs = rbind(weights$decSigs, as.numeric(ds$weights))
  # weights$NMF = rbind(weights$NMF, as.numeric(nmf.weights[sample,]))
}
par(mfrow=c(1,1))
dev.off()
date()
# for(name in names(weights)) {
#   colnames(weights[[name]]) = colnames(cosmic)
#   rownames(weights[[name]]) = colnames(spectra)
#   write.csv(weights[[name]], file=paste0('weights_',name,'.csv'), quote=F)
# }


print("Simulated annealing with increased error by 5%")
date()
set.seed(1234)
pdf('boxplot_samples_error_0.05.pdf', height=11, width=8.5)
par(mfrow=c(3,1))
for(i in seq(ncol(spectra))) {
  print(colnames(spectra)[i])
  qp = QPsig(tumour.ref=spectra[,i], signatures.ref=t(cosmic))
  ds = whichSignatures(tumor.ref=as.data.frame(t(spectra[match(colnames(signatures.cosmic),rownames(spectra)),i,drop=F])), #match mutations order to that in package signatures!
    signatures.ref=signatures.cosmic, contexts.needed=T, signature.cutoff=0.0)
  sat = RunGenSATrialsWithIncreasedError(spectra[,i], cosmic, 100, 0.05)
  CombineBoxplot(sat, qp, as.numeric(ds$weights), as.numeric(nmf.weights[colnames(spectra)[i],]), spectra[,i], colnames(spectra)[i])
}
par(mfrow=c(1,1))
dev.off()
date()


sink()




# set.seed(1234)
# gsa = RunGenSA(spectra[,1], cosmic, list(maxit=1000, nb.stop.im4rovement=10, simple.function=T, verbose=F))
# list(error=gsa$value, fn.count=gsa$count, weights=gsa$par/sum(gsa$par))
# write.table(gsa$trace.mat, file='trace.mat.txt', quote=F, sep='\t')
# set.seed(1234)
# gsa = RunGenSA(spectra[,1], cosmic, list(maxit=1000, nb.stop.improvement=10, temperature=100, simple.function=T, verbose=F))
# list(error=gsa$value, fn.count=gsa$count, weights=gsa$par/sum(gsa$par))
# write.table(gsa$trace.mat, file='trace.mat_temp100.txt', quote=F, sep='\t')

# set.seed(1234)
# date()
# gsa = RunGenSA(spectra[,1], cosmic, list(maxit=1000, threshold.stop=0.008243543*1.01, temperature=100, simple.function=T, verbose=F))
# list(error=gsa$value, fn.count=gsa$count, weights=gsa$par/sum(gsa$par))
# date()








#end
