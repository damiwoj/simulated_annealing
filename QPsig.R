library(quadprog)

QPsig <- function(tumour.ref, signatures.ref) {
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
