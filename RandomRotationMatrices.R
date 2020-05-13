################################################################################
################################################################################
#
# PURPOSE OF THIS SCRIPT
# This script produces random rotation matrices, "Q.new", via QR decomposition.
# These "Q.new" matrices are uniformely distributed with the Haar measure in O(n).
# They will have determinant = -1 with probability 0.5. If rejection sampling is
# applied to exclude "Q.new" matrices with determinant -1 and accept only those
# with determinant 1, then the selected "Q.new" matrices will be SO(n) with a
# uniform distribution relative to the Haar measure. See:
# https://kaba.hilvi.org/pastel-1.5.0/pastel/math/matrix/random_matrix/random_orthogonal_matrix.htm
# https://kaba.hilvi.org/pastel-1.5.0/pastel/math/matrix/random_matrix/random_rotation_matrix.htm
# How to Generate Random Matrices from the Classical Compact Groups, Francesco Mezzadri,
# Notices of the AMS, Volume 54, Number 5, 2007.
# https://math.stackexchange.com/questions/494225/what-is-haar-measure/494247
#
# DATA FILES REQUIRED TO RUN THIS SCRIPT
# This script may run on its own, without any data files.
#
################################################################################
################################################################################

#define dimensionality
d <- 15

#generate random matrix
X <- matrix(rnorm(d^2), ncol=d, nrow=d)

#QR decomposition
Q <- qr.Q(qr(X))
R <- qr.R(qr(X))

#make the QR decomposition unique by fixing
#the signs of diagonal of R
sgn <- sign(diag(R))
R.new <- diag(sgn) %*% R
Q.new <- Q %*% diag(sgn)

#check the results
all.equal(Q.new%*%R.new, X)
summary(Q.new%*%R.new - X)
identical(Q.new%*%R.new, X)

#determinant to be used in
#rejection sampling
det(Q.new)


