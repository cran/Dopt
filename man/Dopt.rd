\name{Dopt}
\alias{Dopt}
\title{
Finding D-optimal experimental designs.
}
\description{
A function to generate approximately D-optimal experimental
designs using a Fortran routine of Alan Miller and Nam-Ky Nguyen
(AS295).  The algorithm is a combination of the exchange
algorithm of Fedorov and a supplementary algorithm of Cook \&
Nachtsheim (1989) for swapping treatments between blocks.
}
\usage{
Dopt(treatments, blocks, 
	data.treatments = sys.frame(sys.parent()), 
	data.blocks     = sys.frame(sys.parent()), 
	initial.points  = rep(0, npoints), 
	fixed.points    = rep(FALSE, npoints), 
	RNG             = as.integer(runif(3) * 30269), 
	eps             = .Machine$double.eps, 
	runs)
}
\arguments{
\item{treatments}{
A number, factor, matrix or formula specifying how the full
design matrix of candidate points is to be constructed.  (If this
argument is a single number it is interpreted as a number of
simple treatments to be applied in a single factor experiment.)  
}
\item{blocks}{
A number, factor or vector specifying the number of runs
in the experiment and the block structure, if any.  (If this
argument is a single number it specifies the number of runs in a
single block experiment.)
}
\item{data.treatments}{
Optional data frame to supply variables for the treatments
specification.
}
\item{data.blocks}{
Optional data frame to supply variables for the blocks
specification.
}
\item{initial.points}{
A vector of numbers specifying, in order, rows of the treatment
design matrix to be used in the initial trial design.  A zero or
missing value denotes an unspecified run.  Any unspecified runs
of the initial design are chosen at random.
}
\item{fixed.points}{
Logical vector; TRUE if the corresponding initial point is to be
kept in the experiment, FALSE if it may be replaced by the algorithm.
}
\item{RNG}{
Vector of three positive integers less than 30269 used to
initialize the Wichmann Hill random number generator incorporated
into the routine.
}
\item{eps}{
Tolerance for zero singular values to detect rank deficiency.
}
\item{runs}{
Alternative name for the \code{blocks} argument.
}}
\value{
A data frame specifying the chosen experimental design.  This
will have
a factor specifying the block structure, (or variables from the
blocks data frame needed to reconstruct it if a formula for
blocks is given).

Variables (or design submatrix) from the treatments specification
specifying the treatments to be used on each run. 

A numeric vector labelled "Point.No" giving the row number of the
treatment design matrix used for each run of the experiment.

If the experiment included fixed runs, a logical vector labelled "Fixed"
specifying which runs were fixed in advance at their initial points.

The data frame has four additional attributes, "log(det)" giving
the log of the final determinant as returned by the algorithm,
which may be useful in comparing two runs for the same design,
"RNG", the three starting integers used for the local random
number generators for this run of the algorithm, "CRE" a
vector of \code{canonical relative efficiencies} and "call", 
a representation of the call used to invoke the algorithm. 
}
\section{Side Effects}{
None.
}
\details{
The program of Miller and Nguyen takes a treatments design
matrix, T, of "candidate points" and forms another matrix, \eqn{X_0}{X0},
from it by selecting rows, possibly with repitition.  If the
experiment has k blocks, and X is the matrix got from \eqn{X_0}{X0} by
subtracting block means from each column, then \eqn{X_0}{X0} is chosen to
make \eqn{det(X'X)} as large as possible.  If this is maximal the
design is D-optimal.  The program uses an exchange algorithm due
to Fedorov, with modifications for computational efficiency, and
exact optimality is not guaranteed so several runs of the
procedure may be useful.

Only simple block structures are allowed.  Zero blocks implies no
intercept term and the matrices X0 and X above are the same.  One
block implies an intercept term, so the columns of X have mean
zero.

The present routine is an R interface to the Miller Nguyen
program.

The \code{blocks} argument has two purposes.  Through its length (or
explicitly, if it is a single number), it prescribes the number
of runs in the experiment.  If it is a vector, it is coerced to
factor and used to define a block structure for the experiment.
If it is a formula it is not interpreted in the ordinary way: any
vectors present are coerced to factors and all terms are then
used to define a "subclasses" block structure.  Thus ~a + b, ~a*b
and ~a/b all define the same block structure in this case.  Note
that non-simple block structures such as for latin square and
row-column designs, or designs for known plot covariates are
not accommodated.  If there are no blocks (see NOTE below) the
blocks argument is only used to define the number of runs.

The treatments argument, together with its data.treatments data
frame, if any, is used to construct the design matrix T of
candidate points.  It either an R formula or effectively coerced
to such, and is interpreted in the usual way.  The rows of the
design matrix, T, or \code{points} are numbered 1 to ncand.

NOTE: In Dopt the only way to specify that the design is not to
have an intercept term is to give the treatments specification as
a formula and explicitly to omit the intercept term, as in 
\eqn{~x + x^2 - 1} for a quadratic regression through the origin.

The optional argument \code{initial.points} is a vector usually of
length equal to the number of runs, with integer entries
specifying rows of T to be used at that position in the initial
trial matrix X0.  Any zero or missing values entries specify rows
that are to be initially chosen at random.

The \code{fixed.points} argument is a logical vector specifying whether
or not the corresponding row of X is to be held fixed at its
initially specified value.  Thus the routine may be used to
augment an existing design with additional points.

The CRE attribute of the result gives the non-zero eigenvalues of
X'X relative to (T'T * npoints/ncand).  Under some circumstances
for an ideal design choice they should all be unity.  In most
cases they will be above and below unity and their harmonic mean
gives one global measure of how well the design succeeds in
capturing potential information, with connexions to A-efficiency. 

NOTE: For some reason, does not works with exactly two blocks, so
this is for present disallowed (at least under windows). Also, in some cases does not
work with exactly 3 blocks, and this may depend on the random seed used.
}

\references{
Cook, R. D. and Nachtsheim, C. J. (1989) Computer-aided blocking
of factorial and response-surface designs.  Technometrics, vol.
31, 339-346.

Miller, A. J. and Nguyen, N-K. (1994). Algorithm AS 295: A
Fedorov exchange algorithm for D-optimal design, Applied
Statistics, vol. 43, 669-678.  (Held in statlib as apstat/295).

Nguyen, N-K. and Miller, A. J. (1992).  A review of some exchange
algorithms for constructing discrete D-optimal designs.
Comput. Statist. \& Data Anal., vol.14, 489-498.
}

\seealso{
S+DOX
}

\examples{
\dontrun{1. Approximate optimal design of a cubic regression with 12 points.}
 x <- 0:100 #
 # Note that the square and cube terms must be protected with I(), in contrast to
 # S-Plus, where this is'nt necessary.
 dsn <- Dopt(~ x + I(x^2) + I(x^3), runs = 12) 
 sort(dsn$x)
 \dontrun{
  1   1   1  28  29  29  73  73  74 101 101 101 
  0   0   0  27  28  28  72  72  73 100 100 100
  } 
 \dontrun{2. A block design for 7 treatments in 7 blocks of size 3:}
 tr <- factor(1:7)
 bl <- c(rep(1:7, rep(3,7)))
 dsn <- Dopt(~tr, bl)
 dsn
\dontrun{
   Blocks tr2 tr3 tr4 tr5 tr6 tr7 Point.No
1       1   0   0   0   1   0   0        5
2       1   1   0   0   0   0   0        2
3       1   0   0   0   0   1   0        6
4       2   0   0   1   0   0   0        4
5       2   0   0   0   1   0   0        5
6       2   0   1   0   0   0   0        3
7       3   0   0   0   0   0   0        1
8       3   0   1   0   0   0   0        3
9       3   1   0   0   0   0   0        2
10      4   0   0   0   0   1   0        6
11      4   0   0   0   0   0   1        7
12      4   0   1   0   0   0   0        3
13      5   0   0   1   0   0   0        4
14      5   1   0   0   0   0   0        2
15      5   0   0   0   0   0   1        7
16      6   0   0   1   0   0   0        4
17      6   0   0   0   0   1   0        6
18      6   0   0   0   0   0   0        1
19      7   0   0   0   0   0   0        1
20      7   0   0   0   1   0   0        5
21      7   0   0   0   0   0   1        7
}
 dsn$tr  <- dsn$tr2
 for(i in seq(along=dsn$tr)) {
    dsn$tr[i] <- if(dsn$tr2[i]==1) 2 else {
                  if(dsn$tr3[i]==1) 3 else {
                   if(dsn$tr4[i]==1) 4 else {
                    if(dsn$tr5[i]==1) 5 else {
                     if(dsn$tr6[i]==1) 6 else {
                      if(dsn$tr7[i]==1) 7 else 1}}}}}}
 # look at the concurrency matrix:
 crossprod(table(dsn$Blocks, dsn$tr))
\dontrun{
  1 2 3 4 5 6 7
1 3 1 1 1 1 1 1
2 1 3 1 1 1 1 1
3 1 1 3 1 1 1 1
4 1 1 1 3 1 1 1
5 1 1 1 1 3 1 1
6 1 1 1 1 1 3 1
7 1 1 1 1 1 1 3
}
 # in this case the design is a BIB and hence optimal.
\dontrun{3. An approximate main effect plan for a 2^11 experiment in 12
   runs.  (In this case a Plackett Burman design could be used
   for comparison).  We also ensure that the first run has all
   factors at their lowest level.
}
 l <- c("-","+")
 fac <- expand.grid(l,l,l,l,l,l,l,l,l,l,l)
 names(fac) <- LETTERS[-6][1:11]      # omit F as factor name
 # dsn <- Dopt(~., 12, fac, initial = 1, fixed = T) # ~. doesn't work in R --- by now.
 dsn <- Dopt(~ A+B+C+D+E+G+H+I+J+K+L,12, fac, initial=1, fixed=T) 
 # sometimes gives: Error: error  8  in dsvdc ???     Just try again.
 dsn   
  # in large cases several tries of the algorithm might be useful.
  # 3 runs of above example gives lndet in range 1.2 -- 1.8 
 # Some examples where there is problems, at least under windows:
 \dontrun{
     dsn <- Dopt(~ A+B+C+D+E+G+H+I+J+K+L, c(rep(1,7), rep(2,7)), fac,
                   initial=1, fixed=T)
     dsn <- Dopt(~ A+B+C+D+E+G+H+I+J+K+L, c(rep(1,8), rep(2,8)), fac,
                   initial=1, fixed=T) 
 }
}
\keyword{design}
% Converted by Sd2Rd version 0.3-2.
