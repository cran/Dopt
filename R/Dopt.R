#  Uses function "is.R()"  to mark changes, as much as possible. That is, (ALMOST) all the  S
#    code is there as well, facilitating comparisons. 

Dopt <- function(treatments, 
		 blocks, 
		 data.treatments = if (is.R() ) sys.frame(sys.parent()) else 
                                     sys.parent(),                                              #
		 data.blocks     = if (is.R() ) sys.frame(sys.parent()) else
                                     sys.parent(),                                              #
		 initial.points  = rep(0, npoints), 
		 fixed.points    = rep(FALSE, npoints), 
		 RNG             = round(runif(3) * 30269),                       
		 eps             = if ( is.R() ) .Machine$double.eps else .Machine$single.eps,  #
		 runs)
{
    if ( is.R() ) eps <- eval(eps)  # Seems to be necessary?

    if(missing(blocks)) {
        if(missing(runs))
            stop("Either blocks or runs must be specified")
        else blocks <- runs
    }

    if(missing(treatments))
        stop("treatments must be specified")

    call <- match.call()


    if(inherits(blocks, "formula")) {
        vars <- attr(terms(blocks, data = data.blocks), "variables")

        block.data.frame <- 
            if(is.R() ){
               stop("Formulas for blocks not implemented in R yet.")
            } else {
            as.data.frame(lapply(as.list(vars), 
	    function(x, dat) 
	    as.factor(eval(as.expression(x), local = dat)), 
            data.blocks))
            }
        names(block.data.frame) <- as.character(vars)
        vars <- lapply(unclass(block.data.frame), 
            function(s) format(as.character(s)))
        bf <- factor(do.call("paste", c(vars, list(sep = ""))))
    } else {  # blocks inn't a formula 
        if(length(blocks) == 1) {
            bf <- rep(1, as.numeric(blocks))
            block.data.frame <- data.frame("(Intercept)" = bf)
        }
        else {
            bf <- as.factor(blocks)
            block.data.frame <- data.frame(Blocks = bf)
        }
    }
    npoints <- length(bf)

    if(!inherits(treatments, "formula")) {
        if(length(treatments) == 1)
            Treatments <- as.factor(1:as.numeric(treatments))
        else {
            Treatments <- eval(treatments, data.treatments)
        }
        treatments <-  ~ Treatments
        data.treatments <- sys.frame()
    }
    X <- model.matrix(treatments, data.treatments)
    intercept <- any(as.logical(match(dimnames(X)[[2]], 
        "(Intercept)", nomatch = 0)))
    if(intercept) {
        X <- scale(X, scale = FALSE)
        bf <- factor(as.numeric(bf))
        nblocks <- length(levels(bf))
        if(nblocks == npoints)
            stop("All blocks have only one plot.")
        blocksizes <- table(bf)
        if(nblocks == 1)
            names(block.data.frame) <- "(Intercept)"
    }
    else {
        nblocks <- 0
        blocksizes <- npoints
        block.data.frame <- data.frame(.Intercept. = bf)
        o <- 1:npoints
    }
    s <- svd(X)
    X <- s$u[, s$d > eps, drop=FALSE]                      # Risks loosing matrixness
    ncand <- dim(X)[1]
    kin <- dim(X)[2]
    k <- kin + nblocks
    nrbar <- (k * (k - 1))/2
    X <- X * sqrt(ncand/npoints)

    # There are problems (at least in windows) with exactly two blocks. --- GPF
    # We disallow this with windows:
    # if(machine()=="Win32" && nblocks==2) stop("Dopt with exactly 2 blocks disallowed, malfunctions.")

    storage.mode(X) <- if ( is.R() ) "double" else "single"      #

    rstart <- missing(initial.points)
    if(!rstart) {
        initial.points <- round(as.numeric(initial.points))
        fixed.points <- as.logical(fixed.points)
        length(initial.points) <- npoints
        length(fixed.points) <- npoints
        initial.points[is.na(initial.points)] <- 0
        fixed.points[is.na(fixed.points)] <- FALSE
        if(any(initial.points < 0 | initial.points > ncand))
            stop("initial points invalidly specified")
    }
    in.nbl <- 0
    m <- sum(where <- (initial.points == 0))
    initial.points[where] <- sample(1:ncand, m, rep = TRUE)
    if(nblocks > 0) {
        o <- order(bf, !fixed.points)
        block.data.frame <- block.data.frame[o,  , drop = FALSE]
        initial.points <- initial.points[o]
        fixed.points <- fixed.points[o]
        in.nbl <- tapply(fixed.points, bf, sum)
    }

#   We translated the fortran in dopt.f to C by f2c, but that keeps fortran calling
#   conventions, AND adds trailing underscores to the names, so we continue calling by 
#   .Fortran below: 

    .Fortran("setran",
        as.integer(RNG))      # without change in R. 

# Strangely enough, the problem with bombing when 2 blocks seem 
# to disappear when call is made through .Fordebug, EVEN when ADD=0!
if (FALSE) {
    fo <- .Fordebug("fedrov", 
                   X        = X,
	           dim1     = as.integer(ncand),
	           ncand    = as.integer(ncand),
	           kin      = as.integer(kin),
	           n        = as.integer(npoints),
	           nblock   = as.integer(nblocks),
	           in.nbl   = as.integer(in.nbl),
	           blocksiz = as.integer(blocksizes),
	           k        = as.integer(k),
	           rstart   = as.logical(rstart),
	           nrbar    = as.integer(nrbar),
	           d        = double(k),
	           rbar     = double(nrbar),
	           picked   = as.integer(initial.points),
	           lndet    = double(1),
	           xx       = double(k),
	           tol      = double(k),
	           zpz      = double(ncand * max(1, nblocks)),
	           xk       = double(k),
	           ifault   = integer(1),
                   ADD=0,  ADD.VALUE=-98989)
} # Then we do this again, through do.call:

    args  <- list(NAME= "fedrov", 
                   X        = X,
	           dim1     = as.integer(ncand),
	           ncand    = as.integer(ncand),
	           kin      = as.integer(kin),
	           n        = as.integer(npoints),
	           nblock   = as.integer(nblocks),
	           in.nbl   = as.integer(in.nbl),
	           blocksiz = as.integer(blocksizes),
	           k        = as.integer(k),
	           rstart   = as.logical(rstart),
	           nrbar    = as.integer(nrbar),
	           d        = double(k),
	           rbar     = double(nrbar),
	           picked   = as.integer(initial.points),
	           lndet    = double(1),
	           xx       = double(k),
	           tol      = double(k),
	           zpz      = double(ncand * max(1, nblocks)),
	           xk       = double(k),
	           ifault   = integer(1),
                   NAOK=FALSE, DUP=TRUE)
     fo <- do.call(".Fortran", args)

    if(fo$ifault < 0)
        stop("No full rank starting design found.")
    if(fo$ifault != 0)
        stop("Esoteric fault in fedorov.  Complain!")

    X <- X[fo$picked,  ]
    if(intercept)
        if(nblocks > 1) {
# This call to qr.resid fails when more than one block. 
#               "qr and y must have the same number of rows".
            block.model.matrix <- model.matrix( ~ bf)   # , data = sys.frame() deleted
            block.model.matrix.qr <- qr(block.model.matrix)
            X <- qr.resid(block.model.matrix.qr, X) }
        else {scale(X, scale = FALSE)}
    cre <- sort(svd(X, 0, 0)$d^2)
    vars <- attr(terms(treatments, data = data.treatments), 
        "variables")
# the following gives problems in R:
# we just do without the logical "simple":
    if ( !is.R() ) {
    simple <- all(sapply(as.character(vars), function(x, dat)
	    { exists(x, where = dat) || exists(x) }, 
            data.treatments))
    if(simple) {
        treatment.data.frame <- as.data.frame(lapply(
            as.list(vars), function(x, dat)
        	eval(as.expression(x), local = dat), 
            	data.treatments))[fo$picked,  , drop = FALSE]
        names(treatment.data.frame) <- as.character(vars)
    }
    else {
        treatment.data.frame <- as.data.frame(
            model.matrix(treatments, 
			 data.treatments)[fo$picked,  , drop = FALSE])
    }
    }
    if ( is.R() ) {
       treatment.data.frame <- as.data.frame(
            model.matrix(treatments, 
                         data.treatments)[fo$picked, , drop=FALSE])
    }


    dat <- cbind(block.data.frame, treatment.data.frame, 
        Point.No = fo$picked)
# The following code have problems if there are two Intercepts in dat:
    if ( !is.R() ){
    if (any(i <- match("(Intercept)", names(dat), nomatch = 0)))
        dat <- dat[,  -i]
    if(any(i <- match(".Intercept.", names(dat), nomatch = 0)))
        dat <- dat[,  -i]} 
# We hack that case:
    if ( is.R() ) {
       i <- match("(Intercept)", names(dat), nomatch = 0)
       while (any(i)) { dat <- dat[, -i]
              i <- match("(Intercept)", names(dat), nomatch = 0)}

       i <- match(".Intercept.", names(dat), nomatch = 0)
       while (any(i)) { dat <- dat[, -i]
             i <- match(".Intercept.", names(dat), nomatch = 0)}
    }

    if(any(fixed.points))
        dat$Fixed <- fixed.points
    attr(dat, "log(det)") <- fo$lndet
    attr(dat, "RNG") <- RNG
    attr(dat, "CRE") <- cre[cre > eps]
    attr(dat, "call") <- call
    dat[order(o),  ]
}

