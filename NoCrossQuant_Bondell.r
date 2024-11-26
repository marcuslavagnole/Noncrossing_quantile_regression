rq.no.cross = function(y, x, taus)
{
    require(quantreg)
     
    if (length(taus)<2)
	  stop("At least 2 quantile levels should be specified. If only using a single quantile, you should use rq.")

    taus = sort(taus)
    if (max(taus)>=1)
	  stop("All quantile levels should be between 0 and 1, not including the boundary.")
    if (min(taus)<=0)
	  stop("All quantile levels should be between 0 and 1, not including the boundary.")
    
    options(warn=-1)
    n = nrow(x)
    p = ncol(x)
    m = length(taus)
    pp = p+1
    Y = rep(y, m)
    xtemp = x

    x = matrix(0,nrow=n,ncol=p)
    shifts = apply(xtemp,2,min)
    scalings = rep(0,p)
    for (i in 1:p)
    {
     	x[,i] = xtemp[,i] - shifts[i]
	scalings[i] = max(x[,i])
	x[,i] = x[,i]/scalings[i]
    }

    x = cbind(rep(1,n),x)

    D = diag(m)
    D[lower.tri(D)] = 1 
    X = D %x% x
    X2 = cbind(X, -X)
    sX = as.matrix.csr(X2)
    K = m*pp
    R1 = (diag(m) %x% rep(1,1) %x% t(c(1, rep(0, p))))[-(1:1),]
    R2 = rep(1,1)
    for (j in 1:p)
    {
    	R2 = cbind(R2, rep(1,1) %x% (diag(1) %x% rep(1,1)))
    }

    R2 = (diag(m) %x% R2)[-(1:1),]
    
    sR = as.matrix.csr(rbind(diag(2*K), cbind(R1, -R2)))
    r2 = rep(0, 2*K + (m-1)*(1))
       
    tau.v = rep(taus, each=n)
    rhs = t(X2)%*%(1-tau.v)
    coeff1 =  myrq.fit.sfnc(sX, Y, sR, r2, tau=tau.v, rhs=rhs,tmpmax=100000)$coef
    coeff = coeff1[1:K]-coeff1[-(1:K)]
    gamma.m = matrix(coeff, ncol=m)
    D = diag(m)
    D[upper.tri(D)]=1
    bhat.temp = gamma.m %*% D
    cov.bhat = array(0,c(pp,pp,m))
    se.bhat = matrix(0,ncol=m,nrow=pp)
    transform.mat = rbind(c(1,-shifts/scalings),cbind(rep(0,p),diag(as.vector(1/scalings), nrow = length(as.vector(1/scalings)), ncol = length(as.vector(1/scalings)))))
    bhat = transform.mat%*%bhat.temp
    for (j in 1:m)
    {
	cov.bhat[,,j] = se.constr(x,y,bhat.temp[,j],taus[j])
   	cov.bhat[,,j] = transform.mat%*%cov.bhat[,,j]%*%t(transform.mat)
	se.bhat[,j] = sqrt(diag(cov.bhat[,,j])) 
    }
	
    vars<-c("intercept",paste("x",1:p,sep=""))
    rownames(bhat)<-rownames(se.bhat)<-vars
    colnames(bhat)<-colnames(se.bhat)<-taus
    dimnames(cov.bhat)<-list(vars,vars,taus) 
	
    constr.qr.fit = NULL
    constr.qr.fit$bhat = bhat
    constr.qr.fit$se.bhat = se.bhat
    constr.qr.fit$cov.bhat = cov.bhat
    constr.qr.fit$taus = taus
    
    return(constr.qr.fit)
}


myrq.fit.sfnc <- function (x, y, R, r, tau, rhs, nsubmax, tmpmax, nnzlmax, cachsz = 64, 
    small = 1e-08, maxiter = 100, warn.mesg = TRUE) 
{
    require(quantreg)
    y <- -y
    r <- -r
    n1 <- length(y)
    m <- x@dimension[2]
    if (n1 != x@dimension[1]) 
        stop("The design matrix A1' and response vector y are not compatible")
    n2 <- length(r)
    if (n2 != R@dimension[1]) 
        stop("The constraint matrix A2' and the constraint right-hand-side are not compatible")
    maxn1n2 <- max(n1, n2)
    u <- rep(1, length = n1)
    if(length(tau)==1)     x1 <- rep(1 - tau, length = n1)
    else x1=1-tau
    x2 <- rep(1, length = n2)
    wwm <- vector("numeric", 6 * m)
    wwm[1:m] <- rhs
    nnzx <- x@ia[x@dimension[1] + 1] - 1
    nnzR <- R@ia[R@dimension[1] + 1] - 1
    nnzdmax <- max(nnzx, nnzR)
    iwmax <- 7 * m + 3
    ao1 <- t(x)
    ao2 <- t(R)
    e <- ao1 %*% x
    g <- ao2 %*% R
    h <- e + g
    nnzemax <- e@ia[e@dimension[1] + 1] - 1
    nnzgmax <- g@ia[g@dimension[1] + 1] - 1
    nnzhmax <- h@ia[h@dimension[1] + 1] - 1
    if (missing(nnzlmax)) 
        nnzlmax <- 4 * nnzdmax
    if (missing(nsubmax)) 
        nsubmax <- nnzhmax
    if (missing(tmpmax)) 
        tmpmax <- 6 * m
    s <- u - x1
    chol.o <- chol(e, tmpmax = tmpmax, nsubmax = nsubmax, nnzlmax = nnzlmax)
    b <- backsolve(chol.o, ao1 %*% y)
    r1 <- y - x %*% b
    z1 <- ifelse(abs(r1) < small, (r1 * (r1 > 0) + small), r1 * 
        (r1 > 0))
    w <- z1 - r1
    z2 <- rep(1, n2)
    wwn1 <- matrix(0, n1, 10)
    wwn1[, 1] <- z1
    wwn1[, 2] <- w
    wwn2 <- matrix(0, n2, 7)
    wwn2[, 2] <- z2
    srqfnc.o <- .Fortran("srqfnc", n1 = as.integer(n1), m = as.integer(m), 
        nnzx = as.integer(nnzx), x = as.double(x@ra), jx = as.integer(x@ja), 
        ix = as.integer(x@ia), ao1 = as.double(ao1@ra), jao1 = as.integer(ao1@ja), 
        iao1 = as.integer(ao1@ia), n2 = as.integer(n2), nnzR = as.integer(nnzR), 
        R = as.double(R@ra), jR = as.integer(R@ja), iR = as.integer(R@ia), 
        ao2 = as.double(ao2@ra), jao2 = as.integer(ao2@ja), iao2 = as.integer(ao2@ia), 
        nnzdmax = as.integer(nnzdmax), d = double(nnzdmax), jd = integer(nnzdmax), 
        id = integer(m + 1), dsub = double(nnzhmax + 1), jdsub = integer(nnzhmax + 
            1), nnzemax = as.integer(nnzemax), e = as.double(e@ra), 
        je = as.integer(e@ja), ie = as.integer(e@ia), nnzgmax = as.integer(nnzgmax), 
        g = double(nnzgmax), jg = integer(nnzgmax), ig = integer(m + 
            1), nnzhmax = as.integer(nnzhmax), h = double(nnzhmax), 
        jh = integer(nnzhmax), ih = integer(m + 1), nsubmax = as.integer(nsubmax), 
        lindx = integer(nsubmax), xlindx = integer(m + 1), nnzlmax = as.integer(nnzlmax), 
        lnz = double(nnzlmax), xlnz = integer(m + 1), iw = integer(m * 
            5), iwmax = as.integer(iwmax), iwork = integer(iwmax), 
        xsuper = integer(m + 1), tmpmax = as.integer(tmpmax), 
        tmpvec = double(tmpmax), maxn1n2 = as.integer(maxn1n2), 
        ww1 = double(maxn1n2), wwm = as.double(wwm), wwn1 = as.double(wwn1), 
        wwn2 = as.double(wwn2), cachsz = as.integer(cachsz), 
        level = as.integer(8), x1 = as.double(x1), x2 = as.double(x2), 
        s = as.double(s), u = as.double(u), c1 = as.double(y), 
        c2 = as.double(r), sol = as.double(b), small = as.double(small), 
        ierr = integer(1), maxiter = as.integer(maxiter), time = double(7), 
        PACKAGE = "quantreg")[c("sol", "ierr", "maxiter", "time")]
    ierr <- srqfnc.o$ierr
    if (ierr == 13) 
        stop("Increase nnzh.factor")
    if (!(ierr == 0) && warn.mesg) 
        warning(sfnMessage(ierr))
    list(coef = -srqfnc.o$sol, ierr = ierr, it = srqfnc.o$maxiter, 
        time = sum(srqfnc.o$time))
}

se.constr = function(x,y,coef,tau)
{
        require(quantreg)
        n = nrow(x)
        p = ncol(x)
        h = bandwidth.rq(tau, n)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
        uhat = c(y - x %*% coef)
        h = (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)), 
            (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
        f = dnorm(uhat/h)/h
        fxxinv = diag(p)
        fxxinv = backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p, drop = FALSE], 
            fxxinv)
        fxxinv = fxxinv %*% t(fxxinv)
        cov = tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
        cov
}
