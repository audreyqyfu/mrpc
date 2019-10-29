mpinv <- function(X)
{
# Moore-Penrose inverse
Eps <- 100000*.Machine$double.eps;# 100*.Machine$double.eps;

# singular value decomposition
s <- svd(X);
#s=SVDmiss(X)
#d <- s$svd$d;
d <- s$d;
m <- length(d);
if (!(is.vector(d)))
return(t(s$v%*%(1/d)%*%t(s$u)));

# remove eigenvalues equal zero
d <- d[d > Eps];
notnull <- length(d);
if (notnull == 1)
{
inv <- 1/d;
} else {
inv <- solve(diag(d));
}

# add rows, columns of zeros if needed 
if (notnull != m)
{
inv <- cbind(inv, matrix(0, nrow=notnull, ncol=(m - notnull)));
inv <- rbind(inv, matrix(0, nrow=(m-notnull), ncol=m));
} 

# compute Moore-Penrose
mp <- s$v%*%inv%*%t(s$u);

# set very small values to zero
mp[abs(mp) < Eps] <- 0;
return(mp);
}

#X <- cbind(1, diag(3)); # singular matrix
#y <- 1:3
#Xp <- qr(X);
#b1 <- qr.coef(Xp, y); # contains NA 
#b2 <- mpinv(X)%*%y # least square fit using Moore-Penrose
#X%*%b2 # == y 



