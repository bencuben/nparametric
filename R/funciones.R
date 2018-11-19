# Function 1 --------------------------------------------------------------

#' Cramer V
#' 
#' description here!!!!
#' @param x first element
#' @param y second element
#' @param digits number of digits you want in the result
#' @param ... Another things that idk
#' @return a numeric value indicating the rao crammer v
#' @export
#' @references librito de no parametrica
CramerV<-function (x, y = NULL, digits = 4, ...) 
{
  CV = NULL
  if (is.factor(x)) {
    x = as.vector(x)
  }
  if (is.factor(y)) {
    y = as.vector(y)
  }
  if (is.vector(x) & is.vector(y)) {
    N = length(x)
    Chi.sq = suppressWarnings(chisq.test(x, y, correct = FALSE, 
                                         ...)$statistic)
    Phi = Chi.sq/N
    R = length(unique(x))
    C = length(unique(y))
    CV = sqrt(Phi/min(R - 1, C - 1))
  }
  if (is.matrix(x)) {
    x = as.table(x)
  }
  if (is.table(x)) {
    N = sum(x)
    Chi.sq = suppressWarnings(chisq.test(x, correct = FALSE, 
                                         ...)$statistic)
    Phi = Chi.sq/N
    R = nrow(x)
    C = ncol(x)
    CV = sqrt(Phi/min(R - 1, C - 1))
  }
  
  CV = signif(as.numeric(CV), digits = digits)
  names(CV) = "Cramer V"
  return(CV)
}


# Function 2 --------------------------------------------------------------

#' Cramer contingence coeficent
#' 
#' description here!!!!
#' @param x first element 
#' @param y second element
#' @param digits number of digits you want in the result
#' @param ... Another things that idk
#' @return a numeric value indicating the contingece coefficient of Crammer 
#' @export
#' @references librito de no parametrica
CoefCont<-function (x, y = NULL, digits = 4, ...) 
{
  Cont = NULL
  if (is.factor(x)) {
    x = as.vector(x)
  }
  if (is.factor(y)) {
    y = as.vector(y)
  }
  if (is.vector(x) & is.vector(y)) {
    N = length(x)
    Chi.sq = suppressWarnings(chisq.test(x, y, correct = FALSE, 
                                         ...)$statistic)
    R = length(unique(x))
    C = length(unique(y))
    Cont = sqrt(Chi.sq/(N+Chi.sq))
  }
  if (is.matrix(x)) {
    x = as.table(x)
  }
  if (is.table(x)) {
    N = sum(x)
    Chi.sq = suppressWarnings(chisq.test(x, correct = FALSE, 
                                         ...)$statistic)
    R = nrow(x)
    C = ncol(x)
    Cont = sqrt(Chi.sq/(N+Chi.sq))
  }
  
  Cont = signif(as.numeric(Cont), digits = digits)
  names(Cont) = "Coef. Cont. de Cramer"
  return(Cont)
}

# Function 3 --------------------------------------------------------------

#' Phi value
#' 
#' description here!!!!
#' @param x first element 
#' @param y second element
#' @param digits number of digits you want in the result
#' @param ... Another things that idk
#' @return a numeric value indicating the Phi value
#' @export
#' @references librito de no parametrica
Phi<-function (x, y = NULL, digits = 4, ...) 
{
  Phi = NULL
  if (is.factor(x)) {
    x = as.vector(x)
  }
  if (is.factor(y)) {
    y = as.vector(y)
  }
  if (is.vector(x) & is.vector(y)) {
    N = length(x)
    Chi.sq = suppressWarnings(chisq.test(x, y, correct = FALSE, 
                                         ...)$statistic)
    Phi = sqrt(Chi.sq/N)
  }
  if (is.matrix(x)) {
    x = as.table(x)
  }
  if (is.table(x)) {
    N = sum(x)
    Chi.sq = suppressWarnings(chisq.test(x, correct = FALSE, 
                                         ...)$statistic)
    Phi = sqrt(Chi.sq/N)
  }
  
  Phi = signif(as.numeric(Phi), digits = digits)
  names(Phi) = "Phi"
  return(Phi)
}

# Function 4 --------------------------------------------------------------

#' Cochranq Hyphotesys test for matrix
#' 
#' description here!!!!
#' @param mat a matrix object
#' @return the Cochranq hyphotesys test
#' @export
#' @references librito de no parametrica
cochranq.test <- function(mat)
{ k <- ncol(mat)
C <- sum(colSums(mat) ^ 2)
R <- sum(rowSums(mat) ^ 2)
T <- sum(rowSums(mat))
num <- (k - 1) * ((k * C) - (T ^ 2))
den <- (k * T) - R
Q <- num / den
df <- k - 1
names(df) <- "df"
names(Q) <- "Cochran's Q"
p.val <- pchisq(Q, df, lower = FALSE)
QVAL <- list(statistic = Q, parameter = df, p.value = p.val,
             method = "Cochran's Q Test for Dependent Samples",
             data.name = deparse(substitute(mat)))
class(QVAL) <- "htest"
return(QVAL)}

# Function 5 --------------------------------------------------------------

#' Monotone no parametric regresion 
#' 
#' description here!!!!
#' @param data a vector, array or data frame
#' @return Ni idea que retorna
#' @export
#' @references librito de no parametrica
monregnp<-function(data) {
  
  x=data[,1]
  y=data[,2]
  n=length(x)
  rxi=rank(x)
  ryi=rank(y)
  ry_ajus=lm(ryi~rxi)
  coef=ry_ajus$coefficients
  names(coef)<-NULL
  ry0=ry_ajus$fitted.values
  names(ry0)<-NULL
  
  orden=order(y)
  yor=y[orden]
  xor=x[orden]
  ry_ajust=ry0[orden]
  ryiord=ryi[orden]
  maxi=max(ryi)
  maxy=max(y)
  mini=min(ryi)
  miny=min(y)
  
  y_ajust=NULL
  
  for (i in 1:n) {
    if ( ry_ajust[i]>maxi) ( y_ajust[i]=maxy) else
      if( ry_ajust[i]<mini) ( y_ajust[i]=miny)
    for (j in 1:n-1)  {
      if (ry_ajust[i]>ryiord[j] && ry_ajust[i]<ryiord[j+1])
        (y_ajust[i]= yor[j]+(yor[j+1]-yor[j])*
           (ry_ajust[i]-ryiord[j])/(ryiord[j+1]-ryiord[j]))
      
    }
  }
  cbind(xor,y_ajust)
  
  
  rx_ajust=(1/coef[2])*(ryi-coef[1])
  maxxa=max(rxi)
  maxx=max(x)
  minxa=min(rxi)
  minx=min(x)
  
  x_ajust=NULL
  for (i in 1:n) {
    if ( rx_ajust[i]>maxxa || rx_ajust[i]<minxa) ( x_ajust[i]=999999)
    for (j in 1:n-1)  {
      if (rx_ajust[i]>=rxi[j] && rx_ajust[i]<=rxi[j+1])
        (x_ajust[i]= x[j]+(x[j+1]-x[j])*
           (rx_ajust[i]-rxi[j])/(rxi[j+1]-rxi[j]))
      
    }
  }
  
  y=unique(y)
  
  x_ajust=unique(x_ajust)
  
  y=y[x_ajust!=999999]
  
  x_ajust=x_ajust[x_ajust!=999999]
  
  x1=c(xor,x_ajust)
  
  y1=c(y_ajust,y)
  
  ord=order(x1)
  
  x1=x1[ord]
  
  y1=y1[ord]
  
  return(cbind(x1,y1))
}
