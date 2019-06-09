# Function 1 --------------------------------------------------------------

#' Cramer Contingence Coefficient 
#' 
#' @description Measure a degree of dependence between two variables, this number should be between 0 and 1.
#' @param x A numeric vector or matrix. x, and y can also both be factors.
#' @param y A numeric vector; ignored if x is a matrix. If x is a factor, y should be a factor of the same length.
#' @param digits Integer indicating the number of decimal places (round) or significant digits (signif) to be used. Negative values are allowed (see ‘Details’).
#' @param ... Additional arguments provided by chisq.test() from stats package.
#' @details Rounding to a negative number of digits means rounding to a power of ten, so for example round(x, digits = -2) rounds to the nearest hundred.
#' @return A numeric value indicating the Cramer Contingence Coefficient. 
#' @export
#' @examples x<- matrix(c(28,5,0,7),ncol=2)
#' CramerV(x)
#' @references Conover, W. J. (1999) "Practical Nonparametric Statistics", 3rd Edition, Jhon Wiley & Sons.Inc
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

#' Contingence Coefficent of Pearson
#' 
#' @description  Measure a degree of dependence between two variables, this number shoud be between 0 and \eqn{\sqrt(q-1/q)}.
#' @param x A numeric vector or matrix. x and y can also both be factors. 
#' @param y A numeric vector; ignored if x is a matrix. If x is a factor, y should be a factor of the same length.
#' @param digits Integer indicating the number of decimal places (round) or significant digits (signif) to be used. Negative values are allowed (see ‘Details’).
#' @param ... Additional arguments provided by chisq.test() from stats package.
#' @return A numeric value indicating the contingece coefficient of Pearson. 
#' @details Rounding to a negative number of digits means rounding to a power of ten, so for example round(x, digits = -2) rounds to the nearest hundred.
#' @export
#' @examples x<- matrix(c(28,5,0,7),ncol=2)
#' CoefCont(x)
#' @references Conover, W. J. (1999) "Practical Nonparametric Statistics", 3rd Edition, Jhon Wiley & Sons.Inc
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
  names(Cont) = "Coef. Cont. de Cramer" #no sería de pearson?
  return(Cont)
}

# Function 3 --------------------------------------------------------------

#' Phi Coefficient
#' 
#' @description  Measure a degree of dependence between two variables, this number shoud be between 0 and \eqn{\sqrt(q-1)}.
#' @param x A numeric vector or matrix. x and y can also both be factors. 
#' @param y A numeric vector; ignored if x is a matrix. If x is a factor, y should be a factor of the same length.
#' @param digits Integer indicating the number of decimal places (round) or significant digits (signif) to be used. Negative values are allowed (see ‘Details’).
#' @param ... Additional arguments provided by chisq.test() from stats package.
#' @return A numeric value indicating the contingece coefficient of Pearson.
#' @details Rounding to a negative number of digits means rounding to a power of ten, so for example round(x, digits = -2) rounds to the nearest hundred.
#' @export
#' @examples x<- matrix(c(28,5,0,7),ncol=2)
#' Phi(x)
#' @references Conover, W. J. (1999) "Practical Nonparametric Statistics", 3rd Edition, Jhon Wiley & Sons.Inc
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

#' Cochran's Q Test for Dependent Samples
#' 
#' @description A mean difference test for depending samples proposed by Marascuilo & McSweeny (1967).
#' @param x A matrix or data.frame 
#' @return A list with class "htest" containing the following components:
#'
#' statistic:	the value of the Cochran's Q statistic.
#' 
#' df: Degree of freedom of the Chi squared distribution corresponding to the test. 
#' 
#' p.value: an approximate p-value for the test.
#' @export
#' @examples x = as.table(matrix(c(1,1,0,1,0,1,1,1,0,0,1,1,1,1,1,1,0,
#' 1,1,1,0,1,1,1,1,1,0,0,0,1,1,0,1,0,1,1),ncol=3))
#' cochranq.test(x)
#' @references Conover, W. J. (1999) "Practical Nonparametric Statistics", 3rd Edition, Jhon Wiley & Sons.Inc
#' 
#' 
cochranq.test <- function(x)
{ k <- ncol(x)
C <- sum(colSums(x) ^ 2)
R <- sum(rowSums(x) ^ 2)
T <- sum(rowSums(x))
num <- (k - 1) * ((k * C) - (T ^ 2))
den <- (k * T) - R
Q <- num / den
df <- k - 1
names(df) <- "df"
names(Q) <- "Cochran's Q"
p.val <- pchisq(Q, df, lower = FALSE)
QVAL <- list(statistic = Q, parameter = df, p.value = p.val,
             method = "Cochran's Q Test for Dependent Samples",
             data.name = deparse(substitute(x)))
class(QVAL) <- "htest"
return(QVAL)}

# Function 5 --------------------------------------------------------------

#' Monotone no parametric regression 
#' 
#' @description This function is used to fit Monotone nonparametric regression for a bi-variated vector of data
#' @param datos a vector, array or data frame
#' @return Return a array, containing the estimated points for the adjusted curve
#' @export
#' @examples x=c(0,0.5,1.0,1.8,2.2,2.7,4.0,4.0,4.9,5.6,6.0,6.5,7.3,8.0,8.8,9.3,9.8)
#' y=c(30,30,30,28,24,19,17,9,12,12,6,8,4,5,6,4,6)
#' d=cbind(x,y)
#' est=monregnp(d)
#' plot(est,type="l",xlab="libras de azúcar", ylab="Días de fermentación")
#' points(x,y)
#' title(main="Curva de regresión no paramétrica")
#' @references Conover, W. J. (1999) "Practical Nonparametric Statistics", 3rd Edition, Jhon Wiley & Sons.Inc
monregnp<-function(datos) {
  
  x=datos[,1]
  y=datos[,2]
  n=length(x)
  rxi=rank(x)
  ryi=rank(y)
  ry_ajus=lm(ryi~rxi)
  coef=ry_ajus$coefficients
  names(coef)<-NULL
  ry0=ry_ajus$fitted.values
  names(ry0)<-NULL
  rx_a=(1/coef[2])*(ryi-coef[1])
  names(rx_a)<-NULL
  
  
  orden=order(y)
  yor=y[orden]
  xor=x[orden]
  ry_ajust=ry0[orden]
  ryiord=ryi[orden]
  maxi=max(ryi)
  maxy=max(y)
  mini=min(ryi)
  miny=min(y)
  
  y_ajust=rep(0,n)
  
  for (i in 1:n) {
    if ( ry_ajust[i]>maxi) ( y_ajust[i]=maxy) else
      if( ry_ajust[i]<mini) ( y_ajust[i]=miny) else
        for (j in 1:n-1)  {
          if (ry_ajust[i]>=ryiord[j] && ry_ajust[i]<=ryiord[j+1])
            (y_ajust[i]= yor[j]+(yor[j+1]-yor[j])*
               (ry_ajust[i]-ryiord[j])/(ryiord[j+1]-ryiord[j]))
          
        }
  }
  cbind(xor,yor,ry_ajust,ryiord,y_ajust)
  
  
  orden1=order(x)
  xord=x[orden1]
  yord=y[orden1]
  rx_ajust=rx_a[orden1]
  rxiord=rxi[orden1]
  maxxa=max(rxi)
  maxx=max(x)
  minxa=min(rxi)
  minx=min(x)
  
  
  
  x_ajust=rep(0,n)
  for (i in 1:n) {
    if (rx_ajust[i]>maxxa) (x_ajust[i]=999999) else
      if (rx_ajust[i]<minxa) (x_ajust[i]=999999) else
        for (j in 1:n-1)  {
          if (rx_ajust[i]>=rxiord[j] && rx_ajust[i]<=rxiord[j+1])
            (x_ajust[i]= xord[j]+(xord[j+1]-xord[j])*(rx_ajust[i]-rxiord[j])/(rxiord[j+1]-rxiord[j]))
        }
  }
  
  cbind(xord,yord,rx_ajust,rxiord,x_ajust)
  
  xor<-unique(xor)
  
  yord=unique(yord)
  
  x_ajust=unique(x_ajust)
  
  yord=yord[x_ajust!=999999]
  
  x_ajust=x_ajust[x_ajust!=999999]
  
  x1=c(xor,x_ajust)
  
  y1=c(y_ajust,yord)
  
  ord=order(x1)
  
  x1=x1[ord]
  
  y1=y1[ord]
  
  return(cbind(x1,y1))
}

# Function 6 --------------------------------------------------------------

#' Variance comparison between two independent samples
#' 
#' @description Compare de variance of two independent samples, depending of the alternative hypothesis
#' @param x A vector or array for the first sample
#' @param y A vector or array for the second sample
#' @param test A character value containing one of the following option: "dos-colas", "inferior", "superior"
#' @return Return a list containing the statistic's value and the P-value for the test
#' @export
#' @examples x<-c(10.8,11.1,10.4,10.1,11.3)
#' y<-c(10.8,10.5,11,10.9,10.8,10.7,10.8)
#' compvar(x,y,test="superior")
#' @references Conover, W. J. (1999) "Practical Nonparametric Statistics", 3rd Edition, Jhon Wiley & Sons.Inc
compvar<-function(x,y,test="dos-colas"){
  
  m1<-mean(x)
  
  m2<-mean(y)
  
  n<-length(x)
  
  m<-length(y)
  
  N<-n+m
  
  k<-n+1
  
  u<-abs(x-m1)
  
  v<-abs(y-m2)
  
  w<-c(u,v)
  
  R<-rank(w)
  
  Ru<-R[1:n]
  
  Rv<-R[k:N]
  
  T<-sum(Ru^2)
  
  R2<-mean(R^2)
  
  R4<-mean(R^4)
  
  T1<-(T-n*R2)/sqrt(n*m/(N-1)*(R4-R2^2))
  
  Vp <- switch(test,
               
               "dos-colas" = 2*(1-pnorm(abs(T1))), "inferior" = pnorm(T1),
               
               "superior" = 1-pnorm(T1) )
  
  
  
  val <- list(Z = T1, Valp = Vp)
  
  return(val)
  
}


