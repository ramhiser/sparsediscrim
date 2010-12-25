\name{LUJID}
\alias{LUJID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LUJID(X, mode = "B", ERR = 1 * 10^-5, RBALANCE = 3, ITER = 200, pseudo = TRUE, shrink.val = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
}
  \item{ERR}{
%%     ~~Describe \code{ERR} here~~
}
  \item{RBALANCE}{
%%     ~~Describe \code{RBALANCE} here~~
}
  \item{ITER}{
%%     ~~Describe \code{ITER} here~~
}
  \item{pseudo}{
%%     ~~Describe \code{pseudo} here~~
}
  \item{shrink.val}{
%%     ~~Describe \code{shrink.val} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(X,mode='B', ERR=1*10^-5, RBALANCE=3, ITER=200, pseudo = TRUE, shrink.val = 0.01){
#LU based Jacbi-like JD; This function minimizes the cost 
#J_{1}(B)=\sum{i=1}^{N} \|BC_{i}B^{T}-diag(BC_{i}B^{T})\|_{F}^{2}
#where \{C_{i}\}_{i=1}^{N} is a set of N, n\times n symmetric matrices 
#and B the joint diagonalizer sought. A related measure that is used
#to measure the error is J_{2}=\sum{i=1}^{N} \|C_{i}-B^{-1}diag(BC_{i}B^{T})B^{-T}\|_{F}^{2}
#
#
#
#Here X is a large matrix of size n\times nN which contains the 
#matrices to be jointly diagonalized such that X=[C1,C2,...,CN], 
#Y contains the jointly diagonalized version of the input 
#matrices, and B is the found diagonalizer.
#
#
#More controlled usage:[Y,B,S,BB]=LUJ1D(X,'mode',ERR or ITER,RBALANCE):
#
#Inputs:
#'mode'='B' or 'E' or 'N':  In the 'B' mode the stopping criteria at each 
#                           step is max(max(abs(LU-I))) which measures 
#                           how much the diagonalizer B has changed
#                           after a sweep. In the 'E' mode 
#                           the stopping criterion is the difference between 
#                           the values of the cost function J2 in two consequtive 
#                           updates.In the 'N' mode the stopping criterion is 
#                           the number of sweeps over L and U phases. 
#
#ERR: In the 'B' mode it specifies the stopping value for the change in B max(max(abs(LU-I))).
#The default value for ERR in this mode and other modes including standard usage 
#is ERR=10^-5. In implementation of the algorithm in order to account 
#for dpendence of accuracy on the dimension n ERR is multiplied 
#by n the size of matrices for JD. In the 'E' mode it ERR specifies the stopping value
#for the relative change of J_{2} in two consequetive sweeps. 
#In the 'B' or 'E' mode or the standard mode
#if the change in B or relative change in J2 does not reach ERR after the default number of 
#iterations (=200) then the program aborts and itreturns the current computed variables.
#
#ITER: Number of iterations in the 'N' mode  
#
##RBALANCE: if given it is the period for row balancing after each sweep.
##################################################################################
#Outputs:
#Y= the diagonalized set of matrices
#B=the found joint diagonalizer
#S=a structure containing some information about the run program:
#          S.iterations: number of iterations
#          S.LUerror: the LU error after each sweep
#          S.J2error: the J2 error after each sweep
#          S.J2RelativeError:the relative J2 error after each sweep
#BB=a three dimensional array containing the joint diagonalizer after each sweep
#Note: S and BB are not required as output
#
#This algorithm is based on a paper presented in ICA2006 conference and published in Springer LNCS
#Bijan Afsari, ''Simple LU and QR based Non-Orthogonal Matrix Joint Diagonalization''
##Coded by Bijan Afsari. Please forward any questions and problem to bijan@glue.umd.edu
#v.1.1
#Acknowledgements: Some data structures and implementation ideas in this code are inspired from the code for JADE
#written by J.F. Cardoso and from the code FFDIAG written by A. Ziehe
#Disclaimer: This code is to be used only for non-commercial research purposes and the author does not
#accept any reponsibility about its performance or fauilure
n=nrow(X)
m=ncol(X)
N=m/n

###
#MODE='B';
#if nargin==0, display('you must enter the data'); B=eye(n); return; end;
#if nargin==1, Err=ERR;Rbalance=RBALANCE;end;
#if nargin> 1, 
#MODE=upper(varargin{1});
#   switch MODE
#   case {'B'} 
#      ERR=varargin{2}; mflag='D'; if ERR >= 1, disp('Error value should be much smaller than unity');B=[];S=[]; return; end;
#   case ('E')
#      ERR=varargin{2};mflag='E'; if ERR >=1, disp('Error value should be much smaller than unity'); B=[];S=[];return;end;
#   case ('N');mflag='N'; ITER=varargin{2}; ERR=0; if ITER <= 1, disp('Number of itternations should be higher than one');B=[];S=[];return;end;
#   end
#end;
#if nargin==4, RBALANCE=varargin{3}; if ceil(RBALANCE)~=RBALANCE | RBALANCE<1, disp('RBALANCE should be a positive integer');B=[];S=[];return;end;end;
BB=array(NA,c(3,3,1))
Xar=array(NA,c(n,m,1))
JJ=c()
EERR=c()
EERRJ2=c()  
X1=X
B=diag(1,n)
Binv=diag(1,n)
J=0
JJ=0

for (t in 1:N){
   cols <- seq(((t-1)*n+1),(t*n))	
   J=J+norm(X1[,cols]-diag(diag(X[,cols])),"f")^2
JJ=cbind(JJ,J)
  }
#err=10^-3;
#the following part implements a sweep 
##########################
err=ERR*n+1
#if (mode=='B') ERR=ERR*n
   k=0
   while (err>ERR & k<ITER){
     
      k=k+1
      
      L=diag(1,n) #Linv=eye(n)
      U=diag(1,n) #Uinv=eye(n)
      Dinv=diag(1,n)
      for (i in 1:(n-1)){
         for (j in (i+1):n){
            cindex=seq(1,m)
            rem<-seq(i,m,by=n)
			cindex=cindex[-rem]			
            a=-(X[i,cindex]\%*\%X[j,cindex])/(X[j,cindex]\%*\%X[j,cindex])
			a=as.numeric(a)
			#a=tanh(a);
            if (abs(a)>1) a=sign(a)*1
            X[i,]=a*X[j,]+X[i,]
            I=seq(i,m,n)
            J=seq(j,m,n)
            X[,I]=a*X[,J]+X[,I]
            #U=eye(n,n);
            #U(i,j)=a;
            #B=U*B;
            U[i,]=U[i,]+a*U[j,]
            #Uinv(i,:)=Uinv(i,:)-a*Uinv(j,:);
            #Now the L phase
         }
     }
      
      for (i in 1:(n-1)){
         #for j=i+1:n
         rindex=c()
         Xj=c()
         for(j in (i+1):n){
            cindex=seq(1,m)
            rem<-seq(j,m,by=n)
			cindex=cindex[-rem]
            a=-(X[i,cindex]\%*\%X[j,cindex])/(X[i,cindex]\%*\%X[i,cindex])
			#coorelation quefficient
            #a=-(X(i,cindex)*X(j,cindex)')/(norm(X(i,cindex))*norm(X(j,cindex)));
            #a=tanh(a);
			a<- as.numeric(a)
            if (abs(a)>1) a=sign(a)*1
            X[j,]=a*X[i,]+X[j,]
            I=seq(i,m,by=n)
            J=seq(j,m,by=n)
            X[,J]=a*X[,I]+X[,J]
            L[j,]=L[j,]+a*L[i,]
            #Linv(j,:)=Linv(j,:)-a*Linv(i,:);
         }
     }
      
      
      B=L\%*\%U\%*\%B
	  #Binv=Binv*Uinv*Linv;
      #err=norm(L*U-eye(n,n),'fro');
      err=max(max(abs(L\%*\%U-diag(1,n))))
	  EERR=c(EERR,err)
      if (k\%\%RBALANCE==0){	#k\%\%RBALANCE gets remainder of k/RBALANCE
         d=rowSums(abs(X))
         D=diag(1/(d*N))
		 Dinv=diag(d*N)
         J=0
         for(t in 1:N){
            X[,seq(((t-1)*n+1),t*n)]=D\%*\%X[,seq(((t-1)*n+1),t*n)]\%*\%D
         }
         B=D\%*\%B; #Binv=Binv*Dinv;
      }
      J=0
	  BB=array(BB,c(n,n,k))
      BB[,,k]=B

	  # If B is singular, we provide the user with two options to invert the matrix:
	  #		1) Use the Moore-Penrose pseudo inverse instead.
	  #		2) Invert B with the shrunken eigenvalues of B with a small shrinkage factor.
      Binv <- try(solve(B), silent = TRUE)
	  if(!is.matrix(Binv)) {
	  	if(pseudo == TRUE) {
			Binv <- ginv(B)
		} else {
			B.eigen <- eigen(B, symmetric = TRUE)
			identity.shrink <- diag(B.eigen$values) + shrink.val * diag(length(B.eigen$vectors))
			Binv <- solve(B.eigen$vectors \%*\% identity.shrink \%*\% t(B.eigen$vectors))
		}
	  }


      for (t in 1:N){
	     cols <- seq(((t-1)*n+1),t*n)
	     J=J+norm(X1[,cols]-Binv\%*\%diag(diag(X[,cols]))\%*\%Binv,"f")^2
         JJ=c(JJ,J)
	  }
      #if MODE=='E', err=abs(JJ(end-1)-JJ(end))/JJ(end-1);EERRJ2=[EERRJ2,err];end
     
	 Xar=array(Xar,c(n,m,k))
     Xar[,,k]=X  
  }
   Y=X
  # return(Y)
return(B)
      #S=struct('iterations',k,'LUerror',EERR,'J2error',JJ,'J2RelativeError',EERRJ2);varargout{1}=S;varargout{2}=BB
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
