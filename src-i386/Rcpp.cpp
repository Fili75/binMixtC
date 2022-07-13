#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>


int whichmaxC(DoubleVector x)
{
  int i;
  int lenx=Rf_length(x);
  int index=0;
  double maxi=x(0);
  for(i=1;i<lenx;i++)
  {
    if(x(i)>maxi) 
    {maxi=x(i);
      index=i;
    }
  }
  return index;
}
// [[Rcpp::export]]
double logsumC(DoubleVector x) 
{
  int l=whichmaxC(x);
  int i;
  
  int lenx=Rf_length(x);
  DoubleVector x1;
  x1=clone(x);
  for(i=0;i<lenx;i++)
  {
    x1(i)=x1(i)-x(l);
  }
  x1=exp(x1);
  for(i=0;i<lenx;i++)
  {
    if(x1(i)==1) x1(i)=0;
  }
  double ret;
  ret=x(l)+log(1+sum(x1));
  return ret;
}

// [[Rcpp::export]]
double logdiffC(DoubleVector x) 
{
  int l=whichmaxC(x);
  int i;
  
  int lenx=Rf_length(x);
  DoubleVector x1;
  x1=clone(x);
  for(i=0;i<lenx;i++)
  {
    x1(i)=x1(i)-x(l);
  }
  x1=exp(x1);
  for(i=0;i<lenx;i++)
  {
    if(x1(i)==1) x1(i)=0;
  }
  double ret;
  ret=x(l)+log(1-sum(x1));
  return ret;
}

#include <random>
#include <array>


// [[Rcpp::export]]
DoubleVector ran(int n, DoubleVector p,DoubleVector mu, DoubleVector v)
{
  int lenp=Rf_length(p);
  DoubleVector ind1(lenp);
  DoubleVector ind(n);
  DoubleVector ind2(n);
  int ind3;
  int i;
  double mu3;
  double v3;
  DoubleVector a;
  for(i=0; i<lenp;i++)
    ind1(i)=i;
  ind=Rcpp::sample(ind1,n,true,p);
  for(i=0;i<n;i++)
  {
    ind3=ind(i);
    mu3=mu(ind3);
    v3=sqrt(v(ind3));
    a=Rcpp::rnorm(1,mu3,v3);
    ind2(i)=a(0);
  }
  return ind2;
}


using namespace arma; 
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

class Bin {
public:
  vec bin;
  vec grid;
  Bin(vec x, vec y);
}; 

Bin::Bin(vec x, vec y) { 
  bin=x; 
  grid = y; 
}




// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = randn(n, ncols);
  arma::mat Z(n,ncols);
  Z=Z.t();
  int i;
  for(i=0;i<n;i++)
    Z.col(i)=mu;
  Z=Z.t();
  return Y * chol(sigma)+Z;
}

// [[Rcpp::export]]
List genmixtC(int n, DoubleVector p,List mu, List sigma) {
  int lenp=Rf_length(p);
  arma::mat Z=sigma[0];
  int ncols=Z.n_cols;
  DoubleVector ind1(lenp);
  DoubleVector ind(n);
  arma::mat Y(n,ncols);
  int i;
  for(i=0; i<lenp;i++)
    ind1(i)=i;
  ind=Rcpp::sample(ind1,n,true,p);
  int j;
  for(i=0;i<n;i++)
  {
    j=ind(i);
    Y.row(i)=mvrnormArma(1, mu[j], sigma[j]);
  }
  return List::create(_["x"]=Y,_["ind"]=ind);
}

// [[Rcpp::export]]
List zerobin1(List bin)
{
  arma::vec gr1=bin[1];
  arma::vec bi=bin[0];
  int lenbin=bi.n_elem;
  int dimgr=gr1.n_elem;
  arma::mat gr(dimgr+2,2);
  int i;
  int j;
  int conta=0;
  for(i=0;i<lenbin;i++)
  {
    if(bi(i)>0) conta=conta+1;
  }
  gr(0,0)=-INFINITY;
  gr(dimgr+1,1)=INFINITY;
  for(i=1;i<dimgr+1;i++)
    gr(i,0)=gr1(i-1);
  for(i=0;i<dimgr;i++)
    gr(i,1)=gr1(i);
  
  arma::mat gr2(conta,2);
  arma::vec bizero(conta);
  conta=0;
  for(i=0;i<lenbin;i++)
  {
    if(bi(i)>0)
    {
      gr2(conta,0)=gr(i,0);
      gr2(conta,1)=gr(i,1);
      bizero(conta)=bi(i);
      conta=conta+1;
    }
  }
  return List::create(_["bin"]=bizero,_["grid"]=gr2);
  
}

// [[Rcpp::export]]
List buildbinC(arma::vec x, int R){
  double maxim=max(x);
  double minim=min(x);
  int len=x.n_elem;
  double interval= (maxim-minim)/(R-1);
  arma::vec z(R);
  z(0)=minim;
  z(R-1)=maxim;
  int i;
  int j;
  for(i=1;i<R-1;i++)
  {
    z(i)=minim+interval*(i);
  }
  arma::vec y(R+1);
  y(R)=0;
  for(i=0;i<R;i++)
  {
    for(j=0;j<len;j++)
    {
      if(i==0)
      {
        if(x(j)<=z(i)) y(i)=y(i)+1;
      } else {
        if( (x(j)<=z(i)) & (x(j)>z(i-1)))
          y(i)=y(i)+1;}
    }
  }
  return zerobin1(List::create(_["bin"]=y,_["grid"]=z));
}

// [[Rcpp::export]]
List buildbinmargC(arma::mat X, arma::vec R)
{
  int i;
  int ncol=X.n_cols;
  List a(ncol);
  for(i=0;i<ncol;i++)
  {
    a[i]=buildbinC(X.col(i),R(i));
  }
  return a;
}

// [[Rcpp::export]]
arma::vec dmixtC(arma::vec x, arma::vec pi, arma::vec mu, arma::vec v)
{
  int len=x.n_elem;
  int lenp=pi.n_elem;
  int i;
  int j;
  arma::vec y(len);
  double z;
  for(i=0;i<len;i++)
  {
    z=0;
    for(j=0;j<lenp;j++)
      z=R::dnorm(x(i),mu(j),sqrt(v(j)),FALSE)*pi(j)+z;
    y(i)=z;
  }
  return y;
}

// [[Rcpp::export]]
double loglimixtC(arma::vec x, arma::vec pi, arma::vec mu, arma::vec v)
{
  int len=x.n_elem;
  int i;
  arma::vec y=dmixtC(x,pi,mu,v);
  double z=0;
  for(i=0;i<len;i++)
  {
    z=z+log(y(i));
  }
  return z;
}

// [[Rcpp::export]]
double loglimultC(Rcpp::List bin, arma::vec pi, arma::vec mu, arma::vec v)
{
  arma::vec bi=bin[0];
  arma::mat gr=bin[1];
  int lenbin=bi.n_elem;
  int lenp=pi.n_elem;
  int j;
  int i;
  arma::vec y(lenbin);
  double z;
  double x;
  for(i=0;i<lenbin;i++)
  {
    z=0;
    for(j=0;j<lenp;j++)
      z=(R::pnorm(gr(i,1),mu(j),sqrt(v(j)),TRUE,FALSE)-
        R::pnorm(gr(i,0),mu(j),sqrt(v(j)),TRUE,FALSE))*pi(j)+z;
    if(z==0)
    {
      for(j=0;j<lenp;j++)
        z=(-R::pnorm(gr(i,1),mu(j),sqrt(v(j)),FALSE,FALSE)+
          R::pnorm(gr(i,0),mu(j),sqrt(v(j)),FALSE,FALSE))*pi(j)+z;
    }
    y(i)=z;
  }
  x=0;
  for(i=0;i<lenbin;i++)
  {
    x=x+bi(i)*log(y(i));
  }
  return x;
}

// [[Rcpp::export]]
List emuniv_nomixtC(List bin,arma::vec mu0,arma::vec sigma0,double eps,int it) 
{
  arma::vec bi=bin[0];
  arma::mat gr=bin[1];
  arma::vec logli(it+1);
  int lenbin=bi.n_elem;
  int j;
  int n=accu(Rcpp::as<arma::vec>((bin[0])));
  arma::vec s(lenbin);
  arma::vec mu1(1);
  arma::vec sigma1(1);
  arma::vec pi(1);
  pi(0)=1;
  logli(0)=loglimultC(bin,pi,mu0,sigma0);
  double z1;
  double z2;
  double z;
  int i;
  int flag=1;
  int k=1;
  while(flag==1)
  {
    for(i=0;i<lenbin;i++)
    {
      z1=(R::pnorm(gr(i,1),mu0(0),sqrt(sigma0(0)),TRUE,FALSE)-
        R::pnorm(gr(i,0),mu0(0),sqrt(sigma0(0)),TRUE,FALSE));
      
      z2=(-R::pnorm(gr(i,1),mu0(0),sqrt(sigma0(0)),FALSE,FALSE)+
        R::pnorm(gr(i,0),mu0(0),sqrt(sigma0(0)),FALSE,FALSE));
      if(z1>=z2)
      {s(i)=z1;} else
      {s(i)=z2;}
    }
    arma::vec sm(lenbin);
    for(i=0;i<lenbin;i++)
    {
      z=(R::dnorm(gr(i,1),mu0(0),sqrt(sigma0(0)),FALSE)-
        R::dnorm(gr(i,0),mu0(0),sqrt(sigma0(0)),FALSE));
      sm(i)=z;
    }
    arma::vec sv(lenbin);
    
    for(i=0;i<lenbin;i++)
    {
      if(gr(i,0)==-INFINITY)
      {
        z=(R::dnorm(gr(i,1),mu0(0),sqrt(sigma0(0)),FALSE))*gr(i,1);
        sv(i)=z;
      } else{
        if(gr(i,1)==INFINITY)
        {
          
          z=-(R::dnorm(gr(i,0),mu0(0),sqrt(sigma0(0)),FALSE))*gr(i,0);
          sv(i)=z;
          
        } else {
          z=(R::dnorm(gr(i,1),mu0(0),sqrt(sigma0(0)),FALSE))*gr(i,1)-
            (R::dnorm(gr(i,0),mu0(0),sqrt(sigma0(0)),FALSE))*gr(i,0);
          sv(i)=z;
        }
      }
    }  
    
    arma::vec intm(lenbin);
    for(i=0;i<lenbin;i++)
      intm(i)=mu0(0)-sigma0(0)*(sm(i)/s(i));
    arma::vec sum2=intm.t()*bi;
    mu1=sum2/n;
    arma::vec intv(lenbin);
    for(i=0;i<lenbin;i++)
      intv(i)=(sigma0(0)*(s(i)+(2*mu1(0)-mu0(0))*sm(i)-sv(i))+s(i)*(mu1(0)-mu0(0))*(mu1(0)-mu0(0)))/s(i);
    arma::vec sum3=intv.t()*bi;
    sigma1=sum3/n;
    logli(k)=loglimultC(bin,pi,mu1,sigma1);
    if(Rcpp::traits::is_nan<REALSXP>(logli(k))==true) 
    {return List::create(_["logli"]=-INFINITY);}
    flag=0;
    if ((abs((logli(k) - logli(k - 1))/logli(k - 1)) > eps) & 
        (k < it)) {
      flag = 1;
      mu0 = mu1;
      sigma0= sigma1;
      k=k+1;
    }
  }
  arma::vec logli1(k+1);
  for(i=0;i<=k;i++)
    logli1(i)=logli(i);
  return List::create(_["bin"]=bin,_["mu"]=mu1,_["v"]=sigma1,_["logli"]=logli1);
}

// [[Rcpp::export]]
List emunivC(List bin,arma::vec pi0,arma::vec mu0,arma::vec sigma0,double eps,int it) 
{
  arma::vec bi=bin[0];
  arma::mat gr=bin[1];
  arma::vec logli(it+1);
  int lenbin=bi.n_elem;
  int lenp=pi0.n_elem;
  int j;
  int n=accu(Rcpp::as<arma::vec>((bin[0])));
  arma::mat s(lenbin,lenp);
  arma::vec pi1(lenp);
  arma::vec mu1(lenp);
  arma::vec sigma1(lenp);
  logli(0)=loglimultC(bin,pi0,mu0,sigma0);
  double z1;
  double z2;
  double z;
  int i;
  int flag=1;
  int k=1;
  while(flag==1)
  {
    for(i=0;i<lenbin;i++)
    {
      for(j=0;j<lenp;j++)
      {
        z1=(R::pnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),TRUE,FALSE)-
          R::pnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),TRUE,FALSE));
        
        z2=(-R::pnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE,FALSE)+
          R::pnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE,FALSE));
        if(z1>=z2)
        {s(i,j)=z1;} else
        {s(i,j)=z2;}
      }
    }
    arma::vec s1(lenbin);
    s1=s*pi0;
    arma::mat sm(lenbin,lenp);
    for(i=0;i<lenbin;i++)
    {
      for(j=0;j<lenp;j++)
      {
        z=(R::dnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE)-
          R::dnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE));
        sm(i,j)=z;
      }
    }
    arma::mat sv(lenbin,lenp);
    
    for(i=0;i<lenbin;i++)
    {
      if(gr(i,0)==-INFINITY)
      {
        for(j=0;j<lenp;j++)
        {
          z=(R::dnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,1);
          sv(i,j)=z;
        }
      } else{
        if(gr(i,1)==INFINITY)
        {
          for(j=0;j<lenp;j++)
          { 
            z=-(R::dnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,0);
            sv(i,j)=z;
          }
        } else {
          for(j=0;j<lenp;j++)
          {
            z=(R::dnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,1)-
              (R::dnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,0);
            sv(i,j)=z;
          }
        }
      }
    }  
    arma::mat inti(lenbin,lenp);
    for(i=0;i<lenbin;i++)
      for(j=0;j<lenp;j++)
        inti(i,j)=s(i,j)*pi0(j)/s1(i);
    arma::mat intm(lenbin,lenp);
    for(i=0;i<lenbin;i++)
      for(j=0;j<lenp;j++)
        intm(i,j)=(pi0(j)/s1(i))*((s(i,j)*mu0(j))-(sm(i,j)*sigma0(j)));
    
    arma::vec sum1=inti.t()*bi;
    pi1=sum1/n;
    arma::vec sum2=intm.t()*bi;
    mu1=sum2/sum1;
    arma::mat intv(lenbin,lenp);
    for(i=0;i<lenbin;i++)
      for(j=0;j<lenp;j++)
        intv(i,j)=(pi0(j)/s1(i))*((s(i,j)+sm(i,j)*(2*mu1(j)-mu0(j))-sv(i,j))*sigma0(j)+s(i,j)*(mu1(j)-mu0(j))*(mu1(j)-mu0(j)));
    arma::vec sum3=intv.t()*bi;
    sigma1=sum3/sum1;
    logli(k)=loglimultC(bin,pi1,mu1,sigma1);
    if(Rcpp::traits::is_nan<REALSXP>(logli(k))==true) 
    {return List::create(_["logli"]=-INFINITY);}
    flag=0;
    if ((abs((logli(k) - logli(k - 1))/logli(k - 1)) > eps) & 
        (k < it)) {
      flag = 1;
      pi0 = pi1;
      mu0 = mu1;
      sigma0 = sigma1;
      k=k+1;
    }
  }
  arma::vec logli1(k+1);
  for(i=0;i<=k;i++)
    logli1(i)=logli(i);
  return List::create(_["bin"]=bin,_["pi"]=pi1,_["mu"]=mu1,_["v"]=sigma1,_["logli"]=logli1);
}

// [[Rcpp::export]]
arma::vec margmu(List mu,int dim)
{
  int listdim=mu.size();
  arma::vec mu1=mu[0];
  int dime=mu1.n_elem;
  int i;
  arma::vec m(listdim);
  arma::vec a(dime);
  for(i=0;i<listdim;i++)
  {
    a=(Rcpp::as<arma::vec>(mu[i]));
    m(i)=a(dim);
  }
  return m;
}

// [[Rcpp::export]]
arma::vec margv(List v,int dim)
{
  int i;
  int listdim=v.size();
  arma::mat v2=v[0];
  int dime=v2.n_rows;
  arma::vec v1(listdim);
  arma::mat a(dime,dime);
  for(i=0;i<listdim;i++)
  {
    a=(Rcpp::as<arma::mat>(v[i]));
    v1(i)=a(dim,dim);
  }
  return v1;
}

// [[Rcpp::export]]
double loglimargC(Rcpp::List bin, arma::vec pi, List mu, List v)
{
  arma::vec mu1=mu[0];
  int dim=mu1.n_elem;
  double z=0;
  int i;
  for(i=0;i<dim;i++)
  {
    arma::vec m=margmu(mu,i);
    arma::vec v0=margv(v,i);
    z=z+loglimultC(Rcpp::as<List>(bin)[i],pi,m,v0);
  }
  return z;
}
// [[Rcpp::export]]
List emdimC(List bin,arma::vec pi0,arma::vec mu0,arma::vec sigma0) 
{
  arma::vec bi=bin[0];
  arma::mat gr=bin[1];
  int lenbin=bi.n_elem;
  int lenp=pi0.n_elem;
  int j;
  arma::mat s(lenbin,lenp);
  double z1;
  double z2;
  double z;
  int i;
  for(i=0;i<lenbin;i++)
  {
    for(j=0;j<lenp;j++)
    {
      z1=(R::pnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),TRUE,FALSE)-
        R::pnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),TRUE,FALSE));
      
      z2=(-R::pnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE,FALSE)+
        R::pnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE,FALSE));
      if(z1>=z2)
      {s(i,j)=z1;} else
      {s(i,j)=z2;}
    }
  }
  arma::vec s1(lenbin);
  s1=s*pi0;
  arma::mat sm(lenbin,lenp);
  for(i=0;i<lenbin;i++)
  {
    for(j=0;j<lenp;j++)
    {
      z=(R::dnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE)-
        R::dnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE));
      sm(i,j)=z;
    }
  }
  arma::mat sv(lenbin,lenp);
  
  for(i=0;i<lenbin;i++)
  {
    if(gr(i,0)==-INFINITY)
    {
      for(j=0;j<lenp;j++)
      {
        z=(R::dnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,1);
        sv(i,j)=z;
      }
    } else{
      if(gr(i,1)==INFINITY)
      {
        for(j=0;j<lenp;j++)
        { 
          z=-(R::dnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,0);
          sv(i,j)=z;
        }
      } else {
        for(j=0;j<lenp;j++)
        {
          z=(R::dnorm(gr(i,1),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,1)-
            (R::dnorm(gr(i,0),mu0(j),sqrt(sigma0(j)),FALSE))*gr(i,0);
          sv(i,j)=z;
        }
      }
    }
  }  
  arma::mat inti(lenbin,lenp);
  for(i=0;i<lenbin;i++)
    for(j=0;j<lenp;j++)
      inti(i,j)=s(i,j)*pi0(j)/s1(i);
  arma::mat intm(lenbin,lenp);
  for(i=0;i<lenbin;i++)
    for(j=0;j<lenp;j++)
      intm(i,j)=(pi0(j)/s1(i))*((s(i,j)*mu0(j))-(sm(i,j)*sigma0(j)));
  
  arma::vec sum1=inti.t()*bi;
  arma::vec sum2=intm.t()*bi;
  arma::vec mu1=sum2/sum1;
  arma::mat intv(lenbin,lenp);
  for(i=0;i<lenbin;i++)
    for(j=0;j<lenp;j++)
      intv(i,j)=(pi0(j)/s1(i))*((s(i,j)+sm(i,j)*(2*mu1(j)-mu0(j))-sv(i,j))*sigma0(j)+s(i,j)*(mu1(j)-mu0(j))*(mu1(j)-mu0(j)));
  arma::vec sum3=intv.t()*bi;
  arma::vec sigma1=sum3/sum1;
  return List::create(_["mu"]=mu1,_["v"]=sigma1,_["sum1"]=sum1);
}
// [[Rcpp::export]]
List embingauscompC(List bin,arma::vec pi0,List mu0,List sigma0,double eps, int it)
{
  int flag=1;
  int n=accu(Rcpp::as<arma::vec>(Rcpp::as<List>(bin[0])[0]));
  int j=0;
  int linp=pi0.n_elem;
  int dim=Rcpp::as<arma::vec>(mu0[0]).n_elem;
  List mu1(linp);
  List sigma1(linp);
  arma::vec pi1(linp);
  int k;
  arma::vec logli(it+1);
  logli(0)=loglimargC(bin,pi0,mu0,sigma0);
  List dim1(dim);
  arma::vec mu2(dim);
  arma::vec smarg(dim);
  arma::vec mmarg(dim);
  arma::mat s2(dim,dim);
  int i;
  while(flag==1) {
    for(i=0;i<dim;i++)
    {
      smarg=margv(sigma0,i);
      mmarg=margmu(mu0,i);
      dim1[i]=emdimC(bin[i],pi0,mmarg,smarg);
    }
    for(k=0;k<linp;k++)
    {
      pi1(k)=0;
      for(i=0;i<dim;i++)
      {
        pi1(k)=Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[i])[2])(k)+pi1(k);
      }
    }
    for(k=0;k<linp;k++) pi1(k)=pi1(k)/(dim*n);
    for(k=0;k<linp;k++)
    {
      for(i=0;i<dim;i++)
      {
        (mu2)(i)=Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[i])[0])(k);
        (s2)(i,i)=Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[i])[1])(k);
      }
      mu1[k]=mu2;
      sigma1[k]=s2;
    }
    j = j + 1;
    logli(j)=loglimargC(bin,pi1,mu1,sigma1);
    if(Rcpp::traits::is_nan<REALSXP>(logli(j))==true) 
    {return List::create(_["logli"]=-INFINITY);}
    flag=0;
    if ((abs((logli(j) - logli(j - 1))/logli(j - 1)) > eps) & 
        (j < it)) {
      flag = 1;
      pi0 = pi1;
      mu0 = mu1;
      sigma0 = sigma1;
    }
  }
  arma::vec logli1(j+1);
  for(i=0;i<=j;i++)
    logli1(i)=logli(i);
  return List::create(_["bin"]=bin,_["pi"]=pi1,_["mu"]=mu1,_["v"]=sigma1,_["logli"]=logli1);
}


// [[Rcpp::export]]
List binmixtC(arma::mat data, int cl, arma::vec R, int it, double eps, int nrep)
{
  int dim=data.n_cols;
  Rcpp::Rcout<<dim<< std::endl;
  
  if(dim==1)
  {
    List bin=buildbinC(data,R(0));
    int i;
    arma::vec pi0(cl);
    arma::vec mu0(cl);
    arma::vec s0(cl);
    double s;
    int k;
    int j;
    Progress p(nrep, true);
    List e(nrep);
    arma::vec logli(nrep);
    double minimo;
    double massimo;
    double varianza;
    minimo=arma::min(data.col(0));
    massimo=arma::max(data.col(0));
    varianza=arma::var(data.col(0));
    for(i=0;i<nrep;i++)
    {
      pi0=arma::randu(cl);
      s=0;
      for(k=0;k<cl;k++)
        s=s+pi0(k);
      pi0=pi0/s;
      for(k=0;k<cl;k++)
      {
        double random=arma::randu();
        mu0(k)=(random*(massimo))+(1-random)*minimo;
      }
      for(k=0;k<cl;k++)
      {
        s0(k)=arma::randu()*varianza;
      }
      List e1;
      if(cl==1)
      {
        e1=emuniv_nomixtC(bin,mu0,s0,eps,it); 
      } else
      {
        e1=emunivC(bin,pi0,mu0,s0,eps,it);
      }
      e[i]=e1;
      logli(i)=Rcpp::as<arma::vec>((e1["logli"])).max();
      p.increment(); 
    }
    double logli1=-INFINITY;
    int indice=0;
    for(i=0;i<nrep;i++)
    {
      if(logli(i)>logli1)
      {
        logli1=logli(i);
        indice=i;
      }
    }
    return e[indice];
  } else {
    List bin=buildbinmargC(data,R);
    int i;
    int dim=data.n_cols;
    arma::vec pi0(cl);
    List mu0(cl);
    List s0(cl);
    arma::vec mu01(dim);
    arma::mat s01(dim,dim);
    double s;
    int k;
    int j;
    Progress p(nrep, true);
    List e(nrep);
    arma::vec logli(nrep);
    arma::vec minimo(dim);
    arma::vec massimo(dim);
    arma::vec varianza(dim);
    for(i=0;i<dim;i++)
    {
      minimo(i)=arma::min(data.col(i));
      massimo(i)=arma::max(data.col(i));
      varianza(i)=arma::var(data.col(i));
    }
    for(i=0;i<nrep;i++)
    {
      pi0=arma::randu(cl);
      s=0;
      for(k=0;k<cl;k++)
        s=s+pi0(k);
      pi0=pi0/s;
      for(k=0;k<cl;k++)
      {
        for(j=0;j<dim;j++)
        {
          double random=arma::randu();
          mu01(j)=(random*(massimo(j)))+(1-random)*minimo(j);
        }
        mu0[k]=mu01;
      }
      for(k=0;k<cl;k++)
      {
        for(j=0;j<dim;j++)
          s01(j,j)=arma::randu()*varianza(j);
        s0[k]=s01;
      }
      List e1;
      e1=embingauscompC(bin,pi0,mu0,s0,eps,it);
      e[i]=e1;
      logli(i)=Rcpp::as<arma::vec>((e1["logli"])).max();
      p.increment(); 
    }
    double logli1=-INFINITY;
    int indice=0;
    for(i=0;i<nrep;i++)
    {
      if(logli(i)>logli1)
      {
        logli1=logli(i);
        indice=i;
      }
    }
    return e[indice];
  }
}


// [[Rcpp::export]]
List binmixtclassicC(arma::mat data, int cl, arma::vec R, int it, double eps, double eps1, int it1,int nrep)
{
  int dim=data.n_cols;
  if(dim==1)
  {
    Rcpp::Rcout<<"Choosing the best initial guess..."<< std::endl;
    
    List migliore=binmixtC(data,cl,R,it,eps,nrep);
    Rcpp::Rcout<<"Estimation phase...";
    if(cl==1)
    {
      List e=emuniv_nomixtC(migliore["bin"],migliore["mu"],migliore["v"],eps1,it1);
      Rcpp::Rcout<<"Done."<< std::endl;
      return e;
    } else
    {
      List e=emunivC(migliore["bin"],migliore["pi"],
                     migliore["mu"],migliore["v"],eps1,it1);
      Rcpp::Rcout<<"Done."<< std::endl;
      return e;
    }
  }
  else
  {
    Rcpp::Rcout<<"Choosing the best initial guess..."<< std::endl;
    
    List migliore=binmixtC(data,cl,R,it,eps,nrep);
    Rcpp::Rcout<<"Estimation phase...";
    
    List e=embingauscompC(migliore["bin"],migliore["pi"],
                          migliore["mu"],migliore["v"],eps1,it1);
    Rcpp::Rcout<<"Done."<< std::endl;
    return e;
  }
}

// [[Rcpp::export]]
NumericVector tabulate2(NumericVector x, int max) {
  NumericVector counts(max);
  std::size_t n = x.size();
  for (std::size_t i = 0; i < n; i++) {
    if (x(i) > 0 && x(i) <= max)
      counts(x(i) - 1)++;
  }
  return counts;
}

// [[Rcpp::export]]
List zerobin(List bin)
{
  arma::mat gr1=bin[1];
  arma::vec bi=bin[0];
  int lenbin=bi.n_elem;
  int dimgr=gr1.n_rows;
  arma::mat gr(dimgr+2,2);
  int i;
  int j;
  int conta=0;
  for(i=0;i<lenbin;i++)
  {
    if(bi(i)>0) conta=conta+1;
  }
  for(j=0;j<2;j++)
  {
    gr(0,j)=-INFINITY;
    gr(dimgr+1,j)=INFINITY;
    for(i=1;i<dimgr+1;i++)
      gr(i,j)=gr1(i-1,j);
  }
  int divisione;
  int resto;
  arma::mat gr2(conta,4);
  arma::vec bizero(conta);
  int R=sqrt(lenbin);
  conta=0;
  for(i=0;i<lenbin;i++)
  {
    divisione=i/R;
    resto=i-divisione*R;
    if(bi(i)>0)
    {
      gr2(conta,0)=gr(resto,0);
      gr2(conta,1)=gr(divisione,1);
      gr2(conta,2)=gr(resto+1,0);
      gr2(conta,3)=gr(divisione+1,1);
      bizero(conta)=bi(i);
      conta=conta+1;
    }
  }
  
  return List::create(_["bin"]=bizero,_["grid"]=gr2);
  
}

// [[Rcpp::export]]
List buildbinmultC(arma::mat x, int R){
  int dim=x.n_cols;
  int i;
  arma::mat b(R,dim);
  for(i=0;i<dim;i++)
  {
    double maxim=max(x.col(i));
    double minim=min(x.col(i));
    double interval= (maxim-minim)/(R-1);
    arma::vec z(R);
    b(0,i)=minim;
    b(R-1,i)=maxim;
    int j;
    for(j=1;j<R-1;j++)
    {
      b(j,i)=minim+interval*(j);
    }
  }
  int len=x.n_rows;
  arma::mat index(len,dim);
  int k;
  int j;
  for(k=0;k<dim;k++)
  {
    for(i=0;i<R;i++)
    {
      for(j=0;j<len;j++)
      {
        if( (x(j,k)>=b(i,k)))
          index(j,k)=index(j,k)+1;
      }
    }
  }
  arma::vec mult(dim);
  NumericVector n(len);
  for(i=0;i<dim;i++)
    mult(i)=std::pow(R+1,i);
  for(i=0;i<len;i++)
    for(j=0;j<dim;j++)
      n(i)=n(i)+index(i,j)*mult(j)+1*(j==0);
  NumericVector indici2(pow(R+1,dim));
  indici2=tabulate2(n,pow(R+1,dim));
  return zerobin(List::create(_["bin"]=indici2,_["grid"]=b));
}

// [[Rcpp::export]]
int choosecpp(int n, int k) {
  if(k == 0) return 1;
  return (n * choosecpp(n - 1, k - 1)) / k;
}

// [[Rcpp::export]]
arma::mat combncpp(int N, int K) {
  if(K > N) Rcpp::stop("K > N");
  std::string bitmask(K, 1);
  bitmask.resize(N, 0);
  
  int n_combos = choosecpp(N,K);
  arma::mat results(n_combos, K);
  uint64_t row_position = 0;
  do {
    uint64_t col_position = 0;
    for (int i = 0; i < N; ++i)  {
      if (bitmask[i]) {
        results(row_position, col_position) = i+1;
        col_position++;
      }
    }
    row_position++;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  return results;
}

// [[Rcpp::export]]
List buildbinbivmargC(arma::mat X, arma::vec R)
{
  int i;
  int ncol=X.n_cols;
  int comb=ncol*(ncol-1);
  comb=comb*0.5;
  arma::mat combn(comb,2);
  combn=combncpp(ncol,2);
  List a(comb);
  
  for(i=0;i<comb;i++)
  {
    a[i]=buildbinmultC(join_horiz(X.col(combn(i,0)-1),X.col(combn(i,1)-1)),R(i));
  }
  return a;
}

//[[Rcpp::export]]
double loglimultmbivC(Rcpp::List bin, arma::vec pi, Rcpp::List mu, Rcpp::List v)
{
  arma::vec bi=bin[0];
  arma::mat gr=bin[1];
  arma::vec mu0=mu[0];
  int dimgr=gr.n_rows;
  int lenbin=bi.n_elem;
  int lenp=pi.n_elem;
  int i;
  int j;
  double z;
  double x;
  double z1;
  double z2;
  arma::vec y(lenbin);
  for(i=0;i<lenbin;i++)
  {
    z=0;
    for(j=0;j<lenp;j++)
    {
      z1=(R::pnorm(gr(i,2),Rcpp::as<arma::vec>(mu[j])(0),sqrt(Rcpp::as<arma::mat>(v[j])(0,0)),TRUE,FALSE)-
        R::pnorm(gr(i,0),Rcpp::as<arma::vec>(mu[j])(0),sqrt(Rcpp::as<arma::mat>(v[j])(0,0)),TRUE,FALSE));
      if(z1==0)
      {
        z1=(-R::pnorm(gr(i,2),Rcpp::as<arma::vec>(mu[j])(0),sqrt(Rcpp::as<arma::mat>(v[j])(0,0)),FALSE,FALSE)+
          R::pnorm(gr(i,0),Rcpp::as<arma::vec>(mu[j])(0),sqrt(Rcpp::as<arma::mat>(v[j])(0,0)),FALSE,FALSE));
      }
      z2=(R::pnorm(gr(i,3),Rcpp::as<arma::vec>(mu[j])(1),sqrt(Rcpp::as<arma::mat>(v[j])(1,1)),TRUE,FALSE)-
        R::pnorm(gr(i,1),Rcpp::as<arma::vec>(mu[j])(1),sqrt(Rcpp::as<arma::mat>(v[j])(1,1)),TRUE,FALSE));
      if(z2==0)
      {
        z2=(-R::pnorm(gr(i,3),Rcpp::as<arma::vec>(mu[j])(1),sqrt(Rcpp::as<arma::mat>(v[j])(1,1)),FALSE,FALSE)+
          R::pnorm(gr(i,1),Rcpp::as<arma::vec>(mu[j])(1),sqrt(Rcpp::as<arma::mat>(v[j])(1,1)),FALSE,FALSE));
      }
      z=z1*z2*pi(j)+z;
    }
    y(i)=z;
  }
  x=0;
  for(i=0;i<lenbin;i++)
  {
    x=x+bi(i)*log(y(i));
  }
  return x;
}




// [[Rcpp::export]]
List margmubiv(List mu,arma::vec dim)
{
  int listdim=mu.size();
  arma::vec mu1=mu[0];
  int dime=mu1.n_elem;
  int i;
  int j;
  List m(listdim);
  arma::vec a(dime);
  arma::vec media(2);
  for(i=0;i<listdim;i++)
  {
    a=(Rcpp::as<arma::vec>(mu[i]));
    for(j=0;j<2;j++)
    {
      media(j)=a(dim(j));
    }
    m[i]=media;
  }
  return m;
}

// [[Rcpp::export]]
List margvbiv(List v,arma::vec dim)
{
  int listdim=v.size();
  arma::mat v2=v[0];
  int dime=v2.n_rows;
  int i;
  int j;
  List v1(listdim);
  arma::mat a(dime,dime);
  arma::mat var(2,2);
  for(i=0;i<listdim;i++)
  {
    a=(Rcpp::as<arma::mat>(v[i]));
    for(j=0;j<2;j++)
    {
      var(j,j)=a(dim(j),dim(j));
    }
    v1[i]=var;
  }
  return v1;
}


// [[Rcpp::export]]
double loglimargbivC(Rcpp::List bin, arma::vec pi, List mu, List v)
{
  arma::vec mu1=mu[0];
  int dim=mu1.n_elem;
  double z=0;
  int i;
  int comb=dim*(dim-1);
  comb=comb*0.5;
  arma::mat combn(comb,2);
  combn=combncpp(dim,2);
  arma::vec ro(2);
  int j;
  
  for(i=0;i<comb;i++)
  {
    for(j=0;j<2;j++)
    {
      ro(j)=combn(i,j)-1;
    }
    
    List m=margmubiv(mu,ro);
    List v0=margvbiv(v,ro);
    z=z+loglimultmbivC(Rcpp::as<List>(bin)[i],pi,m,v0);
  }
  return z;
}


// [[Rcpp::export]]
List emdimbivpimuC(Rcpp::List bin,arma::vec pi0,List mu0,List sigma0)
{
  arma::vec bi=bin[0];
  arma::mat gr=bin[1];
  int dimgr=gr.n_rows;
  arma::vec mu=mu0[0];
  int lenbin=bi.n_elem;
  int lenp=pi0.n_elem;
  int dim=mu.n_elem;
  int i;
  int j;
  double z;
  double x;
  arma::mat p11(lenbin,lenp);
  arma::mat p12(lenbin,lenp);
  
  for(i=0;i<lenbin;i++)
  {
    double z1;
    double z2;
    z=0;
    for(j=0;j<lenp;j++)
    {
      z1=(R::pnorm(gr(i,2),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),TRUE,FALSE)-
        R::pnorm(gr(i,0),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),TRUE,FALSE));
      z2=(-R::pnorm(gr(i,2),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE,FALSE)+
        R::pnorm(gr(i,0),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE,FALSE));
      if(z1>=z2)
      {p11(i,j)=z1;} else
      {p11(i,j)=z2;}
      z1=(R::pnorm(gr(i,3),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),TRUE,FALSE)-
        R::pnorm(gr(i,1),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),TRUE,FALSE));
      
      z2=(-R::pnorm(gr(i,3),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE,FALSE)+
        R::pnorm(gr(i,1),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE,FALSE));
      if(z1>=z2)
      {p12(i,j)=z1;} else
      {p12(i,j)=z2;}
    }
  }
  arma::mat p1(lenbin,lenp);
  for(i=0;i<lenbin;i++)
    for(j=0;j<lenp;j++)
      p1(i,j)=p11(i,j)*p12(i,j);
  arma::vec pcol(lenbin);
  pcol=p1*pi0;
  arma::mat p(lenbin,lenp);
  for(i=0;i<lenbin;i++)
    for(j=0;j<lenp;j++)
      p(i,j)=p1(i,j)*pi0(j)/pcol(i);
  arma::vec tau_pi(lenp);
  tau_pi=p.t()*bi;
  arma::mat a11(lenbin,lenp);
  arma::mat a12(lenbin,lenp);
  for(i=0;i<lenbin;i++)
  {
    double z1;
    double z2;
    z=0;
    for(j=0;j<lenp;j++)
    {
      a11(i,j)=(R::dnorm(gr(i,2),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE)-
        R::dnorm(gr(i,0),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE));
      a12(i,j)=(R::dnorm(gr(i,3),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE)-
        R::dnorm(gr(i,1),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE));
    }
  }
  arma::vec mean1(2);
  arma::vec intm1(lenbin);
  arma::vec intm2(lenbin);
  List mu1(lenp);
  for(j=0;j<lenp;j++)
  {
    for(i=0;i<lenbin;i++)
    {
      intm1(i)=(p12(i,j)*pi0(j)/pcol(i))*((p11(i,j)*Rcpp::as<arma::vec>(mu0[j])(0))-(a11(i,j)*Rcpp::as<arma::mat>(sigma0[j])(0,0)));
      intm2(i)=(p11(i,j)*pi0(j)/pcol(i))*((p12(i,j)*Rcpp::as<arma::vec>(mu0[j])(1))-(a12(i,j)*Rcpp::as<arma::mat>(sigma0[j])(1,1)));
    }
    mean1(0)=0;
    mean1(1)=0;
    for(i=0;i<lenbin;i++)
    {
      mean1(0)=intm1(i)*bi(i)+mean1(0);
      mean1(1)=intm2(i)*bi(i)+mean1(1);
    }
    mu1[j]=mean1;
  }
  arma::mat b11(lenbin,lenp);
  arma::mat b12(lenbin,lenp);
  for(i=0;i<lenbin;i++)
  {
    if(gr(i,0)==-INFINITY)
    { 
      for(j=0;j<lenp;j++)
      {
        z=(R::dnorm(gr(i,2),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE))*gr(i,2);
        b11(i,j)=z;
      }
    } else {
      if(gr(i,2)==INFINITY)
      { 
        for(j=0;j<lenp;j++)
        {
          z=-(R::dnorm(gr(i,0),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE))*gr(i,0);
          b11(i,j)=z;
        }
      } else { 
        for(j=0;j<lenp;j++)
        {
          z=gr(i,2)*(R::dnorm(gr(i,2),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE))-
            gr(i,0)* (R::dnorm(gr(i,0),Rcpp::as<arma::vec>(mu0[j])(0),sqrt(Rcpp::as<arma::mat>(sigma0[j])(0,0)),FALSE));
          b11(i,j)=z;
        }
      }
    }
    
    if(gr(i,1)==-INFINITY)
    { 
      for(j=0;j<lenp;j++)
      {
        z=(R::dnorm(gr(i,3),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE))*gr(i,3);
        b12(i,j)=z;
      }
    } else {
      if(gr(i,3)==INFINITY)
      { 
        for(j=0;j<lenp;j++)
        {
          z=-(R::dnorm(gr(i,1),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE))*gr(i,1);
          b12(i,j)=z;
        }
      } else {
        for(j=0;j<lenp;j++)
        {
          z=gr(i,3)*(R::dnorm(gr(i,3),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE))-
            gr(i,1)* (R::dnorm(gr(i,1),Rcpp::as<arma::vec>(mu0[j])(1),sqrt(Rcpp::as<arma::mat>(sigma0[j])(1,1)),FALSE));
          b12(i,j)=z;
        }
      }       
    }
    
    
  }
  
  return List::create(_["tau_pi"]=tau_pi,_["tau_mu"]=mu1,_["pcol"]=pcol,_["p11"]=p11,_["p12"]=p12,_["a11"]=a11,_["a12"]=a12,_["b11"]=b11,_["b12"]=b12);
}
// [[Rcpp::export]]
List sigmadimfunC(List bin,arma::vec pcol,arma::mat p11,arma::mat p12,arma::mat a11,
                  arma::mat a12,arma::mat b11,arma::mat b12,arma::vec pi0,List mu0,List mu1,List sigma0)
{
  int linp=pi0.n_elem;
  List sigma1(linp); 
  int len=p11.n_rows;
  int i;
  int j;
  arma::mat var1(2,len);
  for(j=0;j<linp;j++)
  {
    
    for(i=0;i<len;i++)
    {
      var1(0,i)=(p12(i,j)*pi0(j)/pcol(i))*(Rcpp::as<arma::mat>(sigma0[j])(0,0)*(p11(i,j)+a11(i,j)*(2*Rcpp::as<arma::vec>(mu1[j])(0)-Rcpp::as<arma::vec>(mu0[j])(0))-b11(i,j))+p11(i,j)*std::pow(Rcpp::as<arma::vec>(mu1[j])(0)-Rcpp::as<arma::vec>(mu0[j])(0),2));
      var1(1,i)=(p11(i,j)*pi0(j)/pcol(i))*(Rcpp::as<arma::mat>(sigma0[j])(1,1)*(p12(i,j)+a12(i,j)*(2*Rcpp::as<arma::vec>(mu1[j])(1)-Rcpp::as<arma::vec>(mu0[j])(1))-b12(i,j))+p12(i,j)*std::pow(Rcpp::as<arma::vec>(mu1[j])(1)-Rcpp::as<arma::vec>(mu0[j])(1),2));
    }           
    sigma1[j]=var1*Rcpp::as<arma::vec>(bin[0]);
    
  }
  return(sigma1);
}

// [[Rcpp::export]]
List embinbivcompC(List bin,arma::vec pi0,List mu0,List sigma0,double eps, int it)
{
  int flag=1;
  int n=accu(Rcpp::as<arma::vec>(Rcpp::as<List>(bin[0])[0]));
  int j=0;
  int linp=pi0.n_elem;
  int dim=Rcpp::as<arma::vec>(mu0[0]).n_elem;
  int comb=dim*(dim-1);
  comb=comb*0.5;
  arma::mat combn(comb,2);
  combn=combncpp(dim,2);
  List mu1(linp);
  List sigma1(linp);
  arma::vec pi1(linp);
  int k;
  arma::vec logli(it+1);
  logli(0)=loglimargbivC(bin,pi0,mu0,sigma0);
  List dim1(comb);
  arma::vec mu2(dim);
  List smarg(comb);
  List mmarg(comb);
  arma::mat s2(dim,dim);
  int i;
  arma::vec ro(2);
  
  while(flag==1) {
    for(i=0;i<comb;i++)
    {
      for(k=0;k<2;k++)
      {
        ro(k)=combn(i,k)-1;
      }
      smarg=margvbiv(sigma0,ro);
      mmarg=margmubiv(mu0,ro);
      dim1[i]=emdimbivpimuC(bin[i],pi0,mmarg,smarg);
    }
    for(k=0;k<linp;k++)
    {
      pi1(k)=0;
      for(i=0;i<dim;i++)
      {
        pi1(k)=Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[i])[0])(k)+pi1(k);
      }
    }
    arma::vec m(dim);
    arma::mat mu1_col(dim,linp);
    for(k=0;k<linp;k++)
    {
      mu1_col(0,k)=(Rcpp::as<arma::vec>(Rcpp::as<List>(Rcpp::as<List>(dim1[0])[1])[k])(0)+Rcpp::as<arma::vec>(Rcpp::as<List>(Rcpp::as<List>(dim1[1])[1])[k])(0))/
        (Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[0])[0])(k)+Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[1])[0])(k));
      mu1_col(1,k)=(Rcpp::as<arma::vec>(Rcpp::as<List>(Rcpp::as<List>(dim1[0])[1])[k])(1)+Rcpp::as<arma::vec>(Rcpp::as<List>(Rcpp::as<List>(dim1[2])[1])[k])(0))/
        (Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[0])[0])(k)+Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[2])[0])(k));
      mu1_col(2,k)=(Rcpp::as<arma::vec>(Rcpp::as<List>(Rcpp::as<List>(dim1[1])[1])[k])(1)+Rcpp::as<arma::vec>(Rcpp::as<List>(Rcpp::as<List>(dim1[2])[1])[k])(1))/
        (Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[1])[0])(k)+Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[2])[0])(k));
    }
    
    for(k=0;k<linp;k++)
    {
      for(i=0;i<dim;i++)
      {
        m(i)=mu1_col(i,k);
      }
      mu1[k]=m;
    }
    for(k=0;k<linp;k++) pi1(k)=pi1(k)/(comb*n);
    List sigmafun(comb);
    List mmarg1(comb);
    for(i=0;i<comb;i++)
    {
      for(k=0;k<2;k++)
      {
        ro(k)=combn(i,k)-1;
      }
      smarg=margvbiv(sigma0,ro);
      mmarg=margmubiv(mu0,ro);
      mmarg1=margmubiv(mu1,ro);
      
      sigmafun[i]=sigmadimfunC(bin[i],Rcpp::as<vec>(Rcpp::as<List>(dim1[i])[2]),Rcpp::as<mat>(Rcpp::as<List>(dim1[i])[3]),Rcpp::as<mat>(Rcpp::as<List>(dim1[i])[4]),Rcpp::as<mat>(Rcpp::as<List>(dim1[i])[5]),Rcpp::as<mat>(Rcpp::as<List>(dim1[i])[6]),Rcpp::as<mat>(Rcpp::as<List>(dim1[i])[7]),Rcpp::as<mat>(Rcpp::as<List>(dim1[i])[8]),pi0, mmarg,mmarg1,smarg);
    } 
    
    
    arma::mat mm(dim,dim);
    arma::mat sigma1_col(dim,linp);
    for(k=0;k<linp;k++)
    {
      sigma1_col(0,k)=(Rcpp::as<arma::vec>(Rcpp::as<List>(sigmafun[0])[k])(0)+Rcpp::as<arma::vec>(Rcpp::as<List>(sigmafun[1])[k])(0))/
        (Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[0])[0])(k)+Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[1])[0])(k));
      sigma1_col(1,k)=(Rcpp::as<arma::vec>(Rcpp::as<List>(sigmafun[0])[k])(1)+Rcpp::as<arma::vec>(Rcpp::as<List>(sigmafun[2])[k])(0))/
        (Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[0])[0])(k)+Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[2])[0])(k));
      sigma1_col(2,k)=(Rcpp::as<arma::vec>(Rcpp::as<List>(sigmafun[1])[k])(1)+Rcpp::as<arma::vec>(Rcpp::as<List>(sigmafun[2])[k])(1))/
        (Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[1])[0])(k)+Rcpp::as<arma::vec>(Rcpp::as<List>(dim1[2])[0])(k));
    }
    for(k=0;k<linp;k++)
    {
      for(i=0;i<dim;i++)
      {
        mm(i,i)=sigma1_col(i,k);
      }
      sigma1[k]=mm;
    }
    j = j + 1;
    
    logli(j)=loglimargbivC(bin,pi1,mu1,sigma1);
    if(Rcpp::traits::is_nan<REALSXP>(logli(j))==true) 
    {return List::create(_["logli"]=-INFINITY);}
    flag=0;
    if ((abs((logli(j) - logli(j - 1))/logli(j - 1)) > eps) & 
        (j < it)) {
      flag = 1;
      pi0 = pi1;
      mu0 = mu1;
      sigma0 = sigma1;
    }
  }
  arma::vec logli1(j+1);
  for(i=0;i<=j;i++)
    logli1(i)=logli(i);
  return List::create(_["bin"]=bin,_["pi"]=pi1,_["mu"]=mu1,_["v"]=sigma1,_["logli"]=logli1);
}





// [[Rcpp::export]]
List binmixtbivC(arma::mat data, int cl, arma::vec R, int it, double eps, int nrep)
{
  List bin=buildbinbivmargC(data,R);
  int i;
  int dim=data.n_cols;
  arma::vec pi0(cl);
  List mu0(cl);
  List s0(cl);
  arma::vec mu01(dim);
  arma::mat s01(dim,dim);
  double s;
  int k;
  int j;
  Progress p(nrep, true);
  List e(nrep);
  arma::vec logli(nrep);
  arma::vec minimo(dim);
  arma::vec massimo(dim);
  arma::vec varianza(dim);
  for(i=0;i<dim;i++)
  {
    minimo(i)=arma::min(data.col(i));
    massimo(i)=arma::max(data.col(i));
    varianza(i)=arma::var(data.col(i));
  }
  for(i=0;i<nrep;i++)
  {
    pi0=arma::randu(cl);
    s=0;
    for(k=0;k<cl;k++)
      s=s+pi0(k);
    pi0=pi0/s;
    for(k=0;k<cl;k++)
    {
      for(j=0;j<dim;j++)
        mu01(j)=(arma::randu()*(massimo(j)))+minimo(j);
      mu0[k]=mu01;
    }
    for(k=0;k<cl;k++)
    {
      for(j=0;j<dim;j++)
        s01(j,j)=arma::randu()*varianza(j);
      s0[k]=s01;
    }
    List e1;
    e1=embinbivcompC(bin,pi0,mu0,s0,eps,it);
    e[i]=e1;
    logli(i)=Rcpp::as<arma::vec>((e1["logli"])).max();
    p.increment(); 
  }
  double logli1=-INFINITY;
  int indice=0;
  for(i=0;i<nrep;i++)
  {
    if(logli(i)>logli1)
    {
      logli1=logli(i);
      indice=i;
    }
  }
  return e[indice];
}


// [[Rcpp::export]]
List binmixtbivclassicC(arma::mat data, int cl, arma::vec R, int it, double eps, double eps1, int it1,int nrep)
{
  Rcpp::Rcout<<"Choosing the best initial guess..."<< std::endl;
  
  List migliore=binmixtbivC(data,cl,R,it,eps,nrep);
  Rcpp::Rcout<<"Estimation phase...";
  
  List e=embinbivcompC(migliore["bin"],migliore["pi"],
                       migliore["mu"],migliore["v"],eps1,it1);
  Rcpp::Rcout<<"Done."<< std::endl;
  return e;
}
