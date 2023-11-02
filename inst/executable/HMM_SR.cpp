
#include <TMB.hpp>
#include <iostream>

// SR function (prediction of log(RPS))
template<class Type>
Type srf(Type alpha, Type beta, Type gamma, Type x, int SRcode){
  Type res=0.0;
  if (SRcode==1) { //Ricker
    res += alpha-beta*x;
  } else {
    if (SRcode==2) { //BH
      res += alpha-log(1+beta*x);
    } else {
      if (SRcode==3) { //HS
        res += CppAD::CondExpLt(beta,x,alpha+log(beta/x),alpha);
      } else {
        if (SRcode==4) { //Mesnil
          res += x+pow(pow(x,Type(2.0))+pow(gamma,Type(2.0))/Type(4.0),Type(0.5))-pow(pow(x-beta,Type(2.0))+pow(gamma,Type(2.0))/Type(4.0),Type(0.5));
          res *= Type(0.5);
          res = log(res); //transformation to log(R)
          res += alpha-log(x); //transformation to log(RPS)
        } else {
          error("SR model code not recognized");
        }
      }
    }
  }
  return res;
}

template<class Type>
vector<Type> segment_1(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type> pi1,vector<Type> alpha,vector<Type> beta,vector<Type> sigma,int t, Type gamma, int SRcode){
  int k_regime = beta.size();
  Type small = pow(10,-300);

  // t = 0
      vector<Type> sr = log(pi1 + small); //log-likelihood of observation
  for(int j = 0;j < k_regime;++j){
     Type f_now = srf(alpha(j), beta(j),gamma,st(0),SRcode);
     // Type f_now = alpha(j) - log( 1 + beta(j)*st(0) );
     sr(j) += dnorm(yt(0), f_now, sigma(j),true);
  }

  // t >= 1
  for(int i = 1;i <= t;++i){

     vector<Type> sr_new = sr;
     for(int j = 0;j < k_regime;++j){
         sr_new(j) = sr(0) +qij(0,j);
         for(int jj = 1;jj < k_regime;++jj){
            Type temp = sr(jj) +qij(jj,j);
            sr_new(j) =  logspace_add(sr_new(j),temp); //log(exp(sr_new)+exp(temp))
         }
     }

     sr = sr_new;

     for(int j = 0;j < k_regime;++j){
        Type f_now = srf(alpha(j), beta(j),gamma,st(i),SRcode);
        // Type f_now = alpha(j) - log( 1 + beta(j)*st(i) );
        sr(j) += dnorm(yt(i), f_now, sigma(j),true); //log-likelihood of regime j
     }

  }

  return sr;
}

template<class Type>
Type segment_2(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type> pi1,vector<Type> alpha,vector<Type> beta,vector<Type> sigma,int rt,int t, Type gamma, int SRcode){
  int k_regime = beta.size();
  int n = yt.size();

  vector<Type> sr = qij.row(rt); // rt: regime ID, log transition probability
  for(int j = 0;j < k_regime;++j){
     Type f_now = srf(alpha(j), beta(j),gamma,st(t+1),SRcode);
     // Type f_now = alpha(j) - log( 1 + beta(j)*st(t+1) ); //t+1
     sr(j) += dnorm(yt(t+1), f_now, sigma(j),true);
  }

  for(int i = t+2;i < n;++i){ //t+2~n-1

     vector<Type> sr_new = sr;
     for(int j = 0;j < k_regime;++j){
         sr_new(j) = sr(0) +qij(0,j);
         for(int jj = 1;jj < k_regime;++jj){
            Type temp = sr(jj) +qij(jj,j);
            sr_new(j) =  logspace_add(sr_new(j),temp);
         }
     }

     sr = sr_new;

     for(int j = 0;j < k_regime;++j){
        Type f_now = srf(alpha(j), beta(j),gamma,st(i),SRcode);
        // Type f_now = alpha(j) - log( 1+ beta(j)*st(i) );
        sr(j) += dnorm(yt(i), f_now, sigma(j),true);
     }

  }

  Type seg2 = sr(0);
  for(int j = 1;j < k_regime;++j){
     seg2 = logspace_add(seg2,sr(j));
  }

  return seg2;
}

//main //////////////////////////////////////////////

template<class Type>
  Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(yt);  //log(R/SSB)
  DATA_VECTOR(st); // SSB
  DATA_SCALAR(alpha_u); //upper value of alpha
  DATA_SCALAR(alpha_l); //lower value of alpha
  DATA_SCALAR(beta_u); //upper value of beta
  DATA_SCALAR(beta_l); //lower value of beta
  DATA_SCALAR(sigma_u); //upper value of sigma
  DATA_INTEGER(SRcode);//1: Ricker, 2: BH, 3: Hockey-stock, 4: Mesnil
  DATA_SCALAR(gamma);

  PARAMETER_VECTOR(lalpha);
  PARAMETER_VECTOR(lbeta);
  PARAMETER_VECTOR(lsigma);
  PARAMETER_VECTOR(pi1_tran); // regime probability of initial states (length = number of regimes)
  PARAMETER_MATRIX(qij_tran); // transition probability (nrow=number of regimes, ncol = number of regimes -1)


  int k_regime = lbeta.size();   // Number of regimes

//  vector<Type> alpha = alpha_tr;
//  for(int i = 1;i < k_regime;++i){
//   alpha(i) = alpha(i-1) + exp(alpha_tr(i));
//  }

  vector<Type> beta = beta_l+(beta_u-beta_l)/(1+exp(-lbeta));// when lbeta is negative infinity, beta=beta_l; when lbeta is positive infinity, beta=beta_u

  vector<Type> alpha(k_regime);
  alpha(0) = (alpha_u-alpha_l)/(1+exp(-lalpha(0)))+alpha_l; // alpha transformation
  for(int i = 1;i < k_regime;++i){
    alpha(i) = alpha(i-1) + (alpha_u-alpha(i-1))/(1+exp(-lalpha(i)));  // alpha(i) > alpha(i-1)
  } // alpha(1) from alpha(0) to alpha_u

  vector<Type> sigma = sigma_u/(1+exp(-lsigma)); // sigmaR

  // initial probabilities
  vector<Type> pi1(k_regime);

  for(int i = 0;i < k_regime-1;++i){
    pi1(i) = exp(pi1_tran(i));
  }
  pi1(k_regime-1) = 1;
  pi1 = pi1/(pi1.sum());

  Type small = pow(10,-300);

  matrix<Type> qij(k_regime,k_regime);
  for(int i = 0;i < k_regime;++i){
    for(int j = 0;j < k_regime-1;++j){
      qij(i,j) = exp(qij_tran(i,j));
    }
    qij(i,k_regime-1) = 1; // fixed at 1 (tentative)
    vector<Type> qij_row = qij.row(i);
    Type row_sum = qij_row.sum();
    for(int j = 0;j < k_regime;++j){
      qij(i,j) = qij(i,j)/row_sum;
      qij(i,j) = log(qij(i,j)+small);
    } //transform to log probability
  }

  int n = yt.size();

  vector<Type> sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,n-1,gamma,SRcode);

  Type nll = sr(0);
  for(int j = 1;j < k_regime;++j){
     nll = logspace_add(nll,sr(j)); //log(exp(nll)+exp(sr(j)))
    // Tang et al. (2021) ICES JMSのequation(3)の最初のΣなので確率の足し算で合っている
    // Supplementary Materials A
  }

  nll = -nll;

// predict r_t ///////////////////////////////////

  matrix<Type> r_pred(k_regime,n);

  for(int i = 0;i < n-1;++i){
     sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,i,gamma,SRcode);
     for(int j = 0;j < k_regime;++j){
         Type tempt = sr(j) + segment_2(yt, st, qij,pi1,alpha,beta,sigma,j,i,gamma,SRcode) + nll;
         r_pred(j,i) = exp(tempt);
     }
  }

  sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,n-1,gamma,SRcode); //for the last year
  for(int j = 0;j < k_regime;++j){
      Type tempt = sr(j) + nll;
      r_pred(j,n-1) = exp(tempt);
  }

  qij = exp(qij.array());

  REPORT(beta);
  REPORT(alpha);
  REPORT(sigma);
  REPORT(pi1);
  REPORT(qij);
  REPORT(r_pred);

  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(sigma);
  ADREPORT(pi1);
  ADREPORT(qij);

  return nll;
  }


