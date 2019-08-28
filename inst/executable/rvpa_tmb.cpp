// ridge vpa

#include <TMB.hpp>
#include <iostream>

/* Parameter transform */
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data Section //
  DATA_INTEGER(Est);  // Est = 0: LS, Est = 1: ML
  DATA_VECTOR(b_fix);  //
  DATA_SCALAR(alpha);
  DATA_SCALAR(lambda);
  DATA_SCALAR(beta);  
  DATA_IVECTOR(Ab_type);  // Ab_type = 1: SSB, Ab_type = 2: number at age, Ab_type = 3: biomass at age
  DATA_IVECTOR(Ab_type_age);
  DATA_VECTOR(w);
  DATA_VECTOR(af);
  DATA_MATRIX(CATCH);
  DATA_MATRIX(WEI);
  DATA_MATRIX(MAT);
  DATA_MATRIX(M);
  DATA_MATRIX(CPUE);
  DATA_MATRIX(MISS);
  DATA_INTEGER(Last_Catch_Zero);
  DATA_IVECTOR(sigma_constraint);
  DATA_SCALAR(eta);
  DATA_IVECTOR(eta_age);

  // Parameter Section //
  PARAMETER_VECTOR(log_F);
  // PARAMETER_VECTOR(log_b);

  int Y=CATCH.rows();
  int A=CATCH.cols();
  int K=CPUE.cols();

  vector<Type> q(K);
  vector<Type> b(K);
  vector<Type> nI(K);
  vector<Type> sigma2(K);
  vector<Type> sum_log_cpue(K);
  vector<Type> mean_log_abund(K);
  vector<Type> sd_log_abund(K);
  vector<Type> sigma(K);
  vector<Type> nI2(K);
  
  q.fill(0.0);
  b.fill(0.0);
  nI.fill(0.0);
  sigma2.fill(0.0);
  sum_log_cpue.fill(0.0);
  mean_log_abund.fill(0.0);
  sd_log_abund.fill(0.0);
  sigma.fill(0.0);
  nI2.fill(0.0);
    
  Type sum_log_ssb, sum_log_N, sum_log_B;
  
  matrix<Type> F(Y,A);
  matrix<Type> N(Y,A);
  matrix<Type> B(Y,A);
  matrix<Type> Z(Y,A);
  F.fill(0.0);
  N.fill(0.0);
      
  matrix<Type> NM(Y,K);
  NM.fill(1.0);

  vector<Type> SSB(Y);
  vector<Type> logSSB(Y);
  SSB.fill(0.0);
            
  //
  
  NM = NM-MISS;

  for (int i=0;i<A;i++){
      if (i < A-1){F(Y-1-Last_Catch_Zero,i) = exp(log_F(i));}else{F(Y-1-Last_Catch_Zero,i)=alpha*F(Y-1-Last_Catch_Zero,A-2);}
      N(Y-1-Last_Catch_Zero,i)=CATCH(Y-1-Last_Catch_Zero,i)*exp(M(Y-1-Last_Catch_Zero,i)/2)/(1-exp(-F(Y-1-Last_Catch_Zero,i)));
  }

  for (int t=(0+Last_Catch_Zero);t<Y-1;t++){
    for (int i=0;i<A-2;i++){
      N(Y-2-t,i)=N(Y-1-t,i+1)*exp(M(Y-2-t,i))+CATCH(Y-2-t,i)*exp(M(Y-2-t,i)/2);
      F(Y-2-t,i)=-log(1-CATCH(Y-2-t,i)*exp(M(Y-2-t,i)/2)/N(Y-2-t,i));
    }
    N(Y-2-t,A-2) = CATCH(Y-2-t,A-2)*exp(M(Y-2-t,A-2)/2)+N(Y-1-t,A-1)*exp(M(Y-2-t,A-2))*CATCH(Y-2-t,A-2)/(CATCH(Y-2-t,A-2)+CATCH(Y-2-t,A-1));
    N(Y-2-t,A-1) = CATCH(Y-2-t,A-1)*exp(M(Y-2-t,A-1)/2)+N(Y-1-t,A-1)*exp(M(Y-2-t,A-1))*CATCH(Y-2-t,A-1)/(CATCH(Y-2-t,A-2)+CATCH(Y-2-t,A-1));
    F(Y-2-t,A-2) = -log(1-CATCH(Y-2-t,A-2)*exp(M(Y-2-t,A-2)/2)/N(Y-2-t,A-2));
    F(Y-2-t,A-1) = alpha*F(Y-2-t,A-2);
  }  
  
  if (Last_Catch_Zero==1){
    N(Y-1,0)=0.001;
    for (int i=0;i<A-2;i++){
      N(Y-1,i+1)=N(Y-2,i)*exp(-F(Y-2,i)-M(Y-2,i));
    }
    N(Y-1,A-1)=N(Y-2,A-2)*exp(-F(Y-2,A-2)-M(Y-2,A-2))+N(Y-2,A-1)*exp(-F(Y-2,A-1)-M(Y-2,A-1));
  }
  
  B = N.array()*WEI.array();
  
  for (int y=0;y<Y;y++){
    for (int i=0;i<A;i++){
      SSB(y) += B(y,i)*MAT(y,i);
    }
    logSSB(y) = log(SSB(y));
  }
  
  //
  
  Type denom;
  
  for (int k=0;k<K;k++){
    denom=0.0;
    nI(k) = NM.col(k).sum();
    if (Ab_type(k)==1){
       sum_log_ssb=0.0;
       mean_log_abund(k)=0.0;
       sd_log_abund(k)=0.0;
       for (int j=0;j<Y;j++){
           sum_log_cpue(k) += log(CPUE(j,k)+MISS(j,k));
           sum_log_ssb += NM(j,k)*logSSB(j);
       }
       if (b_fix(k)==0){
         for (int h=0;h<Y;h++){
           b(k) += NM(h,k)*((log(CPUE(h,k)+MISS(h,k))-sum_log_cpue(k)/nI(k))*(NM(h,k)*logSSB(h)-sum_log_ssb/nI(k)));
           denom += NM(h,k)*((NM(h,k)*logSSB(h)-sum_log_ssb/nI(k))*(NM(h,k)*logSSB(h)-sum_log_ssb/nI(k)));
         }
         b(k) /= denom;          
       } else b(k)=b_fix(k);
       for (int h=0;h<Y;h++){
         q(k) += log(CPUE(h,k)+MISS(h,k))-b(k)*NM(h,k)*logSSB(h);
         mean_log_abund(k) += NM(h,k)*logSSB(h);
       }
       mean_log_abund(k) /= nI(k);
       for (int h=0;h<Y;h++){
         sd_log_abund(k) += pow(NM(h,k)*(log(SSB(h))-mean_log_abund(k)),Type(2.0));
       }
       sd_log_abund(k) /= nI(k);
       sd_log_abund(k) = pow(sd_log_abund(k), Type(0.5));
       q(k) /= nI(k);
       q(k) = exp(q(k));
    }
    if (Ab_type(k)==2){
      sum_log_N=0.0;
      mean_log_abund(k)=0.0;
      sd_log_abund(k)=0.0;
      for (int j=0;j<Y;j++){
           sum_log_cpue(k) += log(CPUE(j,k)+MISS(j,k));
           sum_log_N += NM(j,k)*log(N(j,Ab_type_age(k)));
       }
       if (b_fix(k)==0){
         for (int h=0;h<Y;h++){
           b(k) += NM(h,k)*((log(CPUE(h,k)+MISS(h,k))-sum_log_cpue(k)/nI(k))*(NM(h,k)*log(N(h,Ab_type_age(k)))-sum_log_N/nI(k)));
           denom += NM(h,k)*((NM(h,k)*log(N(h,Ab_type_age(k)))-sum_log_N/nI(k))*(NM(h,k)*log(N(h,Ab_type_age(k)))-sum_log_N/nI(k)));
         }
         b(k) /= denom;
       } else b(k) = b_fix(k);
       for (int h=0;h<Y;h++){
         q(k) += log(CPUE(h,k)+MISS(h,k))-b(k)*NM(h,k)*log(N(h,Ab_type_age(k)));
         mean_log_abund(k) += NM(h,k)*log(N(h,Ab_type_age(k)));
       }
       mean_log_abund(k) /= nI(k);
       for (int h=0;h<Y;h++){
         sd_log_abund(k) += pow(NM(h,k)*(log(N(h,Ab_type_age(k)))-mean_log_abund(k)),Type(2.0));
       }
       sd_log_abund(k) /= nI(k);
       sd_log_abund(k) = pow(sd_log_abund(k), Type(0.5));
       q(k) /= nI(k);
       q(k) = exp(q(k));
    }
    if (Ab_type(k)==3){
      sum_log_B=0.0;
      mean_log_abund(k)=0.0;
      sd_log_abund(k)=0.0;
       for (int j=0;j<Y;j++){
           sum_log_cpue(k) += log(CPUE(j,k)+MISS(j,k));
           sum_log_B += NM(j,k)*log(B(j,Ab_type_age(k)));
       }
       if (b_fix(k)==0){
         for (int h=0;h<Y;h++){
           b(k) += NM(h,k)*((log(CPUE(h,k)+MISS(h,k))-sum_log_cpue(k)/nI(k))*(NM(h,k)*log(B(h,Ab_type_age(k)))-sum_log_B/nI(k)));
           denom += NM(h,k)*((NM(h,k)*log(B(h,Ab_type_age(k)))-sum_log_B/nI(k))*(NM(h,k)*log(B(h,Ab_type_age(k)))-sum_log_B/nI(k)));
         }
         b(k) /= denom;  
       } else b(k)=b_fix(k);
       for (int h=0;h<Y;h++){
         q(k) += log(CPUE(h,k)+MISS(h,k))-b(k)*NM(h,k)*log(B(h,Ab_type_age(k)));
         mean_log_abund(k) += NM(h,k)*log(B(h,Ab_type_age(k)));
       }
       mean_log_abund(k) /= nI(k);
       for (int h=0;h<Y;h++){
         sd_log_abund(k) += pow(NM(h,k)*(log(B(h,Ab_type_age(k)))-mean_log_abund(k)),Type(2.0));
       }
       sd_log_abund(k) /= nI(k);
       sd_log_abund(k) = pow(sd_log_abund(k), Type(0.5));
       q(k) /= nI(k);
       q(k) = exp(q(k));
    }
  }  
  
  Type f=0;
  
  for (int k=0;k<K;k++){
    if(Ab_type(k)==1){
      for (int j=0;j<Y;j++){
        sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(SSB(j))),2);
      }
    }
    if(Ab_type(k)==2){
      for (int j=0;j<Y;j++){
        sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(N(j,Ab_type_age(k)))),2);
      }
    }
    if(Ab_type(k)==3){
      for (int j=0;j<Y;j++){
        sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(B(j,Ab_type_age(k)))),2);
      }
    }
  }
  // sigma2 = Residual of sum of square
  
  if (Est==0){
    for (int k=0;k<K;k++){
      f += w(k)*sigma2(k);
    }
  }
  if (Est==1){
    for (int k=0;k<K;k++){
      for (int j=0;j<K;j++) {
        if (sigma_constraint(j)==sigma_constraint(k)) {
          sigma(k) += sigma2(j);
          nI2(k) += nI(j);
        }
      }
      sigma(k) /= nI2(k); // sigma^2
      // sigma(k) = pow(sigma(k), Type(0.5));
    }
    for (int k=0;k<K;k++){
      // f += 0.5*nI(k)*(1+log(2*PI))+0.5*nI(k)*(log(sigma2(k))-log(nI(k)));
      f += Type(0.5)*nI(k)*log(Type(2.0)*PI*sigma(k))+Type(0.5)*sigma2(k)/sigma(k);
    }
  }
  f *= (1-lambda);
  if (eta==-1) {
    for (int k=0;k<log_F.size();k++) {
      f += lambda*pow(exp(log_F(k)),beta);
    }
  } else {
    for (int k=0;k<log_F.size();k++){
      if (eta_age(k) == 0) {
        f += lambda*eta*pow(exp(log_F(k)),beta);
      } else {
        f += lambda*(1-eta)*pow(exp(log_F(k)),beta);
      }
    }
  }

  // 
  return f;
}
