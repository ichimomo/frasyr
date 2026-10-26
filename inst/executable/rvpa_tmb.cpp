// ridge vpa

#include <TMB.hpp>
#include <iostream>

// https://github.com/ichimomo/frasyr/blob/a38f8bae1aa7f8f84b2c4e1c70ae7d656b7cb4ff/R/rvpa.r#L96-L114
// より移植（石岡・岸田の方法）
template<class Type>
Type ik_est(
    Type catch_n,
    Type next_n,
    Type mortality,
    int maxit)
{
  Type f = Type(1.0);

  for (int it = 0; it < maxit; ++it) {
    f = log(
      Type(1.0)
    + catch_n / next_n * exp(-mortality)
    * (f + mortality) * (Type(1.0) - exp(-f))
    / (f * (Type(1.0) - exp(-f - mortality)))
    );
  }

  return f;
}

// https://github.com/ichimomo/frasyr/blob/a38f8bae1aa7f8f84b2c4e1c70ae7d656b7cb4ff/R/rvpa.r#L116-L134
// より移植（平松の方法）
template<class Type>
Type hira_est(
    Type catch_prev,     // 当年の最高齢−1の漁獲尾数
    Type catch_plus,     // 当年の最高齢の漁獲尾数
    Type next_n,         // 翌年の最高齢の資源尾数
    Type m_prev,         // 当年の最高齢−1のM
    Type m_plus,         // 当年の最高齢のM
    Type alpha,
    int maxit)
{
  Type f = Type(1.0);

  for (int it = 0; it < maxit; ++it) {
    Type factor =
      (Type(1.0) - exp(-f)) * exp(-m_prev) / (next_n * f);

    Type plus_term =
      catch_plus * (alpha * f + m_plus)
      / (alpha * (Type(1.0) - exp(-alpha * f - m_plus)))
      * exp((Type(1.0) - alpha) * f);

    Type prev_term =
      catch_prev * (f + m_prev)
        / (Type(1.0) - exp(-f - m_prev));

    f = log(Type(1.0) + factor * (plus_term + prev_term));
  }

  return f;
}

template<class Type>
Type lower_bound_catch(Type catch_n, Type min_caa)
{
  return CppAD::CondExpGt(
    catch_n, min_caa, catch_n, min_caa
  );
}

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
  DATA_IVECTOR(Ab_type);  // Ab_type = 1: SSB, Ab_type = 2: number at age, Ab_type = 3: biomass at age, AD_type = 4: Bs(biomass x selectivity)
  DATA_INTEGER(sel_def); //selectivity definition (which F or age should be set at 1)
  // 0:max, 1: mean, 2: maxage
  DATA_IVECTOR(Ab_type_age);
  DATA_IVECTOR(Ab_type_max_age);
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
  DATA_INTEGER(Pope); //0: Baranov. 1: Pope
  DATA_INTEGER(maxit); // maximum iteration number for Baranov
  DATA_INTEGER(plus_group); //0: no plus group, 1: with plus group
  DATA_SCALAR(p_pope); // timing of catch within the (fishing) year
  DATA_SCALAR(min_caa); // replace zero caa by this value

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

  Type sum_log_ssb, sum_log_N, sum_log_B, sum_log_Bs;

  matrix<Type> F(Y,A);
  matrix<Type> N(Y,A);
  matrix<Type> B(Y,A);
  matrix<Type> Z(Y,A);
  matrix<Type> S(Y,A);
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
      if (Pope==1) {
        N(Y-1-Last_Catch_Zero,i)=CATCH(Y-1-Last_Catch_Zero,i)*exp(M(Y-1-Last_Catch_Zero,i)*p_pope)/(Type(1.0)-exp(-F(Y-1-Last_Catch_Zero,i)));
      } else { // Baranov
          N(Y-1-Last_Catch_Zero,i)=CATCH(Y-1-Last_Catch_Zero,i)*(1+M(Y-1-Last_Catch_Zero,i)/F(Y-1-Last_Catch_Zero,i))/(Type(1.0)-exp(-F(Y-1-Last_Catch_Zero,i)-M(Y-1-Last_Catch_Zero,i)));
        }
    }

  for (int t=(0+Last_Catch_Zero);t<Y-1;t++){
    for (int i=0;i<A-2;i++){
      if (Pope==1) {
        N(Y-2-t,i)=N(Y-1-t,i+1)*exp(M(Y-2-t,i))+CATCH(Y-2-t,i)*exp(M(Y-2-t,i)*p_pope);
        F(Y-2-t,i)=-log(1-CATCH(Y-2-t,i)*exp(M(Y-2-t,i)*p_pope)/N(Y-2-t,i));
      } else {
        F(Y-2-t,i) = ik_est(
          lower_bound_catch(CATCH(Y-2-t,i),min_caa),
          N(Y-1-t,i+1), M(Y-2-t,i), maxit);
        N(Y-2-t,i)=CATCH(Y-2-t,i)*(Type(1.0)+M(Y-2-t,i)/F(Y-2-t,i))/(Type(1.0)-exp(-F(Y-2-t,i)-M(Y-2-t,i)));
      }
    }

    Type c_prev = CATCH(Y-2-t,A-2);
    Type c_plus = CATCH(Y-2-t,A-1);
    Type denom = alpha*c_prev + c_plus;

    if (Pope==1) {
      if (plus_group==1) {
        N(Y-2-t,A-2) =
          c_prev*exp(M(Y-2-t,A-2)*p_pope)
        + N(Y-1-t,A-1)*exp(M(Y-2-t,A-2))
        *alpha*c_prev/denom;

        N(Y-2-t,A-1) =
        c_plus*exp(M(Y-2-t,A-1)*p_pope)
          + N(Y-1-t,A-1)*exp(M(Y-2-t,A-1))
          *c_plus/denom;
      } else {
        N(Y-2-t,A-2) =
          c_prev*exp(M(Y-2-t,A-2)*p_pope)
        + N(Y-1-t,A-1)*exp(M(Y-2-t,A-2));

        N(Y-2-t,A-1) =
        N(Y-2-t,A-2)*c_plus/c_prev
        *exp((M(Y-2-t,A-1)-M(Y-2-t,A-2))*p_pope);
      }
      F(Y-2-t,A-2) = -log(1-CATCH(Y-2-t,A-2)*exp(M(Y-2-t,A-2)*p_pope)/N(Y-2-t,A-2));
      F(Y-2-t,A-1) = alpha*F(Y-2-t,A-2);
    } else {
      if (plus_group==1) {
        F(Y-2-t,A-2) = hira_est(
          lower_bound_catch(CATCH(Y-2-t,A-2),min_caa),
          lower_bound_catch(CATCH(Y-2-t,A-1),min_caa),
          N(Y-1-t,A-1),         // 翌年の最高齢の資源尾数
          M(Y-2-t,A-2),         // 当年の最高齢−1のM
          M(Y-2-t,A-1),         // 当年の最高齢のM
          alpha,
          maxit);
      } else {
        F(Y-2-t,A-2) = ik_est(lower_bound_catch(CATCH(Y-2-t,A-2),min_caa),
          N(Y-1-t,A-1), M(Y-2-t,A-2), maxit);
      }
      F(Y-2-t,A-1) = alpha*F(Y-2-t,A-2);
      // 最高齢−1と最高齢のNをBaranovの式で計算
      for (int i = A-2; i < A; ++i) {
        N(Y-2-t,i) =
          CATCH(Y-2-t,i)
        * (Type(1.0) + M(Y-2-t,i)/F(Y-2-t,i))
        / (Type(1.0) - exp(-F(Y-2-t,i)-M(Y-2-t,i)));
      }
    }
  }

  if (Last_Catch_Zero==1){
    N(Y-1,0)=0.001;
    for (int i=0;i<A-2;i++){
      N(Y-1,i+1)=N(Y-2,i)*exp(-F(Y-2,i)-M(Y-2,i));
    }

    N(Y-1,A-1) =
      N(Y-2,A-2)*exp(-F(Y-2,A-2)-M(Y-2,A-2));

    if (plus_group == 1) {
      N(Y-1,A-1) +=
        N(Y-2,A-1)*exp(-F(Y-2,A-1)-M(Y-2,A-1));
    }

    // N(Y-1,A-1)=N(Y-2,A-2)*exp(-F(Y-2,A-2)-M(Y-2,A-2))+N(Y-2,A-1)*exp(-F(Y-2,A-1)-M(Y-2,A-1));
  }

  B = N.array()*WEI.array();

  for (int y=0;y<Y;y++){
    for (int i=0;i<A;i++){
      SSB(y) += B(y,i)*MAT(y,i);
    }
    logSSB(y) = log(SSB(y));
  }

  vector<Type> B_total(Y);
  B_total.fill(0.0);
  for (int y=0;y<Y;y++){
    for (int i=0;i<A;i++){
      B_total(y) += B(y,i);
    }
  }

  vector<Type> F_mean(Y);
  F_mean.fill(0.0);
  for (int y=0;y<Y;y++){
    for (int i=0;i<A;i++){
      F_mean(y) += F(y,i);
    }
    F_mean(y) /= A;
  }

  vector<Type> Catch_weight(Y);
  vector<Type> scale_U(Y);
  Catch_weight.fill(0.0);
  vector<Type> U(Y);
  U.fill(0.0);
  scale_U.fill(0.0);
  for (int y=0;y<Y;y++){
    for (int i=0;i<A;i++){
      Catch_weight(y) += CATCH(y,i)*WEI(y,i);
    }
    U(y) += Catch_weight(y)/B_total(y);
    scale_U(y) += U(y)/(1-U(y));
  }

  vector<Type> Sel1(Y);
  vector<Type> FY(A);
  for (int y=0;y<Y;y++){
    for (int i=0;i<A;i++) {
      FY(i)=F(y,i);
      }
    if (sel_def==0) { //max
      Sel1(y)=max(FY);
    }else{
      if (sel_def==1) { //mean
        Sel1(y)=F.row(y).sum();
        Sel1(y)/=A;
      }else{//maxage
        Sel1(y)=F(y,A-1);
        }
      }
    for (int i=0;i<A;i++){
      S(y,i)=F(y,i)/Sel1(y);
    }
    }

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
       for (int j=0;j<Y;j++){
         sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(SSB(j))),2);
       }
    }
    if (Ab_type(k)==2){
      sum_log_N=0.0;
      mean_log_abund(k)=0.0;
      sd_log_abund(k)=0.0;
      vector<Type> sum_N(Y);
      for (int j=0;j<Y;j++){
        sum_N(j)=0.0;
        for (int l=Ab_type_age(k);l<Ab_type_max_age(k);l++){
          sum_N(j) += N(j,l);
        }
        sum_log_cpue(k) += log(CPUE(j,k)+MISS(j,k));
        sum_log_N += NM(j,k)*log(sum_N(j));
       }
       if (b_fix(k)==0){
         for (int h=0;h<Y;h++){
           b(k) += NM(h,k)*((log(CPUE(h,k)+MISS(h,k))-sum_log_cpue(k)/nI(k))*(NM(h,k)*log(sum_N(h))-sum_log_N/nI(k)));
           denom += NM(h,k)*((NM(h,k)*log(sum_N(h))-sum_log_N/nI(k))*(NM(h,k)*log(sum_N(h))-sum_log_N/nI(k)));
         }
         b(k) /= denom;
       } else b(k) = b_fix(k);
       for (int h=0;h<Y;h++){
         q(k) += log(CPUE(h,k)+MISS(h,k))-b(k)*NM(h,k)*log(sum_N(h));
         mean_log_abund(k) += NM(h,k)*log(sum_N(h));
       }
       mean_log_abund(k) /= nI(k);
       for (int h=0;h<Y;h++){
         sd_log_abund(k) += pow(NM(h,k)*(log(sum_N(h))-mean_log_abund(k)),Type(2.0));
       }
       sd_log_abund(k) /= nI(k);
       sd_log_abund(k) = pow(sd_log_abund(k), Type(0.5));
       q(k) /= nI(k);
       q(k) = exp(q(k));
       for (int j=0;j<Y;j++){
         sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(sum_N(j))),2);
       }
    }
    if (Ab_type(k)==3){
      sum_log_B=0.0;
      mean_log_abund(k)=0.0;
      sd_log_abund(k)=0.0;
      vector<Type> sum_B(Y);
      for (int j=0;j<Y;j++){
         sum_B(j)=0.0;
         for (int l=Ab_type_age(k);l<Ab_type_max_age(k);l++){
           sum_B(j) += B(j,l);
         }
         sum_log_cpue(k) += log(CPUE(j,k)+MISS(j,k));
         sum_log_B += NM(j,k)*log(sum_B(j));
       }
       if (b_fix(k)==0){
         for (int h=0;h<Y;h++){
           b(k) += NM(h,k)*((log(CPUE(h,k)+MISS(h,k))-sum_log_cpue(k)/nI(k))*(NM(h,k)*log(sum_B(h))-sum_log_B/nI(k)));
           denom += NM(h,k)*((NM(h,k)*log(sum_B(h))-sum_log_B/nI(k))*(NM(h,k)*log(sum_B(h))-sum_log_B/nI(k)));
         }
         b(k) /= denom;
       } else b(k)=b_fix(k);
       for (int h=0;h<Y;h++){
         q(k) += log(CPUE(h,k)+MISS(h,k))-b(k)*NM(h,k)*log(sum_B(h));
         mean_log_abund(k) += NM(h,k)*log(sum_B(h));
       }
       mean_log_abund(k) /= nI(k);
       for (int h=0;h<Y;h++){
         sd_log_abund(k) += pow(NM(h,k)*(log(sum_B(h))-mean_log_abund(k)),Type(2.0));
       }
       sd_log_abund(k) /= nI(k);
       sd_log_abund(k) = pow(sd_log_abund(k), Type(0.5));
       q(k) /= nI(k);
       q(k) = exp(q(k));
       for (int j=0;j<Y;j++){
         sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(sum_B(j))),2);
       }
    }
    if (Ab_type(k)==4){
      sum_log_Bs=0.0;
      mean_log_abund(k)=0.0;
      sd_log_abund(k)=0.0;
      vector<Type> sum_Bs(Y);
      for (int j=0;j<Y;j++){
        sum_Bs(j)=0.0;
        for (int l=Ab_type_age(k);l<Ab_type_max_age(k);l++){
          sum_Bs(j) += B(j,l)*S(j,l);
        }
        sum_log_cpue(k) += log(CPUE(j,k)+MISS(j,k));
        sum_log_Bs += NM(j,k)*log(sum_Bs(j));
      }
      if (b_fix(k)==0){
        for (int h=0;h<Y;h++){
          b(k) += NM(h,k)*((log(CPUE(h,k)+MISS(h,k))-sum_log_cpue(k)/nI(k))*(NM(h,k)*log(sum_Bs(h))-sum_log_Bs/nI(k)));
          denom += NM(h,k)*((NM(h,k)*log(sum_Bs(h))-sum_log_Bs/nI(k))*(NM(h,k)*log(sum_Bs(h))-sum_log_Bs/nI(k)));
        }
        b(k) /= denom;
      } else b(k)=b_fix(k);
      for (int h=0;h<Y;h++){
        q(k) += log(CPUE(h,k)+MISS(h,k))-b(k)*NM(h,k)*log(sum_Bs(h));
        mean_log_abund(k) += NM(h,k)*log(sum_Bs(h));
      }
      mean_log_abund(k) /= nI(k);
      for (int h=0;h<Y;h++){
        sd_log_abund(k) += pow(NM(h,k)*(log(sum_Bs(h))-mean_log_abund(k)),Type(2.0));
      }
      sd_log_abund(k) /= nI(k);
      sd_log_abund(k) = pow(sd_log_abund(k), Type(0.5));
      q(k) /= nI(k);
      q(k) = exp(q(k));
      for (int j=0;j<Y;j++){
        sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(sum_Bs(j))),2);
      }
      }
  }

  Type f=0;

  // for (int k=0;k<K;k++){
  //   if(Ab_type(k)==1){
  //     for (int j=0;j<Y;j++){
  //       sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(SSB(j))),2);
  //     }
  //   }
  //   if(Ab_type(k)==2){
  //     for (int j=0;j<Y;j++){
  //       sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(sum_N(j))),2);
  //     }
  //   }
  //   if(Ab_type(k)==3){
  //     for (int j=0;j<Y;j++){
  //       sigma2(k) += pow(log(CPUE(j,k)+MISS(j,k))-NM(j,k)*(log(q(k))+b(k)*log(sum_B(j))),2);
  //     }
  //   }
  // }
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
      // f += 0.5*nI(k)*(1+log(2*M_PI))+0.5*nI(k)*(log(sigma2(k))-log(nI(k)));
      f += Type(0.5)*nI(k)*log(Type(2.0)*M_PI*sigma(k))+Type(0.5)*sigma2(k)/sigma(k);
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

  ADREPORT(N);
  ADREPORT(F);
  ADREPORT(SSB);
  ADREPORT(B_total);
  ADREPORT(F_mean);
  ADREPORT(U);
  ADREPORT(scale_U);

  return f;
}
