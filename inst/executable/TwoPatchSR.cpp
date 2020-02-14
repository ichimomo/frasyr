// stock-recruit of two-patch models

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data section
  DATA_VECTOR(rec1);
  DATA_VECTOR(ssb1);
  DATA_VECTOR(rec2);
  DATA_VECTOR(ssb2);
  DATA_IVECTOR(iy); // ID representing years for the patch with immigration (2)
  DATA_IVECTOR(SR); //0:HS,1:BH,2:Ricker,3:RPS
  DATA_INTEGER(mig_type); //0: no migration, 1: constant migration, 2: white noise, 3: AR(1), 4: Random walk
  DATA_INTEGER(one_patch); //0: two patch (two SR), 1: one patch (one SR)
  // DATA_INTEGER(shared_error); //0: two patch (two SR), 1: one patch (one SR)
  DATA_SCALAR(sigma_omega_ratio);
  DATA_VECTOR(mig_occurrence); // year ID representing migration occurrence (1 or 0)

  // Parameter section
  PARAMETER_VECTOR(rec_loga);
  PARAMETER_VECTOR(rec_logb);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(trans_rho);
  PARAMETER_VECTOR(logit_d); // random effects 
  PARAMETER(intercept_d); // logarithmic mean immigrants 
  PARAMETER(trans_phi); //autoregressive coefficient
  PARAMETER(log_omega); // sigma for autocorrelated immigrants
  
  Type nll = 0.0; // negative log-likelihood
  
  vector<Type> rec_a = exp(rec_loga);
  vector<Type> sigma = exp(log_sigma);
  vector<Type> rec_b(rec_logb.size());

  Type phi = (exp(trans_phi)-Type(1.0))/(exp(trans_phi)+Type(1.0));
  vector<Type> rho = (exp(trans_rho)-Type(1.0))/(exp(trans_rho)+Type(1.0));
  
  // Migration Process

  vector<Type> d(logit_d.size());
  d.fill(0.0); // no migration
  vector<Type> tau(logit_d.size()); // process error in migration rate
  tau.fill(0.0);
  Type omega = exp(log_omega);
  if (sigma_omega_ratio > 0.0) omega = sigma_omega_ratio*sigma(0);
  
  if (mig_type==1) { //contant migration
    for (int i=0;i<d.size(); i++) {
      d(i) += Type(1.0)/(Type(1.0)+exp(-intercept_d));
      d(i) *= mig_occurrence(i);
    }
  } else {
    if (mig_type==2) { //white noise
      for (int i=0;i<d.size(); i++) {
        tau(i) += logit_d(i)-intercept_d;
        nll -= dnorm(tau(i),Type(0.0),omega,true);
        d(i) += Type(1.0)/(Type(1.0)+exp(-logit_d(i)));
        d(i) *= mig_occurrence(i);
      }
    } else {
      if (mig_type==3) { //AR(1)
        for (int i=0;i<d.size(); i++) {
          if (i==0) {
            tau(i) += logit_d(i)-intercept_d;
            nll -= dnorm(tau(i),Type(0.0),pow(pow(omega,Type(2.0))/(Type(1.0)-pow(phi,Type(2.0))),Type(0.5)),true);
          } else {
            tau(i) += (logit_d(i)-intercept_d)-phi*(logit_d(i-1)-intercept_d);
            nll -= dnorm(tau(i),Type(0.0),omega,true);
          }
          d(i) += Type(1.0)/(Type(1.0)+exp(-logit_d(i)));
          d(i) *= mig_occurrence(i);
        }
      } else {
          if (mig_type==4) { //Random Walk
            phi = Type(1.0);
            for (int i=1;i<d.size(); i++) { //start at i=1
              tau(i) += logit_d(i) - logit_d(i-1);
              nll -= dnorm(logit_d(i),logit_d(i-1),omega,true);
          }
            d += Type(1.0)/(Type(1.0)+exp(-logit_d));
            d *= mig_occurrence;
          } else {
            if (mig_type != 0) error("mig_type not recognized");
          }
        }
      }  
    }

  // defining stock-recruitment parameters
  if (one_patch==0) {
    for (int i=0;i<2;i++) {
      if (SR(i)==0) {
        if (i==0) {
          rec_b(i)  = min(ssb1)+(max(ssb1)-min(ssb1))/(1+exp(-rec_logb(i)));
        }
        if (i==1) {
          rec_b(i)  = min(ssb2)+(max(ssb2)-min(ssb2))/(1+exp(-rec_logb(i)));
        }
      } else {
        rec_b(i)=exp(rec_logb(i));
      }
    }
  } else {
    vector<Type> ssb_both = ssb1+ssb2;
    vector<Type> rec_both = rec1+rec2;
    for (int i=0;i<1;i++) {
      if (SR(i)==0) {
        rec_b(i)  = min(ssb_both)+(max(ssb_both)-min(ssb_both))/(1+exp(-rec_logb(i)));
      } else {
        rec_b(i)=exp(rec_logb(i));
      }
    }
  }

  // patch dynamics with emmigration (1)
  vector<Type> pred_repro1(rec1.size());
  pred_repro1.fill(0.0);
  vector<Type> delta1(rec1.size());
  delta1.fill(0.0);
  vector<Type> repro1(rec1.size());
  repro1.fill(0.0);
  vector<Type> migrant(rec1.size());
  migrant.fill(0.0);
  vector<Type> pred_repro2(rec2.size());
  pred_repro2.fill(0.0);
  vector<Type> delta2(rec2.size());
  delta2.fill(0.0);
  vector<Type> repro2(rec2.size());
  repro2.fill(0.0);
  
  if (one_patch==0) { // two patches
    for (int i=0;i<rec1.size();i++) {
      if (SR(0)==0) {
        pred_repro1(i) += CppAD::CondExpLt(rec_b(0),ssb1(i),rec_a(0)*rec_b(0),rec_a(0)*ssb1(i));
      }
      // Beverton-Holt
      if (SR(0)==1) {
        pred_repro1(i) += rec_a(0)*ssb1(i)/(Type(1.0)+rec_b(0)*ssb1(i));
      }
      // Ricker
      if (SR(0)==2) {
        pred_repro1(i) += rec_a(0)*ssb1(i)*exp(-rec_b(0)*ssb1(i));
      }
      // RPS = constant (no density-dependent effect)
      if (SR(0)==3) {
        pred_repro1(i) += rec_a(0)*ssb1(i);
      }
      delta1(i) += log(rec1(i))-log(1-d(i))-log(pred_repro1(i));
      if (i==0) {
        nll -= dnorm(delta1(i),Type(0.0),pow(pow(sigma(0),Type(2.0))/(Type(1.0)-pow(rho(0),Type(2.0))),Type(0.5)),true);
      } else {
        nll -= dnorm(delta1(i),rho(0)*delta1(i-1),sigma(0),true);
      }
      repro1(i) += pred_repro1(i)*exp(delta1(i));
      migrant(i) += d(i)*repro1(i);
    }
  } else { // one patch
    vector<Type> ssb_both = ssb1+ssb2;
    vector<Type> rec_both = rec1+rec2;
    for (int i=0;i<rec1.size();i++) {
      if (SR(0)==0) { // Hockey-stick
        pred_repro1(i) += CppAD::CondExpLt(rec_b(0),ssb_both(i),rec_a(0)*rec_b(0),rec_a(0)*ssb_both(i));
      } else {
        if (SR(0)==1) { // Beverton-Holt
          pred_repro1(i) += rec_a(0)*ssb_both(i)/(Type(1.0)+rec_b(0)*ssb_both(i));
        } else {
          // Ricker
          if (SR(0)==2) {
            pred_repro1(i) += rec_a(0)*ssb_both(i)*exp(-rec_b(0)*ssb_both(i));
          } else {
            // RPS = constant (no density-dependent effect)
            if (SR(0)==3) {
              pred_repro1(i) += rec_a(0)*ssb_both(i);
            } else {
              error("sR code not recognized");
            }
          }
        }
      }
      delta1(i) += log(rec_both(i))-log(pred_repro1(i));
      pred_repro2(i) += d(i)*pred_repro1(i);
      pred_repro1(i) *= (Type(1.0)-d(i));
      if (i==0) {
        nll -= dnorm(log(rec1(i)),log(pred_repro1(i)),pow(pow(sigma(0),Type(2.0))/(Type(1.0)-pow(rho(0),Type(2.0))),Type(0.5)),true);
        nll -= dnorm(log(rec2(i)),log(pred_repro2(i)),pow(pow(sigma(0),Type(2.0))/(Type(1.0)-pow(rho(0),Type(2.0))),Type(0.5)),true);
      } else {
        nll -= dnorm(log(rec1(i)),log(pred_repro1(i))+rho(0)*delta1(i-1),sigma(0),true);
        nll -= dnorm(log(rec2(i)),log(pred_repro2(i))+rho(0)*delta1(i-1),sigma(0),true);
      }
      repro1(i) += pred_repro1(i)*exp(delta1(i));
      repro2(i) += pred_repro2(i)*exp(delta2(i));
    }
  }

  // patch dynamics with immigration (2)
  if (one_patch==0) {
    for (int i=0;i<rec2.size();i++) {
      if (SR(1)==0) {
        pred_repro2(i) += CppAD::CondExpLt(rec_b(1),ssb2(i),rec_a(1)*rec_b(1),rec_a(1)*ssb2(i));
      }
      // Beverton-Holt
      if (SR(1)==1) {
        pred_repro2(i) += rec_a(1)*ssb2(i)/(Type(1.0)+rec_b(1)*ssb2(i));
      }
      // Ricker
      if (SR(1)==2) {
        pred_repro2(i) += rec_a(1)*ssb2(i)*exp(-rec_b(1)*ssb2(i));
      }
      // RPS = constant (no density-dependent effect)
      if (SR(1)==3) {
        pred_repro2(i) += rec_a(1)*ssb2(i);
      }
      // if (shared_error == 1) { // 移入と繁殖のプロセスエラーを共通に
      //   pred_repro2(i) += migrant(iy(i));
      //   repro2(i) += rec2(i);
      // } else {
        repro2(i) += CppAD::CondExpLt(rec2(i)-migrant(iy(i)),Type(0.000001),Type(0.000001),rec2(i)-migrant(iy(i)));
      // }
      delta2(i) += log(repro2(i))-log(pred_repro2(i));
      if (i==0) {
        nll -= dnorm(delta2(i),Type(0.0),pow(pow(sigma(1),Type(2.0))/(Type(1.0)-pow(rho(1),Type(2.0))),Type(0.5)),true);
      } else {
        nll -= dnorm(delta2(i),rho(1)*delta2(i-1),sigma(1),true);
      }
    }
  } 

  ADREPORT(phi);
  ADREPORT(rho);
  ADREPORT(d);
  ADREPORT(tau);
  ADREPORT(omega);
  ADREPORT(rec_a);
  ADREPORT(rec_b);
  ADREPORT(sigma);
  ADREPORT(pred_repro1);
  ADREPORT(delta1);
  ADREPORT(repro1);
  ADREPORT(migrant);
  ADREPORT(pred_repro2);
  ADREPORT(delta2);
  ADREPORT(repro2);
  // ADREPORT(rec_both);
  // ADREPORT(ssb_both);
  //  ADREPORT(d_obs);
  //  ADREPORT(logit_d_obs);

  return nll;
}