// Estimate MSY and PGY using TMB
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(naa_mat);
  //  DATA_SCALAR(deviance_init);
  DATA_ARRAY(SR_mat); // 1: HS, 2: BH, 3: RI
  DATA_ARRAY(rec_par_a_mat);
  DATA_ARRAY(rec_par_b_mat); 
  DATA_ARRAY(rec_par_rho_mat);
  DATA_SCALAR(bias_corrected_mean);
  DATA_ARRAY(rec_resid_mat); 
  DATA_ARRAY(waa_mat);
  DATA_ARRAY(maa_mat);  
  DATA_ARRAY(M_mat);
  DATA_ARRAY(faa_mat);  
  DATA_INTEGER(Pope);
  DATA_INTEGER(total_nyear);  
  DATA_INTEGER(future_initial_year);
  DATA_INTEGER(start_F_year);    
  DATA_INTEGER(nsim);
  DATA_INTEGER(nage);
  DATA_INTEGER(recruit_age);    
  DATA_INTEGER(obj_catch); // 0: mean, 1: geomean
  DATA_INTEGER(objective); // 0: MSY, 1: PGY, 2: percentB0 or Bempirical
  DATA_SCALAR(objective_value); // Used for objective 1-2
  //  DATA_VECTOR(Fcurrent);
  //  DATA_INTEGER(Fcurrent_year);
  //  DATA_SCALAR(num_to_mass_scale);

  PARAMETER(x); //x = log(multiplier)
  
  // Type bias_corrected_mean = -Type(0.5)*pow(rec_par_sd,Type(2.0))/(Type(1.0)-pow(rec_par_rho, Type(2.0)));
  array<Type> F_mat(nage,total_nyear,nsim);
  array<Type> N_mat(nage,total_nyear,nsim);
  //  array<Type> rec_deviance_mat(total_nyear,nsim);
  array<Type> spawner_mat(total_nyear,nsim);
  array<Type> catch_mat(nage,total_nyear,nsim);
  F_mat.setZero();  
  N_mat.setZero();
  spawner_mat.setZero();
  //  rec_deviance_mat.setZero();  
  catch_mat.setZero();      
  
  for(int i=0; i<nsim; i++) { // input initial number
    for(int t=0; t<future_initial_year+1; t++){ // この時点で、future_initial_year前までにnaaがデータがあることを確認
      for(int iage=0; iage<nage; iage++) {
	N_mat(iage,t,i) = naa_mat(iage,t,i);
	//	F_mat(iage,t,i) = faa_mat(iage,t,i);
	spawner_mat(t,i) += N_mat(iage,t,i)*waa_mat(iage,t,i)*maa_mat(iage,t,i);
	//	rec_devinace_mat(t,i) = rec_resid_mat(t,i)
      }
    }}

  // Matrix of Fishing mortality (use VPA estimation)
  for(int i=0; i<nsim; i++) { //replication of simulation 
    for(int t=0; t<start_F_year-1; t++) {
       for(int iage=0; iage<nage; iage++) {
	 F_mat(iage,t,i) = faa_mat(iage, t, i);
       }}}  

  // Matrix of Fishing mortality (use future value)
  for(int i=0; i<nsim; i++) { //replication of simulation 
     for(int t=start_F_year-1; t<total_nyear; t++) {
       for(int iage=0; iage<nage; iage++) {
	 F_mat(iage,t,i) = exp(x)*faa_mat(iage,t,i);
       }}}
  
  // Population dynamics
  for(int i=0; i<nsim; i++) {
    for(int t=future_initial_year-1; t<total_nyear; t++) {

      // summing up spawners
      spawner_mat(t,i) = 0; 
      for(int a=1; a<nage; a++) {
	spawner_mat(t,i) += N_mat(a,t,i)*waa_mat(a,t,i)*maa_mat(a,t,i); 
      }
	
      // update recruitment except for t=initial year (t=0)
      if(t>future_initial_year-1){      
	if(SR_mat(t,i) == 1) { //Hockey-stick
	  vector<Type> rec_pred(2);
	  rec_pred(0) = spawner_mat(t-recruit_age,i)*rec_par_a_mat(t,i);
	  rec_pred(1) = rec_par_b_mat(t,i)*rec_par_a_mat(t,i);
	  N_mat(0,t,i) = min(rec_pred);
	}
	if(SR_mat(t,i) == 2) { //Beverton-Holt
	  N_mat(0,t,i) = rec_par_a_mat(t,i)*spawner_mat(t-recruit_age,i)/(1+rec_par_b_mat(t,i)*spawner_mat(t-recruit_age,i));
	}
	if(SR_mat(t,i) == 3) { //Ricker
	  N_mat(0,t,i) = rec_par_a_mat(t,i)*spawner_mat(t-recruit_age,i)*exp(-rec_par_b_mat(t,i)*spawner_mat(t-recruit_age,i));
	}
	//	N_mat(0,t,i) = N_mat(0,t,i)*exp(rec_par_rho_mat(t-1,i)*rec_resid_mat(t-1,i)+
	//					rec_resid_mat(t,i)+bias_corrected_mean);
	N_mat(0,t,i) = N_mat(0,t,i)*exp(rec_resid_mat(t,i));	
      }

      // forward calculation 
      if(t<(total_nyear-1)){
	for(int iage=0; iage<(nage-1); iage++) {
	  N_mat(iage+1,t+1,i) = N_mat(iage,t,i)*exp(-M_mat(iage,t,i)-F_mat(iage,t,i));
	}
        N_mat(nage-1,t+1,i) += N_mat(nage-1,t,i)*exp(-M_mat(nage-1,t,i)-F_mat(nage-1,t,i));
      }
    }
  }
  
  // Catch equation
  for(int i=0; i<nsim; i++) {
    for(int t=0; t<total_nyear; t++) {
      for(int a=0; a<nage; a++){
   	if(Pope) {
   	  catch_mat(a,t,i) = waa_mat(a,t,i)*N_mat(a,t,i)*exp(-Type(0.5)*M_mat(a,t,i))*(1-exp(-F_mat(a,t,i)));
   	}else{
   	  catch_mat(a,t,i) = waa_mat(a,t,i)*N_mat(a,t,i)*(1-exp(-M_mat(a,t,i)-F_mat(a,t,i)))*F_mat(a,t,i)/(M_mat(a,t,i)+F_mat(a,t,i));
	}
      }
    }}
  
  // Get Catch or SSB in the final year
  Type obj = 0;
    for(int i=0; i<nsim; i++) {
      for(int a=0; a<nage; a++) {    
        if(objective < 2) {
    	if(obj_catch == 0) {
    	  obj += catch_mat(a,total_nyear-1,i);
    	}else{
    	  obj += log(catch_mat(a,total_nyear-1,i));
    	}
        }else{
    	if(obj_catch == 0) {
    	  obj += spawner_mat(total_nyear-1,i);
    	}else{
    	  obj += log(spawner_mat(total_nyear-1,i));
    	}
        }
      }}
   obj /= nsim;

  if(obj_catch == 1) { 
    obj = exp(obj); // geomean
  }
  
  if(objective == 0) { //  MSY
    obj = log(obj);
    obj *= -1; // negative log catch
  }else{ // PGY, percentB0, and Bempirical
    obj = pow(log(obj/objective_value), Type(2.0));
  }

  REPORT(F_mat);
  REPORT(N_mat);
  REPORT(spawner_mat);
  REPORT(catch_mat);
  REPORT(future_initial_year);    
  //  ADREPORT(obj);

  return obj;
}

  
