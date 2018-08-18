#include <TMB.hpp>
// Multispecies CMSY

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ========== Inputs ============================
  DATA_INTEGER(n_t); //total number of years
  DATA_INTEGER(n_s); //total number of species
  DATA_MATRIX(C_ts); //catch data with years along rows and species along columns
  // DATA_INTEGER(choose_ref); //which is the reference stock? 1 through n_s
  // DATA_INTEGER(like_type); //0 == original idea, 1==sum of squares for effort

  //prior inputs
  DATA_INTEGER(r_prior); // on == 1, off == 0
  DATA_VECTOR(r_means);
  DATA_VECTOR(r_sds);
  DATA_INTEGER(K_prior); // on == 1, off == 0
  DATA_VECTOR(K_means);
  DATA_VECTOR(K_sds);
  DATA_INTEGER(delta_prior); // on == 1, off == 0
  DATA_VECTOR(delta_means);
  DATA_VECTOR(delta_sds);


  // ======== Parameters ==========================
  PARAMETER_VECTOR(logr);
  PARAMETER_VECTOR(logK);
  PARAMETER_VECTOR(delta_s);
  PARAMETER_VECTOR(logq);
  PARAMETER(logsigmaC);
  PARAMETER_VECTOR(lE_t);

  // ============ Global values ===================

  using namespace density;
  Type jnll = 0;

  //----- transformations
  vector<Type> r_s(n_s);
  vector<Type> K_s(n_s);
  vector<Type> q_s(n_s);
  for(int s=0;s<n_s;s++){
    r_s(s) = exp(logr(s));
    K_s(s) = exp(logK(s));
    q_s(s) = exp(logq(s));
  }
  Type sigmaC = exp(logsigmaC);

  vector<Type> E_t(n_t);
  for(int t=0;t<n_t;t++){
    E_t(t) = exp(lE_t(t));
  }

  //adjust reference stock integer to be from 0 through (n_s - 1)
  // int find_ref = choose_ref - 1;

  //posfun requirements
  Type pen = 0;
  Type eps = 1e-3;

  // ============ Model ===========================

  //Calculate biomass
  matrix<Type> B_ts(n_t,n_s);
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      if(t==0) B_ts(t,s) = K_s(s) * delta_s(s);
      if(t>0) B_ts(t,s) = B_ts(t-1,s) + r_s(s) * B_ts(t-1,s) * (1 - B_ts(t-1,s)/K_s(s)) - C_ts(t-1,s);
      B_ts(t,s) = posfun(B_ts(t,s), eps, pen);
    }
  }

  // Predicted catch
  matrix<Type> Cpred_ts(n_t,n_s);
  matrix<Type> U_ts(n_t,n_s);
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      Cpred_ts(t,s) = q_s(s)*E_t(t)*B_ts(t,s);
    }
  }

  //Derive predicted exploitation fractions
  matrix<Type> Uobs_ts(n_t,n_s);
  matrix<Type> Upred_ts(n_t,n_s);
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      Uobs_ts(t,s) = C_ts(t,s) / B_ts(t,s);
      Upred_ts(t,s) = q_s(s) * E_t(t);  
    }
  }


  //reference stock
  // vector<Type> Uref_t(n_t);
  // Uref_t = U_ts.col(find_ref);

  // ============ Likelihood ===========================

  // likelihood of catch
  matrix<Type> nll_ts(n_t,n_s);
  nll_ts.setZero();
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      nll_ts(t,s) -= dnorm(C_ts(t,s), Cpred_ts(t,s), sigmaC, true);
      // nll_ts(t,s) -= dnorm(Uobs_ts(t,s), Upred_ts(t,s), sigmaC, true);

    }
  }

  //priors
  matrix<Type> nll_sp(n_s,3);
  nll_sp.setZero();
  if(r_prior==1){
    for(int s=0;s<n_s;s++){
      nll_sp(s,0) -= dnorm(r_s(s), r_means(s), r_sds(s), true);
    }
  }
  if(K_prior==1){
    for(int s=0;s<n_s;s++){
      nll_sp(s,1) -= dnorm(K_s(s), K_means(s), K_sds(s), true);
    }
  }
  if(delta_prior==1){
    for(int s=0;s<n_s;s++){
      nll_sp(s,2) -= dnorm(delta_s(s), delta_means(s), delta_sds(s), true);
    }
  }

  //likelihood components
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
    // likelihood penalty if B negative
    jnll_comp(0) += pen;
    // penalized likelihood from priors
    jnll_comp(1) += sum(nll_sp);
    // exploitation ratio likelihood
    jnll_comp(2) += sum(nll_ts);

  jnll = sum(jnll_comp);

  // ============ Derived values ==================

  vector<Type> Bmsy_s(n_s);
  vector<Type> Umsy_s(n_s);
  for(int s=0;s<n_s;s++){
    Bmsy_s(s) = K_s(s)/2;
    Umsy_s(s) = r_s(s)/2;
  }

  matrix<Type> BBmsy_ts(n_t,n_s);
  matrix<Type> UUmsy_cb_ts(n_t,n_s);
  matrix<Type> UUmsy_qE_ts(n_t,n_s);
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      BBmsy_ts(t,s) = B_ts(t,s)/Bmsy_s(s);
      UUmsy_qE_ts(t,s) = Upred_ts(t,s)/Umsy_s(s);
      UUmsy_cb_ts(t,s) = Uobs_ts(t,s)/Umsy_s(s);
    }
  }

  // ============ Reporting section ===============

  // standard errors
  ADREPORT(r_s);
  ADREPORT(K_s);
  ADREPORT(B_ts);
  ADREPORT(Uobs_ts);
  ADREPORT(Upred_ts);
  ADREPORT(E_t);
  ADREPORT(BBmsy_ts);
  ADREPORT(UUmsy_cb_ts);
  ADREPORT(UUmsy_qE_ts);
  ADREPORT(Cpred_ts);

  //-------- check inputs
  REPORT(n_t);
  REPORT(n_s);
  // REPORT(choose_ref);
  // REPORT(find_ref);

  //-------- parameters
  REPORT(r_s);
  REPORT(K_s);
  REPORT(sigmaC);
  REPORT(q_s);
  REPORT(E_t);
  REPORT(delta_s);

  //--------- derived values
  REPORT(B_ts);
  REPORT(Uobs_ts);
  REPORT(Upred_ts);
  REPORT(Cpred_ts);
  REPORT(Bmsy_s);
  REPORT(BBmsy_ts);
  REPORT(Umsy_s);
  REPORT(UUmsy_cb_ts);
  REPORT(UUmsy_qE_ts);

  //--------- likelihood components

  REPORT(K_means);
  REPORT(K_sds);
  REPORT(r_means);
  REPORT(r_sds);

  // REPORT(Z_ts);
  // REPORT(Zavg_s);
  // REPORT(qavg_s);
  // REPORT(sq_ts);
  // REPORT(ssq_s);
  REPORT(nll_ts);
  REPORT(nll_sp);
  REPORT(pen);
  REPORT(jnll_comp);

  REPORT(jnll);
  return(jnll);

  //Calculate log-ratio of secondary stocks to reference stock
  // matrix<Type> Z_ts(n_t,n_s);
  // vector<Type> Zavg_s(n_s);
  // vector<Type> Zsum_s(n_s);
  // Zsum_s.setZero();
  // vector<Type> qavg_s(n_s);
  // matrix<Type> E_ts(n_t,n_s);

  // for(int s=0;s<n_s;s++){
  //     for(int t=0;t<n_t;t++){ 
  //       Z_ts(t,s) = log(U_ts(t,s) / Uref_t(t));
  //       Zsum_s(s) += Z_ts(t,s);
  //     }
  //     //calculate average for each stock
  //     Zavg_s(s) = Zsum_s(s)/Type(n_t);
  //     qavg_s(s) = exp(Zavg_s(s));

  //     for(int t=0;t<n_t;t++){
  //       E_ts(t,s) = U_ts(t,s)/qavg_s(s);
  //     }
  // }

  // vector<Type> Eref_t(n_t);
  // Eref_t = E_ts.col(find_ref);

  // //Calculate SSQ of log-ratio and NLL
  // matrix<Type> sq_ts(n_t,n_s);
  // // vector<Type> ssq_s(n_s);
  // // ssq_s.setZero();
  // matrix<Type> nll_ts(n_t,n_s);
  // nll_ts.setZero();
  // for(int s=0;s<n_s;s++){
  //   for(int t=0;t<n_t;t++){

  //     if(like_type==0){
  //       sq_ts(t,s) = pow((E_ts(t,s) - U_ts(t,s)), 2);
  //       if(s != find_ref) nll_ts(t,s) = Type(-1) * dnorm(sq_ts(t,s), Type(0.0), sigma, true);
  //     }
  //     if(like_type==1){
  //       sq_ts(t,s) = pow((E_ts(t,s) - Eref_t(t)), 2);
  //       if(s != find_ref) nll_ts(t,s) = sq_ts(t,s);
  //     }
  //     // sq_ts(t,s) = pow((E_ts(t,s) - Uref_t(t)),2);
  //   }
  // }


}