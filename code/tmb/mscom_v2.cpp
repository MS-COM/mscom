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
  PARAMETER_VECTOR(logr_s);
  PARAMETER_VECTOR(logK_s);
  PARAMETER_VECTOR(logp_s);
  PARAMETER_VECTOR(delta_s);
  PARAMETER_VECTOR(logq_s);
  PARAMETER(logsigmaC);

  // ======== Random effects =====================
  // PARAMETER_VECTOR(lE_t);
  PARAMETER(logmuE);
  PARAMETER(logsigmaE);
  PARAMETER_VECTOR(Eps_t);

  // ============ Global values ===================

  using namespace density;
  Type jnll = 0;

  //----- transformations
  vector<Type> r_s(n_s);
  vector<Type> K_s(n_s);
  vector<Type> p_s(n_s);
  vector<Type> q_s(n_s);
  vector<Type> msy_s(n_s);
  vector<Type> Bmsy_s(n_s);
  vector<Type> Umsy_s(n_s);
  for(int s=0;s<n_s;s++){
    r_s(s) = exp(logr_s(s));
    K_s(s) = exp(logK_s(s));
    p_s(s) = exp(logp_s(s));
    q_s(s) = exp(logq_s(s));
    msy_s(s) = (r_s(s)*K_s(s))/(pow((p_s(s)+1), (p_s(s)+1)/p_s(s)));
    Bmsy_s(s) = K_s(s) * pow((1/(p_s(s)+1)), 1/p_s(s));
    Umsy_s(s) = msy_s(s)/Bmsy_s(s);
  }
  Type sigmaC = exp(logsigmaC);
  Type sigmaE = exp(logsigmaE);

  vector<Type> lE_t(n_t);
  vector<Type> E_t(n_t);
  for(int t=0;t<n_t;t++){
    lE_t(t) = logmuE + Eps_t(t) - pow(sigmaE,2)/Type(2);
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
      if(t>0) B_ts(t,s) = B_ts(t-1,s) + (r_s(s)/p_s(s))*B_ts(t-1,s)*(1-pow(B_ts(t-1,s)/K_s(s), p_s(s))) - C_ts(t-1,s);
      B_ts(t,s) = posfun(B_ts(t,s), eps, pen);
    }
  }

  // Predicted catch
  matrix<Type> Cpred_ts(n_t,n_s);
  matrix<Type> Upred_ts(n_t,n_s);
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      // Cpred_ts(t,s) = q_s(s)*E_t(t)*B_ts(t,s);
      Upred_ts(t,s) = q_s(s) * E_t(t);
      Cpred_ts(t,s) = Upred_ts(t,s) * B_ts(t,s);
    }
  }

  //Derive predicted exploitation fractions
  matrix<Type> Uobs_ts(n_t,n_s);
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      Uobs_ts(t,s) = C_ts(t,s) / B_ts(t,s);
    }
  }


  //reference stock
  // vector<Type> Uref_t(n_t);
  // Uref_t = U_ts.col(find_ref);

  // ============ Likelihood ===========================

  // likelihood of catch
  matrix<Type> NLL_Catch_ts(n_t,n_s);
  NLL_Catch_ts.setZero();
  for(int s=0;s<n_s;s++){
    for(int t=0;t<n_t;t++){
      NLL_Catch_ts(t,s) -= dnorm(C_ts(t,s), Cpred_ts(t,s), sigmaC, true);
      // NLL_Catch_ts(t,s) -= dnorm(Uobs_ts(t,s), Upred_ts(t,s), sigmaC, true);
    }
  }

  //priors
  matrix<Type> NLL_Priors_sp(n_s,3);
  NLL_Priors_sp.setZero();
  if(r_prior==1){
    for(int s=0;s<n_s;s++){
      NLL_Priors_sp(s,0) -= dnorm(r_s(s), r_means(s), r_sds(s), true);
    }
  }
  if(K_prior==1){
    for(int s=0;s<n_s;s++){
      NLL_Priors_sp(s,1) -= dnorm(K_s(s), K_means(s), K_sds(s), true);
    }
  }
  if(delta_prior==1){
    for(int s=0;s<n_s;s++){
      NLL_Priors_sp(s,2) -= dnorm(delta_s(s), delta_means(s), delta_sds(s), true);
    }
  }

 //random effect
  vector<Type> NLL_RandEff_t(n_t);
  NLL_RandEff_t(0) = dnorm(Eps_t(0), Type(0.0), sigmaE, true);
  for(int t=1;t<n_t;t++){
    NLL_RandEff_t(t) = dnorm(Eps_t(t), Eps_t(t-1), sigmaE, true);
  }

  //likelihood components
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
    // likelihood penalty if B negative
    jnll_comp(0) += pen;
    // penalized likelihood from priors
    jnll_comp(1) += sum(NLL_Priors_sp);
    // exploitation ratio likelihood
    jnll_comp(2) += sum(NLL_Catch_ts);
    // effort random effect
    jnll_comp(3) += sum(NLL_RandEff_t);

  jnll = sum(jnll_comp);

  // // ============ Derived values ==================

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
  // ADREPORT(U_t);
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
  REPORT(p_s);
  REPORT(sigmaC);
  REPORT(q_s);
  REPORT(E_t);
  // REPORT(U_t);
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
  REPORT(NLL_Catch_ts);
  REPORT(NLL_Priors_sp);
  REPORT(pen);
  REPORT(jnll_comp);

  REPORT(jnll);
  return(jnll);


}
