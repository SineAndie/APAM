#define TMB_LIB_INIT R_init_APAM

#include <TMB.hpp>
#include "pnorm4.hpp" //Atomic functions for censored likelihoods
#include <iostream>


template<class Type>
  Type objective_function<Type>::operator() ()
{

  //input data;
    DATA_MATRIX(M);
    DATA_MATRIX(weight);
    DATA_MATRIX(mat);
    DATA_MATRIX(midy_weight);
    DATA_VECTOR(index);
    DATA_VECTOR(olandings);
    DATA_IVECTOR(iyear);
    DATA_IVECTOR(iage);
    DATA_IVECTOR(isurvey);
    DATA_IVECTOR(isd);
    DATA_IVECTOR(is_year);
    DATA_VECTOR(fs);
    DATA_INTEGER(A);
    DATA_INTEGER(Y);
    DATA_INTEGER(Ns);
    DATA_INTEGER(NsF);
    DATA_INTEGER(NsSpan);
    DATA_ARRAY(crl);
    DATA_IVECTOR(isurvey1);
    DATA_VECTOR(landings_L);
    DATA_VECTOR(landings_U);
    DATA_VECTOR(log_landings);
    DATA_VECTOR(log_landings_L);
    DATA_VECTOR(log_landings_U);
    DATA_SCALAR(std_log_landings);
    DATA_VECTOR(log_lowerM);
    DATA_VECTOR(log_upperM);
	  DATA_VECTOR(d);
	  DATA_VECTOR(nll_wt);

    int n = index.size();
    Type one = 1.0;
    Type zero = 0.0;

  //define fixed parameters;
    PARAMETER_VECTOR(log_R);
    PARAMETER_ARRAY(m_q);
    PARAMETER_VECTOR(log_cv_index);
    PARAMETER(log_std_log_R);
    PARAMETER_ARRAY(log_F_mean);
    PARAMETER_ARRAY(log_std_log_F);
    PARAMETER_ARRAY(log_std_pe);
    PARAMETER_ARRAY(log_std_crl);
    PARAMETER_VECTOR(logit_ar_index_age);
    PARAMETER_VECTOR(logit_ar_logF);
    PARAMETER_VECTOR(logit_ar_crl);
    PARAMETER_VECTOR(logit_ar_pe);
  	PARAMETER(logit_ar_logRec);

  //define random effects
    PARAMETER_ARRAY(log_F_devt);
    PARAMETER_ARRAY(log_Nt);

  //to use when calculating local influence
	PARAMETER(h);

    array<Type> log_F_dev = log_F_devt.transpose();
    array<Type> log_N = log_Nt.transpose();

//  set bounds on parameters
    Type std_log_R = exp(log_std_log_R);
    vector<Type> cv_index = exp(log_cv_index);

    Type ar_logF_age = invlogit(logit_ar_logF(0));
    Type ar_logFone_year = invlogit(logit_ar_logF(1));
    Type ar_logF_year = invlogit(logit_ar_logF(2));
    Type ar_crl_age = invlogit(logit_ar_crl(0));
    Type ar_crl_year = invlogit(logit_ar_crl(1));
    Type ar_pe_age = invlogit(logit_ar_pe(0));
    Type ar_pe_year = invlogit(logit_ar_pe(1));
  	Type ar_logRec = invlogit(logit_ar_logRec);
    vector<Type> ar_index_age = exp(logit_ar_index_age)/(one + exp(logit_ar_index_age));

  //containers
    matrix<Type> N(Y,A);
    array<Type> pe(Y-1,A-1);
    vector<Type> log_Rec_dev(Y);
    vector<Type> log_Rec(Y);

    matrix<Type> EC(Y,A-4);
    matrix<Type> ECW(Y,A-4);

    vector<Type> C_tot(Y);
    vector<Type> CW_tot(Y);
    vector<Type> log_landings_pred(Y);
    vector<Type> std_landings_resid(Y);
    vector<Type> landings_resid(Y);

    vector<Type> Elog_index(n);
    vector<Type> Eindex(n);
    vector<Type> resid_index(n);
    vector<Type> std_resid_index(n);

    matrix<Type> p_ya(A-4,Y);
    matrix<Type> pya_sum(A-4,Y);
    matrix<Type> pi_ya(A-4,Y);
  	matrix<Type> pred_crl(Y,A-5);
    matrix<Type> resid_crl(Y,A-5);
    matrix<Type> std_resid_crl(Y,A-5);
    array<Type> std_crl(Y,A-5);
    array<Type> std_pe(Y-1,A-1);
    array<Type> std_log_F(Y,A-4);

  //SD report objects
    matrix<Type> B_matrix(Y,A);
    matrix<Type> SSB_matrix(Y,A-4);
    vector<Type> biomass(Y);
    vector<Type> log_biomass(Y);
    vector<Type> ssb(Y);
    vector<Type> log_ssb(Y);
    vector<Type> aveF_914(Y);
    vector<Type> log_aveF_914(Y);


    //initialize the negative log likelihood
    Type nll = zero;
    vector<Type> jnll(8);
    jnll.setZero();
   using namespace density;

	matrix<Type> F(Y,A-4);
	matrix<Type> Z(Y,A);
	array<Type> F_dev(Y,A-4);
	array<Type> F_mean(Y,2);
	matrix<Type> log_F = log(F.array());
	matrix<Type> M_new(Y,A);

  int i,j;
  for(i = 0;i < Y;++i){
    for(j = 0;j < A-5;++j){std_crl(i,j) = exp(log_std_crl(i,j));}
  }
  for(i = 0;i < Y-1;++i){
    for(j = 0;j < A-1;++j){std_pe(i,j) = exp(log_std_pe(i,j));}
  }

  for(i = 0;i < Y;++i){
    for(j = 0;j < A-4;++j){std_log_F(i,j) = exp(log_std_log_F(i,j));
      F_dev(i,j) = exp(log_F_dev(i,j));}
  }

  for(i = 0;i < Y;++i){
     log_F(i,0) =  log_F_mean(i,0)+log_F_dev(i,0);
     log_F(i,1) =  log_F_mean(i,1)+log_F_dev(i,1);
     F(i,0) = exp(log_F(i,0));
     F(i,1) = exp(log_F(i,1));

    for(j = 2;j < A-4;++j){
      F(i,j) = F(i,j-1)+F_dev(i,j);
      log_F(i,j) = log(F(i,j));
    }
  }
  F_mean = exp(log_F_mean);

  int pc = 0;

  for(int j = 0;j < A;++j){
  for(i = 0;i < Y;++i){
    M_new(i,j) = M(i,j)+(h*d(pc));
    pc = pc+1;
    }
    }

	for(int i = 0;i < Y;++i){
    for(int j = 0;j < A;++j){
      if(j<4){Z(i,j) =  M_new(i,j);}
      if(j>=4){Z(i,j)= F(i,j-4)+M_new(i,j);}
    }}

  //compute process errors
    for(int i = 1;i < Y;++i){
      for(int j = 1;j < A-1;++j){
        pe(i-1,j-1) = log_N(i,j) - log_N(i-1,j-1) + Z(i-1,j-1);
      }
      int j=A-1;
      pe(i-1,j-1) = log_N(i,j) - log(exp(log_N(i-1,j-1)-Z(i-1,j-1))+exp(log_N(i-1,j)-Z(i-1,j)));
    }

      //calculate recruitment deviations;
    for(int i = 0;i < Y;++i){
      if(i<=32){log_Rec_dev(i)= log_N(i,0) - log_R(0);}
      if(i>32){log_Rec_dev(i)= log_N(i,0) - log_R(1);}
      log_Rec(i) = log_N(i,0);
    }

    for(int i = 0;i < Y;++i){
      for(int j = 0;j < A;++j){
        N(i,j)=exp(log_N(i,j));}}

  //Baranov catch equation predictions and residuals
    for(int i = 0;i < Y;++i){
      C_tot(i) = zero;
  	  CW_tot(i) = zero;
        for(int j = 0;j < A-4;++j){
          EC(i,j) = N(i,j+4)*((one - exp(-one*Z(i,j+4)))*F(i,j)/Z(i,j+4));
          ECW(i,j) = EC(i,j)*midy_weight(i,j);
          C_tot(i) += EC(i,j);
          CW_tot(i) += ECW(i,j);
    } }

    log_landings_pred = log(CW_tot);
    landings_resid = log_landings - log_landings_pred;
    std_landings_resid = landings_resid/std_log_landings;

   //age composition catch
    for(int i = 0;i < Y;++i){
        for(int j = 0;j < A-4;++j){
          p_ya(j,i) = EC(i,j)/C_tot(i);
      }}

    Type total;
    for(int i = 0;i <Y;++i){
      total=zero;
  	    pya_sum(0,i) = one;
        for(int j = 1;j < A-4;++j){
          total+=p_ya(j-1,i);
          pya_sum(j,i) = (one-total);
      }}

  	  for(int i = 0;i <Y;++i){
        for(int j = 0;j < A-4;++j){
  			pi_ya(j,i) = p_ya(j,i)/pya_sum(j,i);}}


    for(int i = 0;i<Y;i++){
      for(int j = 0;j<A-5;j++){
        pred_crl(i,j)= log(pi_ya(j,i)/(one - pi_ya(j,i)));
        resid_crl(i,j) = crl(i,j) - pred_crl(i,j);
        std_resid_crl(i,j) = resid_crl(i,j)/std_crl(i,j);}}

 //Survey index predictions, and residuals;
  vector<Type> std_index_vec(n);
  matrix<Type> mresid_index(Ns,A);
  matrix<Type> msd_index(Ns,A);
  matrix<Type> q(Ns,A);


  int ia,iy,iy1,is;

  for(i = 0;i < Ns;++i){
    for(j = 1;j < A-1;++j){
      m_q(i,j) = exp(m_q(i,j));
    }}


  for(i = 0;i < Ns;++i){
     q(i,0) = m_q(i,0);

    for(j=1;j<A-1;++j){

    //Fall trawl split
    if(i<5&&j<4){q(i,j) = q(i,j-1)+m_q(i,j);}
    if(i>=5&&i<NsF){q(i,j) = q(i,j-1)+m_q(i,j);}

    //Spanish
    if(i>=NsF&&i<(NsF+NsSpan)){q(i,j) = q(i,j-1)+m_q(i,j);}

    //Spring trawl split
    if(i>=(NsF+NsSpan)&&i<(NsF+NsSpan+11)&&j<4){q(i,j) = q(i,j-1)+m_q(i,j);}
    if(i>=(NsF+NsSpan+11)){q(i,j) = q(i,j-1)+m_q(i,j);}
    }
    int j=A-1;
    q(i,j) = q(i,j-1);}

  for(i = 0;i < Ns;++i){
    for(j=1;j<A-1;++j){
    if(i<5&&j>=4){q(i,j) = q(i+5,j-1)+m_q(i+5,j);}
    if(i>=(NsF+NsSpan)&&i<(NsF+NsSpan+11)&&j>=4){q(i,j) = q(i+11,j-1)+m_q(i+11,j);}
    }
    int j=A-1;
    if(i<5&&j==A-1){q(i,j) = q(i+5,j-1);}
    if(i>=(NsF+NsSpan)&&i<(NsF+NsSpan+11)&&j==A-1){q(i,j) = q(i+11,j-1);}
    }

  for(i = 0;i < n;++i){
    ia = iage(i);
    iy = iyear(i);
    is = isd(i);
    iy1 = is_year(i);

    Elog_index(i) = q(iy1,ia)  + log_N(iy,ia) - fs(i)*Z(iy,ia);
    Eindex(i) = exp(Elog_index(i));
    std_index_vec(i) = cv_index(is)*Eindex(i);
    resid_index(i) = index(i) - Eindex(i);
    std_resid_index(i) = resid_index(i)/std_index_vec(i);

    mresid_index(iy1,ia) = resid_index(i);
    msd_index(iy1,ia) = std_index_vec(i);
  }

   //NEGATIVE LOGLIKELIHOODS
   //Index OBSERVATION MODEL

    vector<Type> del(A);
    vector<Type> sd_del(A);

    for(i = 0;i < Ns;++i){
      iy = isurvey1(i);
      del = vector<Type>(mresid_index.row(i));
      sd_del = vector<Type>(msd_index.row(i));

      if(isurvey1(i)==0){jnll(0) += VECSCALE(AR1(ar_index_age(iy)),sd_del)(del);}
      if(isurvey1(i)==1){jnll(1) += VECSCALE(AR1(ar_index_age(iy)),sd_del)(del);}
      if(isurvey1(i)==2){jnll(2) += VECSCALE(AR1(ar_index_age(iy)),sd_del)(del);}
    }

    jnll(0) = jnll(0)*nll_wt(0);
    jnll(1) = jnll(1)*nll_wt(1);
    jnll(2) = jnll(2)*nll_wt(2);

    nll += jnll(0);
    nll += jnll(1);
    nll += jnll(2);

   //Landings censored nll;
      for(int i = 0;i < Y;++i){
        jnll(3) -= censored_bounds(log_landings(i),log_landings_pred(i),std_log_landings,
        -log_lowerM(i),log_upperM(i));
    }

    jnll(3) = jnll(3)*nll_wt(3);
    nll += jnll(3);

  	array<Type> temp(Y,A-5);
    temp = resid_crl.array();
    jnll(4)  += VECSCALE(SEPARABLE(AR1(ar_crl_age),AR1(ar_crl_year)),std_crl)(temp);

    jnll(4) = jnll(4)*nll_wt(4);
    nll += jnll(4);
   //PROCESS MODEL

  //Log recruitS;
   jnll(5) += SCALE(AR1(ar_logRec),std_log_R)(log_Rec_dev);

    nll += jnll(5);
      //Log F
      //yearxage correlation on first:last ages; F split 1994

      vector<Type> delF1(Y);
      vector<Type> delstdF1(Y);

      for(int j = 0;j < Y;++j){delF1(j) = log_F_dev(j,0);
        delstdF1(j) = std_log_F(j,0);}

      jnll(6) += VECSCALE(AR1(ar_logFone_year),delstdF1)(delF1);

      array<Type> log_F1_dev1(Y,A-5);
      array<Type> std_F_dev1(Y,A-5);

      for(int i = 1;i < A-4;++i){
        for(int j = 0;j < Y;++j){log_F1_dev1(j,i-1) = log_F_dev(j,i);
          std_F_dev1(j,i-1) =  std_log_F(j,i);}
      }

      jnll(6) += VECSCALE(SEPARABLE(AR1(ar_logF_age),AR1(ar_logF_year)),std_F_dev1)(log_F1_dev1);

     nll += jnll(6);
  //Process error

   jnll(7) += VECSCALE(SEPARABLE(AR1(ar_pe_age),AR1(ar_pe_year)),std_pe)(pe);
   nll += jnll(7);

   //Useful output
  //Biomass and SSB
    for(int i=0;i<Y;++i){
      for(int j=0;j<A;++j){
        B_matrix(i,j) = weight(i,j)*N(i,j);
    }}

    for(int i=0;i<Y;i++){
      for(int j=0;j<A-4;++j){
        SSB_matrix(i,j) = mat(i,j)*B_matrix(i,j+4);
    }}

    for(int i = 0;i < Y;++i){
      biomass(i) = zero;
      ssb(i) = zero;
        for(int j = 0;j < A;++j){
          biomass(i) += B_matrix(i,j);
            if(j>=4){ssb(i) += SSB_matrix(i,j-4);}
    }}

    log_biomass = log(biomass);
    log_ssb = log(ssb);

  //pop size weighted ave F(9-14)
    Type tni;

    for(int i = 0;i < Y;++i){
      aveF_914(i) = zero;
      tni = zero;
        for(int j =8 ;j < 14;++j){
          aveF_914(i) += F(i,j-4)*N(i,j);
          tni += N(i,j);
        }
      aveF_914(i) = aveF_914(i)/tni;
    }

    log_aveF_914 = log(aveF_914);


  	REPORT(N);
    REPORT(F);
    REPORT(Z);
	REPORT(M_new);
    REPORT(B_matrix);
    REPORT(SSB_matrix);
    REPORT(biomass);
    REPORT(ssb);
    REPORT(aveF_914);
    REPORT(log_Rec);
    REPORT(C_tot);
    REPORT(pred_crl);
    REPORT(resid_crl);
    REPORT(std_resid_crl);
    REPORT(p_ya);
    REPORT(pya_sum);
    REPORT(pi_ya);
    REPORT(CW_tot);
    REPORT(EC);
  	REPORT(ECW);
    REPORT(landings_resid);
    REPORT(log_landings_pred);
    REPORT(std_landings_resid);
    REPORT(std_log_F);
    REPORT(std_pe);
    REPORT(std_crl);
    REPORT(cv_index);
    REPORT(std_index_vec);
	REPORT(std_log_R);

    REPORT(ar_crl_age);
    REPORT(ar_crl_year);
    REPORT(ar_logF_age);
    REPORT(ar_logF_year);
    REPORT(ar_logFone_year);
    REPORT(ar_index_age);
    REPORT(ar_logRec);
    REPORT(ar_pe_year);
    REPORT(ar_pe_age);

    REPORT(logit_ar_index_age);
    REPORT(logit_ar_logF);
    REPORT(logit_ar_crl);
    REPORT(logit_ar_pe);
  	REPORT(logit_ar_logRec);

    REPORT(Elog_index);
    REPORT(Eindex);
    REPORT(resid_index);
    REPORT(std_resid_index);
    REPORT(mresid_index);
    REPORT(msd_index);
    REPORT(m_q);

    REPORT(log_N);
    REPORT(log_Rec_dev);
    REPORT(log_R);
    REPORT(log_F_mean);
    REPORT(log_F_dev);
    REPORT(F_mean);
    REPORT(F_dev);
    REPORT(log_F1_dev1);
    REPORT(log_F);
    REPORT(pe);
    REPORT(q);
    REPORT(jnll);

    ADREPORT(log_landings_pred);
    ADREPORT(log_biomass);
    ADREPORT(log_ssb);
    ADREPORT(log_aveF_914);
    ADREPORT(log_Rec);
    ADREPORT(log_F_mean);
    ADREPORT(F_mean);
  	ADREPORT(log_N);
  	ADREPORT(ECW);
  	ADREPORT(m_q);

	ADREPORT(std_log_F);
	ADREPORT(std_pe);
	ADREPORT(std_crl);
	ADREPORT(std_log_R);
	ADREPORT(cv_index);
	ADREPORT(ar_logF_age);
	ADREPORT(ar_logF_year);
	ADREPORT(ar_logFone_year);
	ADREPORT(ar_pe_year);
	ADREPORT(ar_pe_age);
	ADREPORT(ar_logRec);
	ADREPORT(ar_index_age);
	ADREPORT(ar_crl_year);
	ADREPORT(ar_crl_age);
	ADREPORT(ssb);
	ADREPORT(aveF_914);

  return nll;
  }
