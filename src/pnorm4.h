//Log of dnorm
template<class Type>
Type ldnorm(Type x){
  Type ret = 0.5*(-x*x-log(2.0*M_PI));
  return ret;
}

//Atomic function for log(pnorm)
//Needed to get good stable derivatives for use with pnorm4
//Atomic double specifes the forward mode derivative, which is just pnorm(x,log.p=TRUE) = log(pnorm)
//Reverse is the explicit derivative of log(pnorm) = dnorm/pnorm = exp(log(dnorm)-log(pnorm)), more stable derivatives at
//extreme values.
TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   lpnorm1
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   ty[0] = atomic::Rmath::Rf_pnorm5(tx[0],0,1,1,1);
			   ,
			   // ATOMIC_REVERSE
			   Type pn = ty[0]; //Get the log-supporting pnorm
			   px[0] = exp(ldnorm(tx[0])-pn)*py[0]; //derivative of log(pnorm) in a more computatially stable way
			   )

//Type function to link to atomic
template<class Type>
Type lpnorm1(Type x){
  CppAD::vector<Type> tx(1);
  tx[0] = x;
  return lpnorm1(tx)[0];
}

//vectorize it
VECTORIZE1_t(lpnorm1);

//A version of pnorm taking 4 arguments,
//q quantile of the normal distribution
//mean the mean of the normal distribution
//sd the standard deviation
//give_log if 1, returns log(pnorm(q,mean,sd))
template<class Type>
Type pnorm4(Type q, Type mean = 0.,Type sd = 1.,int give_log = 1){
  if(give_log == 0){
    return pnorm(q,mean,sd);
  }else{
    Type m = (q-mean)/sd;
    return lpnorm1(m);
  }
}

VECTORIZE1_t(pnorm4);
VECTORIZE4_ttti(pnorm4);


//Version for censored upper and lower bounds

//logspace_sub to do log(exp(log(x))-exp(log(y))) safely
extern "C" {
    /* See 'R-API: entry points to C-code' (Writing R-extensions) */
    double Rf_logspace_sub (double logx, double logy);
}

//Here x is std. normal
//This is the forward double function needed by the atomic function
//Returns log(exp(log(pnorm(x+upper)))-exp(log(pnorm(x-lower)))) but in a much more compuationally stable way. 
double censored_bou(double x,double lowerb,double upperb){
  double upper,lower;
  if(x > 0){
    lower = atomic::Rmath::Rf_pnorm5(-(x+upperb),0,1,1,1);
    upper = atomic::Rmath::Rf_pnorm5(-(x-lowerb),0,1,1,1);
  }else{
    upper = atomic::Rmath::Rf_pnorm5(x+upperb,0,1,1,1);
    //Rcout << "Upper bound: " << upper << std::endl; 
    lower = atomic::Rmath::Rf_pnorm5(x-lowerb,0,1,1,1);
    //Rcout << "Lower bound: " << lower << std::endl;
  }
  double ret = Rf_logspace_sub(upper,lower);
  //Rcout << "What we ret: " << ret << std::endl;
  return ret;
}


TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   censored_b
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   ty[0] = censored_bou(tx[0],tx[1],tx[2]);
			   ,
			   // ATOMIC_REVERSE
			   Type cb = ty[0]; //Get the censored bounds from the forward mode.
			   Type d1 = ldnorm(tx[0]+tx[2]);
			   Type d2 = ldnorm(tx[0]-tx[1]);
			   Type f1 = d1-cb;
			   //Rcout << "This is f1: " << f1 << std::endl;
			   Type f2 = d2-cb;
			   //Rcout << "This is f2: " << f2 << std::endl;
			   px[0] = (exp(f1)-exp(f2))*py[0];
			   )

template<class Type>
Type censored_b(Type x,Type lowerb,Type upperb){
  CppAD::vector<Type> tx(3);
  tx[0] = x;
  tx[1] = lowerb;
  tx[2] = upperb;
  return censored_b(tx)[0];
}

//Function to calculate log(pnorm(x+upper,mu,sd)-pnorm(x-lower,mu,sd)) with stable derivatives
template<class Type>
Type censored_bounds(Type x,Type mu,Type sd,Type lower,Type upper){
  Type z = (x-mu)/sd;
  Type upperb = upper/sd;
  Type lowerb = lower/sd;
  return censored_b(z,lowerb,upperb);
}


		    
	 
  
  
