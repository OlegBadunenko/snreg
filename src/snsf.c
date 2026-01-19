// this is a package version
#include <R.h>
#include <Rmath.h>

# include "owen.h"

# include "asa076.h"

#ifndef MY_OMP_H_
#define MY_OMP_H_

#ifdef _OPENMP
#include <omp.h>
#if _OPENMP >= 201307
#define OMP_VER_4
#endif
#endif

// Insert SIMD pragma if supported
#ifdef OMP_VER_4
#define SAFE_SIMD _Pragma("omp simd")
#define SAFE_FOR_SIMD _Pragma("omp for simd")
#else
#define SAFE_SIMD 
#define SAFE_FOR_SIMD
#endif

#endif

// prod == 1  => production function
// prod == -1 => cost function

double pnormstd(double x)
{
  return 0.5 * erfc(-x * M_SQRT1_2);
}

double dnormstd(double x)
{
  return 0.5 * M_2_SQRTPI * M_SQRT1_2 * exp(- x * x * 0.5);
}

void pnormstdC(double *x)
{
  *x = 0.5 * erfc(-*x * M_SQRT1_2);
}

void dnormstdC(double *x)
{
  *x = 0.5 * M_2_SQRTPI * M_SQRT1_2 * exp(- *x * *x * 0.5);
}

/*
TOwen.gr.2 <- function(x,y){
  0.5/pi/(1+y^2.0)*exp(-0.5*x^2.0*(1.0+y^2))
}

TOwen.gr.1 <- function(x,y){
# -0.5/sqrt(2.0*pi)*(erf(x*y/sqrt(2.0)))*exp(-0.5*x^2.0)
  -0.5/sqrt(2.0*pi)*(2.0*pnorm(x*y)-1.0)*exp(-0.5*x^2.0)
}
 
*/

double TOwenGr1 ( double x, double y )
  
  // the derivative of the T Owen function w.r.t. the second argument
{
  double value;
  
  //value = -0.5/sqrt(2.0*M_PI)*(2.0*pnorm(x*y, 0.0, 1.0, 1.0, 0.0)-1.0)*exp(-0.5*x*x);
  value = -0.5 / sqrt(2.0*M_PI) * erf(x*y/sqrt(2.0)) * exp(-0.5*x*x);
  
  return value;
}

double TOwenGr2 ( double x, double y )
  
  // the derivative of the T Owen function w.r.t. the first argument
{
  double value;
  
  value = 0.5 / M_PI / (1.0+y*y) * exp(-0.5*x*x*(1.0+y*y));
  
  return value;
}

void ll_sn_exp(
    int *Nthreds,
    double *prod, double *y, double *x,
    double *zsv, double *zsk, double *zsu,
    int *n, int *k, int *ksv, int *ksk, int *ksu,
    double *theta, 
    double *lnls, double *lnl){
  
  double lnl0 = 0.0;
  
  // http://jakascorner.com/blog/2016/06/omp-for-reduction.html about `reduction`
  #pragma omp parallel for num_threads(*Nthreds) reduction(+: lnl0)
  for(int i = 0; i < *n; i++){
    // init
    double eps0 = y[i];
    double sv   = 0.0;
    double al0  = 0.0;
    double lam  = 0.0;
    //eps0  <- y - x %*%    theta[1                        :k]
    for(int q = 0; q < *k; q++){
      eps0 -= x[ q**n + i ] * theta[q];
    }
    //sv    <- sqrt(exp( zsv %*% theta[(k + 1)                  :(k + ksv)]))
    for(int q = 0; q < *ksv; q++){
      sv += zsv[ q**n + i ] * theta[*k+ q];
    }
    sv = sqrt(exp(sv));
    // Skewness
    //al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
    for(int q = 0; q < *ksk; q++){
      al0 += zsk[ q**n + i ] * theta[*k + *ksv + q];
    }
    //lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] ))
    for(int q = 0; q < *ksu; q++){
      lam += zsu[ q**n + i ] * theta[*k + *ksv + *ksk + q];
    }
    lam = exp(-lam/2.0); //1.0 / sqrt(exp(lam));
    
    //epsr  = eps0 + sv * sqrt(2.0/M_PI) * al0 / sqrt(1 + pow(al0,2.0));
    double epsr  = eps0 + sv * M_2_SQRTPI * M_SQRT1_2 * al0 / sqrt(1.0 + al0*al0);
    
    //double u1    = (prod[0]*epsr + lam*pow(sv,2.0))/sv;
    double u1    = prod[0]*epsr/sv + lam*sv;
    double b1    = *prod * al0;
    //double b2    = sqrt(1.0 + b1*b1);
    double a1    = -b1 * lam * sv;
    double a2    = a1 / sqrt(1.0 + b1*b1);
    
    //double term1 = -tfn( u1, a2 / u1 );
    //double term2 = -tfn( a2, u1 / a2 );
    //double term3 =  tfn( u1, b1 + a1 / u1 );
    //double term4 =  tfn( a2, b1 + u1 * (1.0 + b1*b1) / a1 );    
    
    double term1 = -t( u1, a2 / u1 );
    double term2 = -t( a2, u1 / a2 );
    double term3 =  t( u1, b1 + a1 / u1 );
    double term4 =  t( a2, b1 + u1 * (1.0 + b1*b1) / a1 );
    
    //double term1 = -tha( u1, 1.0, a2 / u1, 1.0 );
    //double term2 = -tha( a2, 1.0, u1 / a2, 1.0 );
    //double term3 =  tha( u1, 1.0, b1 + a1 / u1, 1.0 );
    //double term4 =  tha( a2, 1.0, b1 + u1 * (1.0 + b1*b1) / a1, 1.0 );
    
    //double term5 = pnorm(a2, 0.0, 1.0, 1.0, 0.0) * pnorm(-u1, 0.0, 1.0, 1.0, 0.0);
    double term5 = pnormstd(a2) * pnormstd(-u1);
    double t1_5  = term1 + term2 + term3 + term4 + term5;
    //Rprintf("i=%2i, term5=%4.8f, t1_5=%4.8e \n", i, term5, t1_5);
    //if(t1_5 <= 0.0){
      //t1_5 = sqrt(DBL_MIN);
      //double lt1_5 = log(t1_5);
      //Rprintf("Was negative; i=%2i, lt1_5=%4.8f, t1_5=%4.8f \n", i, lt1_5, t1_5);
    //}
    

    //lnls[i] = log( t1_5 ) + log( 2.0*lam ) + *prod*epsr*lam + pow(lam*sv,2.0)/2.0;
    //lnls[i] = log( t1_5 * 2.0 * lam ) + *prod*epsr*lam + 0.5*pow(lam*sv, 2.0);
    lnls[i] = log( t1_5 * 2.0*lam )  + *prod*epsr*lam + 0.5 * lam*sv*lam*sv;
    //lnls[i] = log( t1_5 * 2.0 * lam * exp(*prod*epsr*lam) * exp( 0.5 * lam*sv*lam*sv) );
    lnl0 +=lnls[i];
    //Rprintf(" End of i=%2i, t1_5=%4.8e, lnls[i]=%4.8f, lnl0=%4.8f \n", i, t1_5, lnls[i], lnl0);
  }
  *lnl = lnl0;
}

void sn_exp_ll_grad(
    int *Nthreds,
    double *prod, double *y, double *x,
    double *zsv, double *zsk, double *zsu,
    int *n, int *k, int *ksv, int *ksk, int *ksu,
    double *theta, 
    double *lnls, double *lnl, 
    double *grad, double *grad2, double *t1_5){
  
  double lnl0 = 0.0;
  // Allocate memory
  //double * eb   = (double *) malloc (*k * sizeof(double));
  
  double * SVUv = (double *) malloc (*ksv * sizeof(double));
  double * SVUu = (double *) malloc (*ksu * sizeof(double));
  double * As   = (double *) malloc (*ksk * sizeof(double));
  
  double * av   = (double *) malloc (*ksv * sizeof(double));
  double * au   = (double *) malloc (*ksu * sizeof(double));
  //double * alps = (double *) malloc (*ksk * sizeof(double));
  double * as   = (double *) malloc (*ksk * sizeof(double));
  
  double * d_u1_b    = (double *) malloc (*k * sizeof(double));
  double * d_a2_b    = (double *) malloc (*k * sizeof(double));
  double * d_a2u1_b  = (double *) malloc (*k * sizeof(double));
  double * d_u1a2_b  = (double *) malloc (*k * sizeof(double));
  double * d_bau1_b  = (double *) malloc (*k * sizeof(double));
  double * d_bu1a_b  = (double *) malloc (*k * sizeof(double));
  double * d_rest_b  = (double *) malloc (*k * sizeof(double));
  double * d_term1_b = (double *) malloc (*k * sizeof(double));
  double * d_term2_b = (double *) malloc (*k * sizeof(double));
  double * d_term3_b = (double *) malloc (*k * sizeof(double));
  double * d_term4_b = (double *) malloc (*k * sizeof(double));
  double * d_term5_b = (double *) malloc (*k * sizeof(double));

  double * d_u1_v    = (double *) malloc (*ksv * sizeof(double));
  double * d_a2_v    = (double *) malloc (*ksv * sizeof(double));
  double * d_a2u1_v  = (double *) malloc (*ksv * sizeof(double));
  double * d_u1a2_v  = (double *) malloc (*ksv * sizeof(double));
  double * d_bau1_v  = (double *) malloc (*ksv * sizeof(double));
  double * d_bu1a_v  = (double *) malloc (*ksv * sizeof(double));
  double * d_rest_v  = (double *) malloc (*ksv * sizeof(double));
  double * d_term1_v = (double *) malloc (*ksv * sizeof(double));
  double * d_term2_v = (double *) malloc (*ksv * sizeof(double));
  double * d_term3_v = (double *) malloc (*ksv * sizeof(double));
  double * d_term4_v = (double *) malloc (*ksv * sizeof(double));
  double * d_term5_v = (double *) malloc (*ksv * sizeof(double));
  
  double * d_u1_s    = (double *) malloc (*ksk * sizeof(double));
  double * d_a2_s    = (double *) malloc (*ksk * sizeof(double));
  double * d_a2u1_s  = (double *) malloc (*ksk * sizeof(double));
  double * d_u1a2_s  = (double *) malloc (*ksk * sizeof(double));
  double * d_bau1_s  = (double *) malloc (*ksk * sizeof(double));
  double * d_bu1a_s  = (double *) malloc (*ksk * sizeof(double));
  double * d_rest_s  = (double *) malloc (*ksk * sizeof(double));
  double * d_term1_s = (double *) malloc (*ksk * sizeof(double));
  double * d_term2_s = (double *) malloc (*ksk * sizeof(double));
  double * d_term3_s = (double *) malloc (*ksk * sizeof(double));
  double * d_term4_s = (double *) malloc (*ksk * sizeof(double));
  double * d_term5_s = (double *) malloc (*ksk * sizeof(double));
  
  double * d_u1_u    = (double *) malloc (*ksu * sizeof(double));
  double * d_a2_u    = (double *) malloc (*ksu * sizeof(double));
  double * d_a2u1_u  = (double *) malloc (*ksu * sizeof(double));
  double * d_u1a2_u  = (double *) malloc (*ksu * sizeof(double));
  double * d_bau1_u  = (double *) malloc (*ksu * sizeof(double));
  double * d_bu1a_u  = (double *) malloc (*ksu * sizeof(double));
  double * d_rest_u  = (double *) malloc (*ksu * sizeof(double));
  double * d_term1_u = (double *) malloc (*ksu * sizeof(double));
  double * d_term2_u = (double *) malloc (*ksu * sizeof(double));
  double * d_term3_u = (double *) malloc (*ksu * sizeof(double));
  double * d_term4_u = (double *) malloc (*ksu * sizeof(double));
  double * d_term5_u = (double *) malloc (*ksu * sizeof(double));

  // http://jakascorner.com/blog/2016/06/omp-for-reduction.html about `reduction`
  #pragma omp parallel for num_threads(*Nthreds) reduction(+: lnl0)
  for(int i = 0; i < *n; i++){
    // init
    double eps0 = y[i];
    double sv   = 0.0;
    double alp  = 0.0;
    double lam  = 0.0;
    
    //eps0  <- y - x %*%    theta[1                        :k]
    for(int q = 0; q < *k; q++){
      eps0 -= x[ q**n + i ] * theta[q];
      //eb[q] = -x[ q**n + i ];
    }
    //sv    <- sqrt(exp( zsv %*% theta[(k + 1)                  :(k + ksv)]))
    for(int q = 0; q < *ksv; q++){
      sv += zsv[ q**n + i ] * theta[*k+ q];
    }
    sv = sqrt(exp(sv));
    // Skewness
    //al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
    for(int q = 0; q < *ksk; q++){
      alp += zsk[ q**n + i ] * theta[*k + *ksv + q];
    }
    //lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] ))
    for(int q = 0; q < *ksu; q++){
      lam += zsu[ q**n + i ] * theta[*k + *ksv + *ksk + q];
    }
    lam = exp(-lam/2.0); //1.0 / sqrt(exp(lam));
    double su   = 1.0 / lam;
    
    double SVU = sv * lam;
 
    for(int q = 0; q < *ksv; q++){
      SVUv[q] = zsv[ q**n + i ] * SVU * 0.5;
      av[q]   = -SVUv[q] * *prod * alp;
    }
    for(int q = 0; q < *ksu; q++){
      SVUu[q] = -zsu[ q**n + i ] * SVU * 0.5;
      au[q]   = -SVUu[q] * *prod * alp;
    }
    double A = alp / sqrt(1+alp*alp);
    for(int q = 0; q < *ksk; q++){
      As[q]   = zsk[ q**n + i ] * pow(1.0 + alp*alp, -1.5);
      //alps[q] = zsk[ q**n + i ];
      //as[q]   = -alps[q] * *prod * SVU;
      as[q]   = -zsk[ q**n + i ] * *prod * SVU;
    }
    
    //epsr  = eps0 + sv * sqrt(2.0/M_PI) * al0 / sqrt(1 + pow(al0,2.0));
    double epsr  = eps0 + sv * M_2_SQRTPI * M_SQRT1_2 * alp / sqrt(1.0 + alp*alp);
    
    //double u1    = (prod[0]*epsr + lam*pow(sv,2.0))/sv;
    double u1    = prod[0] * epsr/sv + lam*sv;
    double b     = *prod * alp;
    //double b2    = sqrt(1.0 + b1*b1);
    double a     = -b * lam * sv;
    double a2    = a / sqrt(1.0 + b*b);

    double term1 = -t( u1, a2 / u1 );
    double term2 = -t( a2, u1 / a2 );
    double term3 =  t( u1, b + a / u1 );
    double term4 =  t( a2, b + u1 * (1.0 + b*b) / a );
    
    //double term1 = -tha( u1, 1.0, a2 / u1, 1.0 );
    //double term2 = -tha( a2, 1.0, u1 / a2, 1.0 );
    //double term3 =  tha( u1, 1.0, b + a / u1, 1.0 );
    //double term4 =  tha( a2, 1.0, b + u1 * (1.0 + b*b) / a, 1.0 );

    //double term1 = -tfn( u1, a2 / u1 );
    //double term2 = -tfn( a2, u1 / a2 );
    //double term3 =  tfn( u1, b + a / u1 );
    //double term4 =  tfn( a2, b + u1 * (1.0 + b*b) / a );
    
    //double term5 = pnorm(a2, 0.0, 1.0, 1.0, 0.0) * pnorm(-u1, 0.0, 1.0, 1.0, 0.0);
    double term5 = pnormstd(a2) * pnormstd(-u1);
    t1_5[i]  = term1 + term2 + term3 + term4 + term5;
    //if(t1_5[i] <= 0.0){
    //  t1_5[i] = DBL_MIN;
    //}

    //lnls[i] = log( t1_5[i] * 2.0 * lam * exp(*prod*epsr*lam) * exp( 0.5 * lam*sv*lam*sv) );
    lnls[i] = log( t1_5[i]*2.0*lam) + *prod*epsr*lam + pow(lam*sv,2.0)/2.0;
    
    //log( t1_5[i] * 2.0 * lam) + exp(*prod*epsr*lam) * exp( 0.5 * lam*sv*lam*sv) );
    lnl0 +=lnls[i];
    
    double T_gr_1_u1_a2u1 = TOwenGr1(u1,a2/u1); 
    double T_gr_2_u1_a2u1 = TOwenGr2(u1,a2/u1);
    
    double T_gr_1_a2_u1a2 = TOwenGr1(a2,u1/a2);
    double T_gr_2_a2_u1a2 = TOwenGr2(a2,u1/a2);
    
    double T_gr_1_u1_bau1 = TOwenGr1(u1,b+a/u1);
    double T_gr_2_u1_bau1 = TOwenGr2(u1,b+a/u1);
    
    double T_gr_1_a2_bu1a = TOwenGr1(a2,b+u1*(1.0+b*b)/a);
    double T_gr_2_a2_bu1a = TOwenGr2(a2,b+u1*(1.0+b*b)/a);
    
    double Phis_gr_1 = pnormstd(a2) * dnormstd(-u1);
    double Phis_gr_2 = pnormstd(-u1) * dnormstd(a2);
    //double Phis_gr_1 = pnorm(a2, 0.0, 1.0, 1.0, 0.0) * dnorm(-u1, 0.0, 1.0, 0.0);
    //double Phis_gr_2 = pnorm(-u1, 0.0, 1.0, 1.0, 0.0) * dnorm(a2, 0.0, 1.0, 0.0);
    
    // w.r.t. beta
    for(int q = 0; q < *k; q++){
      d_a2_b[q]    = 0.0;
      //d_u1_b[q]    = eb[q] * *prod / sv;
      d_u1_b[q]    = -x[ q**n + i ] * *prod / sv;
      d_a2u1_b[q]  = -d_u1_b[q] * a2 / u1 / u1;
      d_u1a2_b[q]  = d_u1_b[q] / a2;
      d_bau1_b[q]  = -d_u1_b[q] * a / u1 / u1;
      d_bu1a_b[q]  = d_u1_b[q] * (1.0+b*b) / a;
      //d_rest_b[q]  = eb[q] * *prod/su;
      d_rest_b[q]  = -x[ q**n + i ] * *prod / su;
      d_term1_b[q] = d_u1_b[q] * T_gr_1_u1_a2u1 + d_a2u1_b[q] * T_gr_2_u1_a2u1;
      d_term2_b[q] =                              d_u1a2_b[q] * T_gr_2_a2_u1a2;
      d_term3_b[q] = d_u1_b[q] * T_gr_1_u1_bau1 + d_bau1_b[q] * T_gr_2_u1_bau1;
      d_term4_b[q] =                              d_bu1a_b[q] * T_gr_2_a2_bu1a;
      d_term5_b[q] = -d_u1_b[q] * Phis_gr_1 + d_a2_b[q] * Phis_gr_2;
      // give the result back
      grad[q* *n + i] = d_rest_b[q] + (-d_term1_b[q]-d_term2_b[q]+d_term3_b[q]+d_term4_b[q]+d_term5_b[q]) / t1_5[i];
      //grad[q* *n + i] = d_term4_b[q];
    }
    
    // w.r.t. gamma_v
    for(int q = 0; q < *ksv; q++){
      d_a2_v[q]    = -SVUv[q] * *prod * A;
      d_u1_v[q]    = -zsv[ q**n + i ] * *prod * eps0 * 0.5 / sv + SVUv[q];
      d_a2u1_v[q]  = d_a2_v[q] / u1 - d_u1_v[q] * a2 / u1 / u1;
      d_u1a2_v[q]  = d_u1_v[q] / a2 - d_a2_v[q] * u1 / a2 / a2;
      d_bau1_v[q]  = av[q] / u1     - d_u1_v[q] * a / u1 / u1;
      d_bu1a_v[q]  = d_u1_v[q] * (1.0+b*b) / a - av[q] * u1*(1.0+b*b) / a / a;
      d_rest_v[q]  = SVUv[q] * (*prod*sqrt(2.0/M_PI )*A + SVU);
      d_term1_v[q] = d_u1_v[q] * T_gr_1_u1_a2u1 + d_a2u1_v[q] * T_gr_2_u1_a2u1;
      d_term2_v[q] = d_a2_v[q] * T_gr_1_a2_u1a2 + d_u1a2_v[q] * T_gr_2_a2_u1a2;
      d_term3_v[q] = d_u1_v[q] * T_gr_1_u1_bau1 + d_bau1_v[q] * T_gr_2_u1_bau1;
      d_term4_v[q] = d_a2_v[q] * T_gr_1_a2_bu1a + d_bu1a_v[q] * T_gr_2_a2_bu1a;
      d_term5_v[q] = -d_u1_v[q] * Phis_gr_1 + d_a2_v[q] * Phis_gr_2;
      // give the result back
      grad[(q+*k)* *n + i] = d_rest_v[q] + (-d_term1_v[q]-d_term2_v[q]+d_term3_v[q]+d_term4_v[q]+d_term5_v[q]) / t1_5[i];
      //grad[(q+*k)* *n + i] = d_term4_v[q];
    }
    
    // w.r.t. gamma_s
    for(int q = 0; q < *ksk; q++){
      d_a2_s[q]    = -As[q] * *prod * SVU;
      d_u1_s[q]    = *prod * sqrt(2.0/M_PI) * As[q];
      d_a2u1_s[q]  = d_a2_s[q] / u1 - d_u1_s[q] * a2 / u1 / u1;
      d_u1a2_s[q]  = d_u1_s[q] / a2 - d_a2_s[q] * u1 / a2 / a2;
      //d_bau1_s[q]  = as[q] / u1     - d_u1_s[q] * a / u1 / u1 + *prod * alps[q];
      d_bau1_s[q]  = as[q] / u1     - d_u1_s[q] * a / u1 / u1 + *prod * zsk[ q**n + i ];
      //d_bu1a_s[q]  = d_u1_s[q] * (1.0+b*b) / a - as[q] * u1/SVU/SVU*(1.0/b/b-1.0)+ *prod * alps[q];
      d_bu1a_s[q]  = d_u1_s[q] * (1.0+b*b) / a - as[q] * u1/SVU/SVU*(1.0/b/b-1.0)+ *prod * zsk[ q**n + i ];
      d_rest_s[q]  = As[q] * *prod*sqrt(2.0/M_PI)*SVU;
      d_term1_s[q] = d_u1_s[q] * T_gr_1_u1_a2u1 + d_a2u1_s[q] * T_gr_2_u1_a2u1;
      d_term2_s[q] = d_a2_s[q] * T_gr_1_a2_u1a2 + d_u1a2_s[q] * T_gr_2_a2_u1a2;
      d_term3_s[q] = d_u1_s[q] * T_gr_1_u1_bau1 + d_bau1_s[q] * T_gr_2_u1_bau1;
      d_term4_s[q] = d_a2_s[q] * T_gr_1_a2_bu1a + d_bu1a_s[q] * T_gr_2_a2_bu1a;
      d_term5_s[q] = -d_u1_s[q] * Phis_gr_1 + d_a2_s[q] * Phis_gr_2;
      // give the result back
      grad[(q+*k+*ksv)* *n + i] = d_rest_s[q] + (-d_term1_s[q]-d_term2_s[q]+d_term3_s[q]+d_term4_s[q]+d_term5_s[q]) / t1_5[i];
      //grad[(q+*k+*ksv)* *n + i] = d_term4_s[q];
    }
    
    // w.r.t. gamma_u
    for(int q = 0; q < *ksu; q++){
      d_a2_u[q]    = -SVUu[q] * *prod * A;
      d_u1_u[q]    = SVUu[q];
      d_a2u1_u[q]  = d_a2_u[q] / u1 - d_u1_u[q] * a2 / u1 / u1;
      d_u1a2_u[q]  = d_u1_u[q] / a2 - d_a2_u[q] * u1 / a2 / a2;
      d_bau1_u[q]  = au[q] / u1     - d_u1_u[q] * a / u1 / u1;
      d_bu1a_u[q]  = d_u1_u[q] * (1.0+b*b) / a - au[q] * u1*(1.0+b*b) / a / a;
      d_rest_u[q]  = SVUu[q] * (*prod*sqrt(2.0/M_PI)*A + SVU) + zsu[ q**n + i ] * (*prod * eps0 * -0.5 / su - 0.5);
      d_term1_u[q] = d_u1_u[q] * T_gr_1_u1_a2u1 + d_a2u1_u[q] * T_gr_2_u1_a2u1;
      d_term2_u[q] = d_a2_u[q] * T_gr_1_a2_u1a2 + d_u1a2_u[q] * T_gr_2_a2_u1a2;
      d_term3_u[q] = d_u1_u[q] * T_gr_1_u1_bau1 + d_bau1_u[q] * T_gr_2_u1_bau1;
      d_term4_u[q] = d_a2_u[q] * T_gr_1_a2_bu1a + d_bu1a_u[q] * T_gr_2_a2_bu1a;
      d_term5_u[q] = -d_u1_u[q] * Phis_gr_1 + d_a2_u[q] * Phis_gr_2;
      // give the result back
      grad[(q+*k+*ksv+*ksk)* *n + i] = d_rest_u[q] + (-d_term1_u[q]-d_term2_u[q]+d_term3_u[q]+d_term4_u[q]+d_term5_u[q]) / t1_5[i];
      //grad[(q+*k+*ksv+*ksk)* *n + i] = d_term4_u[q];
    }
    
  }
  
  for(int i = 0; i < *n; i++){
    for(int q = 0; q < *k+*ksv+*ksk+*ksu; q++){
      if( i == 0){
        grad2[q] = 0.0;
      }
      grad2[q] += grad[q* *n + i];
    }    
  }
  
  *lnl = lnl0;
  
  //free(eb);
  
  free(SVUv);
  free(SVUu);
  free(As);
  
  free(av);
  free(au);
  //free(alps);
  free(as);
  
  free(d_u1_b);
  free(d_a2_b);
  free(d_a2u1_b);
  free(d_u1a2_b);
  free(d_bau1_b);
  free(d_bu1a_b);
  free(d_rest_b);
  free(d_term1_b);
  free(d_term2_b);
  free(d_term3_b);
  free(d_term4_b);
  free(d_term5_b);
  free(d_u1_v);
  free(d_a2_v);
  free(d_a2u1_v);
  free(d_u1a2_v);
  free(d_bau1_v);
  free(d_bu1a_v);
  free(d_rest_v);
  free(d_term1_v);
  free(d_term2_v);
  free(d_term3_v);
  free(d_term4_v);
  free(d_term5_v);
  free(d_u1_s);
  free(d_a2_s);
  free(d_a2u1_s);
  free(d_u1a2_s);
  free(d_bau1_s);
  free(d_bu1a_s);
  free(d_rest_s);
  free(d_term1_s);
  free(d_term2_s);
  free(d_term3_s);
  free(d_term4_s);
  free(d_term5_s);
  free(d_u1_u);
  free(d_a2_u);
  free(d_a2u1_u);
  free(d_u1a2_u);
  free(d_bau1_u);
  free(d_bu1a_u);
  free(d_rest_u);
  free(d_term1_u);
  free(d_term2_u);
  free(d_term3_u);
  free(d_term4_u);
  free(d_term5_u);
  
}

void u_sn_exp(
    int *Nthreds,
    double *prod, double *y, double *x,
    double *zsv, double *zsk, double *zsu,
    int *n, int *k, int *ksv, int *ksk, int *ksu,
    double *theta, 
    double *u, double *eps01, double *sv1, double *al01, double *lam1){
  
  
  
  #pragma omp parallel for num_threads(*Nthreds)
  for(int i = 0; i < *n; i++){
    // init
    double eps0 = y[i];
    double sv   = 0.0;
    double al0  = 0.0;
    double lam  = 0.0;
    //eps0  <- y - x %*%    theta[1                        :k]
    for(int q = 0; q < *k; q++){
      eps0 -= x[ q**n + i ] * theta[q];
    }
    eps01[i] = eps0;
    //sv    <- sqrt(exp( zsv %*% theta[(k + 1)                  :(k + ksv)]))
    for(int q = 0; q < *ksv; q++){
      sv += zsv[ q**n + i ] * theta[*k+ q];
    }
    sv = sqrt(exp(sv));
    sv1[i] = sv;
    // Skewness
    //al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
    for(int q = 0; q < *ksk; q++){
      al0 += zsk[ q**n + i ] * theta[*k + *ksv + q];
    }
    al01[i] = al0;
    //lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] ))
    for(int q = 0; q < *ksu; q++){
      lam += zsu[ q**n + i ] * theta[*k + *ksv + *ksk + q];
    }
    lam = exp(-lam/2.0); //1.0 / sqrt(exp(lam));
    lam1[i] = lam;
    
    //epsr  = eps0 + sv * sqrt(2.0/M_PI) * al0 / sqrt(1 + pow(al0,2.0));
    double epsr  = eps0 + sv * M_2_SQRTPI * M_SQRT1_2 * al0 / sqrt(1.0 + pow(al0,2.0));
    
    //double u1    = (prod[0]*epsr + lam*pow(sv,2.0))/sv;
    double u1    = prod[0]*epsr/sv + lam*sv;
    double b1    = *prod * al0;
    double b2    = sqrt(1.0 + b1*b1);
    double a1    = -b1 * lam * sv;
    double a2    = a1 / b2;
    
    //double term1 = -tfn( u1, a2 / u1 );
    //double term2 = -tfn( a2, u1 / a2 );
    //double term3 =  tfn( u1, b1 + a1 / u1 );
    //double term4 =  tfn( a2, b1 + u1 * (1.0 + b1*b1) / a1 );    
    
    double term1 = -t( u1, a2 / u1 );
    double term2 = -t( a2, u1 / a2 );
    double term3 =  t( u1, b1 + a1 / u1 );
    double term4 =  t( a2, b1 + u1 * (1.0 + b1*b1) / a1 );
    
    //double term1 = -tha( u1, 1.0, a2 / u1, 1.0 );
    //double term2 = -tha( a2, 1.0, u1 / a2, 1.0 );
    //double term3 =  tha( u1, 1.0, b1 + a1 / u1, 1.0 );
    //double term4 =  tha( a2, 1.0, b1 + u1 * (1.0 + b1*b1) / a1, 1.0 );
    
    //double term5 = pnorm(a2, 0.0, 1.0, 1.0, 0.0) * pnorm(-u1, 0.0, 1.0, 1.0, 0.0);
    double term5 = pnormstd(a2) * pnormstd(-u1);
    double t122  = term1 + term2 + term3 + term4 + term5;
    //double t135  = b1/b2 * dnormstd(a2) * pnormstd(-u1*b2 - b1*a2) + 
    //  dnormstd(u1)*pnormstd(a1+b1*u1);
    double t135  = b1/b2 * dnorm(a2, 0.0, 1.0, 0.0) * pnorm(-u1*b2 - b1*a2, 0.0, 1.0, 1.0, 0.0) + 
      dnorm(u1, 0.0, 1.0, 0.0) * pnorm(a1+b1*u1, 0.0, 1.0, 1.0, 0.0);
    
    u[i] = sv * ( t135 / t122 - u1);
    if(u[i] < 0){
      u[i] = 0;
    }
  }
}

void ll_sn_tn(
    int *Nthreds,
    double *prod, double *y, double *x,
    double *zsv, double *zsk, double *zsu, double *zmu,
    int *n, int *k, int *ksv, int *ksk, int *ksu, int *kmu,
    double *theta, double *tn,
    double *lnls, double *lnl){
  
  double lnl0 = 0.0;
  
  #pragma omp parallel for num_threads(*Nthreds) reduction(+: lnl0)
  for(int i = 0; i < *n; i++){
    double eps0 = y[i];
    double sv2 = 0.0;
    double al0 = 0.0;
    double su2 = 0.0;
    double mu0 = 0.0;
    //eps0  <- y - x %*%    theta[1                        :k]
    for(int q = 0; q < *k; q++){
      eps0 -= x[ q* *n + i ] * theta[q];
    }
    //sv    <- sqrt(exp( zsv %*% theta[(k + 1)                  :(k + ksv)]))
    for(int q = 0; q < *ksv; q++){
      sv2 += zsv[ q* *n + i ] * theta[*k+ q];
    }
    sv2 = exp(sv2);
    // Skewness
    //al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
    for(int q = 0; q < *ksk; q++){
      al0 += zsk[ q* *n + i ] * theta[*k + *ksv + q];
    }
    //lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] ))
    for(int q = 0; q < *ksu; q++){
      su2 += zsu[ q* *n + i ] * theta[*k + *ksv + *ksk + q];
    }
    su2 = exp(su2);
    // mu
    if(*tn == 1.0){
      for(int q = 0; q < *kmu; q++){
        mu0 += zsu[ q* *n + i ] * theta[*k + *ksv + *ksk + *ksu + q];
      }
    }
    //Rprintf("i = %2i, *muIsZero = %4.8e, mu0 = %4.8e \n", i, *muIsZero, mu0);
    
    double sv    = sqrt(sv2);
    double su    = sqrt(su2);
    double sig   = sqrt(sv2 + su2);
    double sstar = sv*su/sig;
    //epsr  = eps0 + sv * sqrt(2.0/M_PI) * al0 / sqrt(1 + pow(al0,2.0));
    double epsr  = eps0 + sv * M_2_SQRTPI * M_SQRT1_2 * al0 / sqrt(1.0 + al0*al0);
    double mu1   = (mu0*sv2 - epsr* *prod*su2) / pow(sig,2.0);
    double b1    = al0 / sv * *prod * sstar ;
    double a1    = al0 / sv * (epsr + *prod * mu1);
    double a2    = a1 / sqrt(1.0 + pow(b1,2.0));
    double u1    = -mu1 / sstar;
    
    double term1 = -t( u1, a2 / u1 );
    double term2 = -t( a2, u1 / a2 );
    double term3 = t( u1, b1 + a1 / u1 );
    double term4 = t( a2, b1 + u1 * (1.0 + pow(b1,2.0)) / a1 );
    double term5 = pnormstd(a2) * pnormstd(-u1);
    
    double t1_5  = term1 + term2 + term3 + term4 + term5;
    
    lnls[i] = log( t1_5*2.0/sig) - log(pnormstd( mu0/su )) + log(dnormstd( (mu0 + *prod*epsr)/sig ));
    lnl0 +=lnls[i];
  }
  *lnl = lnl0;
}

// fun0:	ordered probit
// N:		y0, xb0
// K:		yl0
// K-1: 	k0
// 1:		v0
// rezval1: N

void sf_oprobit_ll(int *Nthreds,
                   double *prod, double *y, double *y0, 
                   double *x, double *zsu,
                   int *n, int *k, int *ksu, int *kmu,
                   double *theta, 
                   double *su0, double *xb0, 
                   double *lnlY, double *lnl)
{
  //for(int j = 0; j < *k + *ksu + *kmu - 2; j++){
  //  Rprintf("j = %1i, k_j = %4.4f, k_j+1 = %4.4f\n", j, theta[j], theta[j + 1]);
  //  }
  // check if k0 is sorted
  for (int j = 0; j < *kmu - 2; j++)
  {
    if (theta[*k + *ksu + j] > theta[*k + *ksu + j + 1])
    {
      Rprintf("j = %1i, k_j = %4.4f, k_j+1 = %4.4f\n", j, theta[*k + *ksu + j], theta[*k + *ksu + j + 1]);
      Rprintf("The supplied vector m0 is not sorted\n");
      return;
    }
  }
  
  
  double lnl0 = 0.0;
  
  // http://jakascorner.com/blog/2016/06/omp-for-reduction.html about `reduction`
#pragma omp parallel for num_threads(*Nthreds) reduction(+: lnl0)
  for (int i = 0; i < *n; i++)
  {
    // We have K+1 possible values for y, 
    // we need to figure which this current y[i] is
    int yi;
    for (int j = 0; j < *kmu; j++)
    {
      //if(i < 20){
      //  Rprintf("i = %1i, j = %1i, y[i] = %4.4f, y0[j] = %4.4f\n", 
      //          i, j, y[i], y0[j]);
      //}
      if (y[i] == y0[j])
      {
        yi = j;
        //if(i < 20){
        //  Rprintf("CHosen: i = %1i, j = %1i, yi = %1i, y[i] = %4.4f, y0[j] = %4.4f\n", 
        //          i, j, yi, y[i], y0[j]);
        //}
        break;
      }
    }
    // Thus, current y[i] is yi's order;
    // xb  <- x %*% theta[1:k]
    double xb = 0.0;
    for(int q = 0; q < *k; q++){
      xb += x[ q**n + i ] * theta[q];
    }
    xb0[i] = xb;
    // su2 <- exp( zsu %*% theta[(k + 1) :(k + ksu)])
    double lnsu2 = 0.0;
    for(int q = 0; q < *ksu; q++){
      lnsu2 += zsu[ q**n + i ] * theta[*k + q];
    }
    double su   = sqrt(exp(lnsu2));
    double ze = 1.0 / sqrt(1.0 + exp(lnsu2));
    su0[i] = su;
    //if(i < 20){
    //  Rprintf("i = %1i, su = %4.4f, xb = %4.4f, \n", 
    //        i, su, xb);
    //}
    //if(i < 20){
    //  Rprintf("Before the loop: i = %1i, yi = %1i, y[i] = %4.4f\n", 
    //          i, yi, y[i]);
    //}
    // If it the first one
    if(yi == 0)	{
      double h = ( theta[*k + *ksu] - xb ) * ze;
      double a = prod[0] * su;
      lnlY[i] = log( pnormstd(h) + 2.0 * t(h, a) );		  
      // lnlY[i] = log( pnorm(am, 0.0, 1.0, 1.0, 0.0) + 2.0 * t(am, -prod[0] * su) );
      // if(i < 2000){
      // Rprintf("i = %1i, first, yi = %1i, su = %4.4f, h = %4.4f, pnormstd(h) = %4.4f, t(h, a) = %4.4e, lnls[i]+1 = %4.4f\n",
      // i, yi, su, h, pnormstd(h), t(h, a), lnlY[i]);
      // }
    } else
      // If it is the last one
      if(yi == *kmu - 1) {
        double h = ( theta[*k + *ksu + *kmu - 2] - xb ) * ze;
        double a = prod[0] * su;
        lnlY[i] = log( pnormstd(-h) + 2.0 * t(-h, -a) );      
        //lnlY[i] = log( pnorm(am, 0.0, 1.0, 1.0, 0.0) + 2.0 * t(am, prod[0] * su) );
        //if(i < 20){
        //  Rprintf("i = %1i, last, yi = %1i, su = %4.4f, am = %4.4f, pnormstd(am) = %4.4f, t(am, -prod[0] * su) = %4.4e, lnls[i]+1 = %4.4f\n", 
        //          i, yi, su, am, pnormstd(am), t(am, -prod[0] * su), theta[*k + *ksu + *kmu - 2]);
        //}
      } else {
        // The rest
        double h1 = ( theta[*k + *ksu + yi] - xb ) * ze;
        double h2 = ( theta[*k + *ksu + yi - 1] - xb ) * ze;
        double a = prod[0] * su;
        lnlY[i] = log( pnormstd(h1) + 2.0 * t(h1, a) -
          pnormstd(h2) - 2.0 * t(h2, a) );
        //lnlY[i] = log( pnorm(am, 0.0, 1.0, 1.0, 0.0) + 2.0 * t(am, -prod[0] * su) -
        //  pnorm(am1, 0.0, 1.0, 1.0, 0.0) - 2.0 * t(am1, -prod[0] * su) );
        //if(i < 20){
        //  Rprintf("i = %1i, middle, yi = %1i, su = %4.4f, am = %4.4f, pnormstd(am) = %4.4e, t(am, -prod[0] * su) = %4.4e, lnls[i]+1 = %4.4f\n", 
        //          i, yi, su, am, pnormstd(am), t(am, -prod[0] * su), theta[*k + *ksu + yi]);
        //}
      }
      lnl0 += lnlY[i];
  }
  *lnl = lnl0;
}


// fun0:	ordered probit
// N:		y0, xb0
// K:		yl0
// K-1: 	k0
// 1:		v0
// rezval1: N

void sf_oprobit_gr_i(int *Nthreds,
                     double *prod, double *y, double *y0, 
                     double *x, double *zsu,
                     int *n, int *k, int *ksu, int *kmu,
                     double *theta, 
                     double *lnlY, double *lnl,
                     double *gradY, double *grad)
{
  
  //for(int j = 0; j < *k + *ksu + *kmu - 2; j++){
  //  Rprintf("j = %1i, k_j = %4.4f, k_j+1 = %4.4f\n", j, theta[j], theta[j + 1]);
  //  }
  // check if k0 is sorted
  for (int j = 0; j < *kmu - 2; j++)
  {
    if (theta[*k + *ksu + j] > theta[*k + *ksu + j + 1])
    {
      //Rprintf("j = %1i, k_j = %4.4f, k_j+1 = %4.4f\n", j, theta[*k + *ksu + j], theta[*k + *ksu + j + 1]);
      Rprintf("The supplied vector m0 is not sorted\n");
      return;
    }
  }
  
  
  double lnl0 = 0.0;
  
  // http://jakascorner.com/blog/2016/06/omp-for-reduction.html about `reduction`
#pragma omp parallel for num_threads(*Nthreds) reduction(+: lnl0)
  for (int i = 0; i < *n; i++)
  {
    // We have K+1 possible values for y, 
    // we need to figure which this current y[i] is
    int yi;
    for (int j = 0; j < *kmu; j++)
    {
      //Rprintf("i = %1i, start of i, yi = %1i\n",  i, yi);
      
      //if(i < 20){
      //  Rprintf("i = %1i, j = %1i, y[i] = %4.4f, y0[j] = %4.4f\n", 
      //          i, j, y[i], y0[j]);
      //}
      if (y[i] == y0[j])
      {
        yi = j;
        //if(i < 20){
        //  Rprintf("CHosen: i = %1i, j = %1i, yi = %1i, y[i] = %4.4f, y0[j] = %4.4f\n", 
        //          i, j, yi, y[i], y0[j]);
        //}
        break;
      }
    }
    // Thus, current y[i] is yi's order;
    // xb  <- x %*% theta[1:k]
    
    double xb = 0.0;
    for(int q = 0; q < *k; q++){
      xb += x[ q**n + i ] * theta[q];
    }
    // xb0[i] = xb;
    // su2 <- exp( zsu %*% theta[(k + 1) :(k + ksu)])
    double lnsu2 = 0.0;
    for(int q = 0; q < *ksu; q++){
      lnsu2 += zsu[ q**n + i ] * theta[*k + q];
    }
    double su2  = exp(lnsu2);
    double su   = sqrt(su2);
    double ze = 1.0 / sqrt(1.0 + exp(lnsu2));
    double ze2 = ze*ze;
    //if(i < 20){
    //  Rprintf("i = %1i, su = %4.4f, xb = %4.4f, \n", 
    //        i, su, xb);
    //}
    //if(i < 20){
    //  Rprintf("Before the loop: i = %1i, yi = %1i, y[i] = %4.4f\n", 
    //          i, yi, y[i]);
    //}
    // If it the first one
    if(yi == 0)	{
      //double am = ( theta[*k + *ksu] - xb ) / su2s;
      //lnlY[i] = log( pnormstd(am) + 2.0 * t(am, -prod[0] * su) );		  
      //lnlY[i] = log( pnorm(am, 0.0, 1.0, 1.0, 0.0) + 2.0 * t(am, -prod[0] * su) );
      
      double h = ( theta[*k + *ksu] - xb ) * ze;
      double b  = prod[0] * su;
      double Lm = pnormstd(h) + 2.0 * t(h, b);
      double fh =  dnormstd(h) * pnormstd(-b*h);
      double D1m = 2.0 * ze * fh;
      double D2m = ze2 * (su2 * h * fh - exp(-h * h / (2.0 * ze2)) * M_1_PI / 2.0 * b );
      // if(i < 500){
      //   Rprintf("i = %1i, first, yi = %1i, su = %4.4f, h = %4.4f, pnormstd(h) = %4.4f, t(am, -prod[0] * su) = %4.4e, lnls[i]+1 = %4.4f\n",
      //           i, yi, su, h, pnormstd(h), t(h, a), theta[*k + *ksu]);
      // }
      lnlY[i] = log( Lm ) ;
      // w.r.t. beta
      for(int q = 0; q < *k; q++){
        gradY[q * *n + i] = -x[ q**n + i ] * D1m / Lm ;
      }
      // w.r.t. gamma
      for(int q = 0; q < *ksu; q++){
        gradY[(q+*k) * *n + i] = -zsu[ q**n + i ] * D2m / Lm ;
      }      
      // w.r.t. mu
      for(int q = 0; q < *kmu-1; q++){
        gradY[(q+*k+*ksu) * *n + i] = 0 ;
      }
      gradY[(yi+*k+*ksu) * *n + i] = D1m / Lm ;
    } else
      // If it is the last one
      if(yi == *kmu - 1) {
        //double am = -( theta[*k + *ksu + *kmu - 2] - xb ) / su2s;
        //lnlY[i] = log( pnormstd(am) + 2.0 * t(am, prod[0] * su) );      
        //lnlY[i] = log( pnorm(am, 0.0, 1.0, 1.0, 0.0) + 2.0 * t(am, prod[0] * su) );
        //if(i < 500){
        //  Rprintf("i = %1i, last, yi = %1i, su = %4.4f, am = %4.4f, pnormstd(am) = %4.4f, t(am, -prod[0] * su) = %4.4e, lnls[i]+1 = %4.4f\n", 
        //          i, yi, su, am, pnormstd(am), t(am, -prod[0] * su), theta[*k + *ksu + *kmu - 2]);
        //}
        
        double h = ( theta[*k + *ksu + *kmu - 2] - xb ) * ze;
        double b  = prod[0] * su;
        double Lm = pnormstd(-h) + 2.0 * t(-h, -b);
        double fh =  dnormstd(-h) * pnormstd(-b*h) ;
        double D1m = -2.0 * ze * fh;
        double D2m = -ze2 * ( su2 * h * fh - exp(-h * h / (2.0 * ze2)) * M_1_PI / 2.0 * b );
        // Rprintf("i = %1i, 123456, h = %4.4f, b = %4.4f, Lm = %4.4f, fh = %4.4e, D1m = %4.4e, D2m = %4.4e, p = %4.0f\n",
                // i, h, a, Lm, fh, D1m, D2m, prod[0]);
        // if(i < 500){
          // Rprintf("i = %1i, 123456, yi = %1i, su = %4.4f, h = %4.4f, pnormstd(h) = %4.4f, t(h, -a) = %4.4e, D1m = %4.4e, D2m = %4.4e, p = %4.04f\n",
          //         i, yi, su, h, pnormstd(h), t(-h, -a), D1m, D2m, prod[0]);
        // }
        lnlY[i] = log( Lm ) ;
        // w.r.t. beta
        for(int q = 0; q < *k; q++){
          //if(i < 500){
          //  Rprintf("i = %1i, beta, q * *n + i = %1i, x[ q**n + i ] * f0 / su2s = %4.4f\n", 
          //         i, q * *n + i, x[ q**n + i ] * f0 / su2s);
          //}
          gradY[q * *n + i] = -x[ q**n + i ] * D1m / Lm  ;
        }
        // w.r.t. gamma
        for(int q = 0; q < *ksu; q++){
          gradY[(q+*k) * *n + i] = -zsu[ q**n + i ] * D2m / Lm ;
        }      
        // w.r.t. mu
        for(int q = 0; q < *kmu-1; q++){
          gradY[(q+*k+*ksu) * *n + i] = 0 ;
        }
        gradY[(yi-1+*k+*ksu) * *n + i] = D1m / Lm ;
      } else {
        // The rest
        //double am = ( theta[*k + *ksu + yi] - xb ) / su2s;
        //double am1 = ( theta[*k + *ksu + yi - 1] - xb ) / su2s;
        //lnlY[i] = log( pnormstd(am) + 2.0 * t(am, -prod[0] * su) -
        //  pnormstd(am1) - 2.0 * t(am1, -prod[0] * su) );
        //lnlY[i] = log( pnorm(am, 0.0, 1.0, 1.0, 0.0) + 2.0 * t(am, -prod[0] * su) -
        //  pnorm(am1, 0.0, 1.0, 1.0, 0.0) - 2.0 * t(am1, -prod[0] * su) );
        // if(i < 500){
        //  Rprintf("i = %1i, middle, yi = %1i, su = %4.4f, am = %4.4f, pnormstd(am) = %4.4e, t(am, -prod[0] * su) = %4.4e, lnls[i]+1 = %4.4f\n", 
        //          i, yi, su, am, pnormstd(am), t(am, -prod[0] * su), theta[*k + *ksu + yi]);
        // }
        double h1 = ( theta[*k + *ksu + yi] - xb ) * ze;
        double h2 = ( theta[*k + *ksu + yi - 1] - xb ) * ze;
        double a  = prod[0] * su;
        // double a  = su;
        double fh1 = dnormstd(h1) * pnormstd(-a*h1) ;
        double fh2 = dnormstd(h2) * pnormstd(-a*h2) ;
        double D1m = 2.0 * ze * (fh1 - fh2);
        double D2m = ze2 * (
          su2 * fh1 * h1 - exp(-h1 * h1 / (2.0 * ze2)) * M_1_PI / 2.0 * a -
            su2 * fh2 * h2 + exp(-h2 * h2 / (2.0 * ze2)) * M_1_PI / 2.0 * a 
        );
        // double D2m = ze2 * (su2 * h * fh - exp(-h * h / (2.0 * ze2)) * M_1_PI / 2.0 * a );
        
        // double fu = 0.5 * su2 / (1.0 + su2) * f0 * h0 - 
        //   exp(-h0 * h0 * (1.0 + a*a) / 2.0) * M_1_PI / 2.0 / (1.0 + a*a) * a -
        //   0.5 * su2 / (1.0 + su2) * f1 * h1 + 
        //   exp(-h1 * h1 * (1.0 + a*a) / 2.0) * M_1_PI / 2.0 / (1.0 + a*a) * a;
        double Lm = pnormstd(h1) + 2.0 * t(h1, a) - pnormstd(h2) - 2.0 * t(h2, a) ;
        lnlY[i] = log( Lm ) ;
        // w.r.t. beta
        for(int q = 0; q < *k; q++){
          gradY[q * *n + i] = -x[ q**n + i ] * D1m / Lm ;
          //if(i < 500){
          //  Rprintf("i = %1i, beta, q * *n + i = %1i, -x[ q**n + i ] * (f0-f1) / su2s = %4.4f\n", 
          //        i, q * *n + i, -x[ q**n + i ] * (f0-f1) / su2s);
          //}
        }
        // w.r.t. gamma
        for(int q = 0; q < *ksu; q++){
          gradY[(q+*k) * *n + i] = -zsu[ q**n + i ] * D2m / Lm ;
          //if(i < 500){
          //  Rprintf("i = %1i, gamma, (q+*k) * *n + i = %1i, q * *n + i = %1i, -x[ q**n + i ] * (f0-f1) / su2s = %4.4f\n", 
          //          i, (q+*k) * *n + i, q * *n + i, -zsu[ q**n + i ] * fu);
          //}
        }      
        // w.r.t. mu
        for(int q = 0; q < *kmu-1; q++){
          //if(i < 500){
          //  Rprintf("i = %1i, mu, q = %1i, (*k+*ksu+*kmu-1)**n = %1i, (q+*k+*ksu) * *n + i = %1i, gradY[(q+*k+*ksu) * *n + i] = %4.4f\n", 
          //          i, q, (*k+*ksu+*kmu-1)**n, (q+*k+*ksu) * *n + i, q * *n + i, gradY[(q+*k+*ksu) * *n + i]);
          //}
          gradY[(q+*k+*ksu) * *n + i] = 0.0 ;
        }
        if(i < 500){
          //Rprintf("i = %1i, mu, (yi+*k+*ksu) * *n + i = %1i, (yi-1+*k+*ksu) * *n + i = %1i, -f1 / su2s = %4.4f\n", 
          //        i, (yi+*k+*ksu) * *n + i, (yi-1+*k+*ksu) * *n + i, -f1 / su2s);
        }
        // gradY[(yi+*k+*ksu) * *n + i] =  D1m / Lm ;
        gradY[(yi+*k+*ksu) * *n + i] =  2.0 * ze * fh1 / Lm ;
        gradY[(yi-1+*k+*ksu) * *n + i] = -2.0 * ze * fh2 / Lm ; // this is from previous mu
      }
      lnl0 += lnlY[i];
      //if(i < 500){
      //Rprintf("i = %1i, end of i, yi = %1i, su = %4.4f\n",  i, yi, su);
      //}
  } // end of the `n` loop
  // return the log-likelihood
  *lnl = lnl0;
  // collect the gradient
  for(int i = 0; i < *n; i++){
    //if(i < 50){
    //  Rprintf("i = %1i, collect gradient start of i\n",  i);
    //}
    for(int q = 0; q < *k+*ksu+*kmu-1; q++){
      if( i == 0 ){
        grad[q] = 0.0;
        //Rprintf("i = %1i, q = %1i, collect gradient end of i\n",  i, q);
        //Rprintf("i = %1i, q = %1i, gradY[q * *n + i] = %4.4f, grad[q] = %4.4f\n",  i, q, gradY[q * *n + i], grad[q] ) ;
      }
      grad[q] += gradY[q * *n + i];
    }
    //Rprintf("i = %1i, collect gradient end of i\n",  i);
  } // end of the collecting
  //Rprintf("All is done\n");
}


void sf_oprobit_u(int *Nthreds,
                  double *prod, double *y, double *y0, 
                  double *x, double *zsu,
                  int *n, int *k, int *ksu, int *kmu,
                  double *theta, 
                  double *uY)
{
  //for(int j = 0; j < *k + *ksu + *kmu - 2; j++){
  //  Rprintf("j = %1i, k_j = %4.4f, k_j+1 = %4.4f\n", j, theta[j], theta[j + 1]);
  //  }
  // check if k0 is sorted
  for (int j = 0; j < *kmu - 2; j++)
  {
    if (theta[*k + *ksu + j] > theta[*k + *ksu + j + 1])
    {
      Rprintf("j = %1i, k_j = %4.4f, k_j+1 = %4.4f\n", j, theta[*k + *ksu + j], theta[*k + *ksu + j + 1]);
      Rprintf("The supplied vector m0 is not sorted\n");
      return;
    }
  }
  
  
  // http://jakascorner.com/blog/2016/06/omp-for-reduction.html about `reduction`
#pragma omp parallel for num_threads(*Nthreds) //reduction(+: lnl0)
  for (int i = 0; i < *n; i++)
  {
    // We have K+1 possible values for y, 
    // we need to figure which this current y[i] is
    int yi;
    for (int j = 0; j < *kmu; j++)
    {
      if (y[i] == y0[j])
      {
        yi = j;
        break;
      }
    }
    // Thus, current y[i] is yi's order;
    // xb  <- x %*% theta[1:k]
    double xb = 0.0;
    for(int q = 0; q < *k; q++){
      xb += x[ q**n + i ] * theta[q];
    }
    // su2 <- exp( zsu %*% theta[(k + 1) :(k + ksu)])
    double lnsu2 = 0.0;
    for(int q = 0; q < *ksu; q++){
      lnsu2 += zsu[ q**n + i ] * theta[*k + q];
    }
    double su2s = sqrt(1.0 + exp(lnsu2));
    double su   = sqrt(exp(lnsu2));
    double ze = 1.0 / sqrt(1.0 + exp(lnsu2));
    //if(i < 20){
    //  Rprintf("Before the loop: i = %1i, yi = %1i, y[i] = %4.4f\n", 
    //          i, yi, y[i]);
    //}
    // If it the first one
    if(yi == 0)	{
      double a = theta[*k + *ksu] - xb;
      double h = a * ze;
      double b = prod[0] * su;
      double L = pnormstd(h) + 2.0 * t(h, b);
      double fh =  dnormstd(h) * pnormstd(-b*h);
      uY[i] = 2.0 * su / L * (
        b * ze * fh + pnormstd(a) / sqrt(M_PI * 2.0)
      );
      // Rprintf("i = %1i, a = %4.4f, b = %4.4f, h = %4.4f, L = %4.4f, B = %4.4f, \n", i, a, b, h, L, b* ze * fh + pnormstd(a) / sqrt(M_PI * 2.0) );
      
    } else
      // If it is the last one
      if(yi == *kmu - 1) {
        // double am = -( theta[*k + *ksu + *kmu - 2] - xb ) / su2s;
        // double b = prod[0] * su;
        // double L = pnormstd(am) + 2.0 * t(am, b);
        // uY[i] = 2.0 * su / L * (
        //   b/su2s * dnormstd(am) * pnormstd(-b*am) + 
        //     0.5*M_2_SQRTPI*M_SQRT1_2 * pnormstd(am*su2s) 
        // );
        double a = -1.0*(theta[*k + *ksu + yi - 1] - xb);
        double h = a * ze;
        double b = -1.0*(prod[0] * su);
        double L = pnormstd(h) + 2.0 * t(h, b);
        double fh =  dnormstd(-h) * pnormstd(-b*h);
        uY[i] = 2.0 * su / L * (
          b * ze * fh + pnormstd(a) / sqrt(M_PI * 2.0)
        );

      } else {
        // The rest

        double a1 = theta[*k + *ksu + yi] - xb;
        double a2 = theta[*k + *ksu + yi - 1] - xb;
        double h1 = a1 * ze;
        double h2 = a2 * ze;
        double b = prod[0] * su;
        double L = pnormstd(h1) + 2.0 * t(h1, b) -
          pnormstd(h2) - 2.0 * t(h2, b);
        double fh1 =  dnormstd(h1) * pnormstd(-b*h1);
        double fh2 =  dnormstd(h2) * pnormstd(-b*h2);
        uY[i] = 2.0 * su / L * (
          b * ze * fh1 + pnormstd(a1) / sqrt(M_PI * 2.0) -
            b * ze * fh2 - pnormstd(a2) / sqrt(M_PI * 2.0)
        );
      }
  }
}

/*
 vectorized t
 */

void TOwen(int* n, double* h, double* a, double* TOwenValues, int *Nthreds){
  //int omp_int_t  = omp_get_num_procs();
  //Rprintf("omp_int_t=%2i\n", omp_int_t);
  #pragma omp parallel for num_threads(*Nthreds)
  for(int i = 0; i < *n; i++){
    TOwenValues[i] = t(h[i], a[i]);
  }
}

void TOwen1(int* n, double* h, double* a, double* TOwenValues, int *Nthreds){
  //int omp_int_t  = omp_get_num_procs();
  //Rprintf("omp_int_t=%2i\n", omp_int_t);
  #pragma omp parallel for num_threads(*Nthreds)
  for(int i = 0; i < *n; i++){
    TOwenValues[i] = tha(h[i], 1.0, a[i], 1.0);
  }
}

double t ( double h, double a )
  
  /******************************************************************************/
  /*
   Purpose:
   
   T computes Owen's T function for arbitrary H and A.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license.
   
   Modified:
   
   13 April 2012
   
   Author:
   
   Original FORTRAN77 version by Mike Patefield, David Tandy.
   C version by John Burkardt.
   
   Reference:
   
   Mike Patefield, David Tandy,
   Fast and Accurate Calculation of Owen's T Function,
   Journal of Statistical Software,
   Volume 5, Number 5, 2000, pages 1-25.
   
   Parameters:
   
   Input, double H, A, the arguments.
   
   Output, double T, the value of Owen's T function.
   */
{
  double absa;
  double absh;
  double ah;
  const double cut = 0.67;
  double normah;
  double normh;
  double value;
  
  absh = fabs ( h );
  absa = fabs ( a );
  ah = absa * absh;
  
  if ( absa <= 1.0 )
  {
    value = tfun ( absh, absa, ah );
  }
  else if ( absh <= cut )
  {
    value = 0.25 - znorm1 ( absh ) * znorm1 ( ah ) 
    - tfun ( ah, 1.0 / absa, absh );
  }
  else
  {
    normh = znorm2 ( absh );
    normah = znorm2 ( ah );
    value = 0.5 * ( normh + normah ) - normh * normah 
      - tfun ( ah, 1.0 / absa, absh );
  }
  
  if ( a < 0.0 )
  {
    value = - value;
  }
  
  return value;
}
/******************************************************************************/


double tfun ( double h, double a, double ah )
  
  /******************************************************************************/
  /*
   Purpose:
   
   TFUN computes Owen's T function for a restricted range of parameters.
   
   Discussion:
   
   This routine computes Owen's T-function of H and A.
   
   This routine, originally named "TF", was renamed "TFUN" to avoid
   a conflict with a built in MATLAB function.
   
   Thanks to Marko Jarvenpaa for pointing out an incorrect modification
   to the code for algorithms T2 and T3 which multiplied ZNORM1 by 1/2.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license.
   
   Modified:
   
   13 July 2017
   
   Author:
   
   Original FORTRAN77 version by Mike Patefield, David Tandy.
   C version by John Burkardt.
   
   Reference:
   
   Mike Patefield, David Tandy,
   Fast and Accurate Calculation of Owen's T Function,
   Journal of Statistical Software,
   Volume 5, Number 5, 2000, pages 1-25.
   
   Parameters:
   
   Input, double H, the H argument of the function.
   0 <= H.
   
   Input, double A, the A argument of the function.
   0 <= A <= 1.
   
   Input, double AH, the value of A*H.
   
   Output, double TF, the value of Owen's T function.
   */
{
  double ai;
  double aj;
  double arange[7] = {
    0.025, 0.09, 0.15, 0.36, 0.5,
    0.9, 0.99999 };
  double as;
  double c2[21] = {
    0.99999999999999987510,
    -0.99999999999988796462,      0.99999999998290743652,
    -0.99999999896282500134,      0.99999996660459362918,
    -0.99999933986272476760,      0.99999125611136965852,
    -0.99991777624463387686,      0.99942835555870132569,
    -0.99697311720723000295,      0.98751448037275303682,
    -0.95915857980572882813,      0.89246305511006708555,
    -0.76893425990463999675,      0.58893528468484693250,
    -0.38380345160440256652,      0.20317601701045299653,
    -0.82813631607004984866E-01,  0.24167984735759576523E-01,
    -0.44676566663971825242E-02,  0.39141169402373836468E-03 };
    double dhs;
    double dj;
    double gj;
    double hrange[14] = {
      0.02, 0.06, 0.09, 0.125, 0.26,
      0.4,  0.6,  1.6,  1.7,   2.33,
      2.4,  3.36, 3.4,  4.8 };
    double hs;
    int i;
    int iaint;
    int icode;
    int ihint;
    int ii;
    int j;
    int jj;
    int m;
    int maxii;
    int meth[18] = {
      1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6 };
    double normh;
    int ord[18] = {
      2, 3, 4, 5, 7,10,12,18,10,20,30,20, 4, 7, 8,20,13, 0 };
    double pts[13] = {
      0.35082039676451715489E-02,
      0.31279042338030753740E-01,  0.85266826283219451090E-01,
      0.16245071730812277011,      0.25851196049125434828,
      0.36807553840697533536,      0.48501092905604697475,
      0.60277514152618576821,      0.71477884217753226516,
      0.81475510988760098605,      0.89711029755948965867,
      0.95723808085944261843,      0.99178832974629703586 };
    double r;
    const double rrtpi = 0.39894228040143267794;
    const double rtwopi = 0.15915494309189533577;
    int select[15*8] = {
      1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9,
      1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9,
      2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10,
      2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10,
      2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11,
      2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12,
      2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12,
      2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12 };
    double value;
    double vi;
    double wts[13] = {
      0.18831438115323502887E-01,
      0.18567086243977649478E-01,  0.18042093461223385584E-01,
      0.17263829606398753364E-01,  0.16243219975989856730E-01,
      0.14994592034116704829E-01,  0.13535474469662088392E-01,
      0.11886351605820165233E-01,  0.10070377242777431897E-01,
      0.81130545742299586629E-02,  0.60419009528470238773E-02,
      0.38862217010742057883E-02,  0.16793031084546090448E-02 };
    double y;
    double yi;
    double z;
    double zi;
    /*
     Determine appropriate method from t1...t6
     */
    ihint = 15;
    
    for ( i = 1; i <= 14; i++ )
    {
      if ( h <= hrange[i-1] )
      {
        ihint = i;
        break;
      }
    }
    
    iaint = 8;
    
    for ( i = 1; i <= 7; i++ )
    {
      if ( a <= arange[i-1] )
      {
        iaint = i;
        break;
      }
    }
    
    icode = select[ihint-1+(iaint-1)*15];
    m = ord[icode-1];
    /*
     t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18
     jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)**j / j
     aj = a**(2j-1) / (2*pi)
     */
    if ( meth[icode-1] == 1 )
    {
      hs = - 0.5 * h * h;
      dhs = exp ( hs );
      as = a * a;
      j = 1;
      jj = 1;
      aj = rtwopi * a;
      value = rtwopi * atan ( a );
      dj = dhs - 1.0;
      gj = hs * dhs;
      
      for ( ; ; )
      {
        value = value + dj * aj / ( double ) ( jj );
        
        if ( m <= j )
        {
          return value;
        }
        j = j + 1;
        jj = jj + 2;
        aj = aj * as;
        dj = gj - dj;
        gj = gj * hs / ( double ) ( j );
      }
    }
    /*
     t2(h, a, m) ; m = 10, 20 or 30
     z = (-1)^(i-1) * zi ; ii = 2i - 1
     vi = (-1)^(i-1) * a^(2i-1) * exp[-(a*h)^2/2] / sqrt(2*pi)
     */
    else if ( meth[icode-1] == 2 )
    {
      maxii = m + m + 1;
      ii = 1;
      value = 0.0;
      hs = h * h;
      as = - a * a;
      vi = rrtpi * a * exp ( - 0.5 * ah * ah );
      z = znorm1 ( ah ) / h;
      y = 1.0 / hs;
      
      for ( ; ; )
      {
        value = value + z;
        
        if ( maxii <= ii )
        {
          value = value * rrtpi * exp ( - 0.5 * hs );
          return value;
        }
        z = y * ( vi - ( double ) ( ii ) * z );
        vi = as * vi;
        ii = ii + 2;
      }
    }
    /*
     t3(h, a, m) ; m = 20
     ii = 2i - 1
     vi = a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi)
     */
    else if ( meth[icode-1] == 3 )
    {
      i = 1;
      ii = 1;
      value = 0.0;
      hs = h * h;
      as = a * a;
      vi = rrtpi * a * exp ( - 0.5 * ah * ah );
      zi = znorm1 ( ah ) / h;
      y = 1.0 / hs;
      
      for ( ; ; )
      {
        value = value + zi * c2[i-1];
        
        if ( m < i )
        {
          value = value * rrtpi * exp ( - 0.5 * hs );
          return value;
        }
        zi = y  * ( ( double ) ( ii ) * zi - vi );
        vi = as * vi;
        i = i + 1;
        ii = ii + 2;
      }
    }
    /*
     t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1
     ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)**i / (2*pi)
     */
    else if ( meth[icode-1] == 4 )
    {
      maxii = m + m + 1;
      ii = 1;
      hs = h * h;
      as = - a * a;
      value = 0.0;
      ai = rtwopi * a * exp ( - 0.5 * hs * ( 1.0 - as ) );
      yi = 1.0;
      
      for ( ; ; )
      {
        value = value + ai * yi;
        
        if ( maxii <= ii )
        {
          return value;
        }
        ii = ii + 2;
        yi = ( 1.0 - hs * yi ) / ( double ) ( ii );
        ai = ai * as;
      }
    }
    /*
     t5(h, a, m) ; m = 13
     2m - point gaussian quadrature
     */
    else if ( meth[icode-1] == 5 )
    {
      value = 0.0;
      as = a * a;
      hs = - 0.5 * h * h;
      for ( i = 1; i <= m; i++ )
      {
        r = 1.0 + as * pts[i-1];
        value = value + wts[i-1] * exp ( hs * r ) / r;
      }
      value = a * value;
    }
    /*
     t6(h, a);  approximation for a near 1, (a<=1)
     */
    else if ( meth[icode-1] == 6 )
    {
      normh = znorm2 ( h );
      value = 0.5 * normh * ( 1.0 - normh );
      y = 1.0 - a;
      r = atan ( y / ( 1.0 + a ) );
      
      if ( r != 0.0 )
      {
        value = value - rtwopi * r * exp ( - 0.5 * y * h * h / r );
      }
    }
    return value;
}
/******************************************************************************/


double znorm1 ( double z )
  
  /******************************************************************************/
  /*
   Purpose:
   
   ZNORM1 evaluates the normal CDF from 0 to Z.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license.
   
   Modified:
   
   13 April 2012
   
   Author:
   
   John Burkardt
   
   Parameters:
   
   Input, double Z, the upper limit.
   
   Output, double ZNORM1, the probability that a standard
   normal variable will lie between 0 and Z.
   */
{
  double value;
  
  value = 0.5 * erf ( z / sqrt ( 2.0 ) );
  //value = pnorm(z, 0.0, 1.0, 1.0, 0.0) - 0.5;
  
  return value;
}
/******************************************************************************/


double znorm2 ( double z )
  
  /******************************************************************************/
  /*
   Purpose:
   
   ZNORM2 evaluates the normal CDF from Z to +oo.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license.
   
   Modified:
   
   13 April 2012
   
   Author:
   
   John Burkardt
   
   Parameters:
   
   Input, double Z, the lower limit.
   
   Output, double ZNORM2, the probability that a standard
   normal variable will lie between Z and +oo.
   */
{
  double value;
  
  value = 0.5 * erfc ( z / sqrt ( 2.0 ) );
  //value = pnorm(z, 0.0, 1.0, 0.0, 0.0);
  
  return value;
}
/******************************************************************************/

// Some other code

/******************************************************************************/

double alnorm ( double x, int upper )
  
  /******************************************************************************/
  /*
   Purpose:
   
   ALNORM computes the cumulative density of the standard normal distribution.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license. 
   
   Modified:
   
   17 January 2008
   
   Author:
   
   Original FORTRAN77 version by David Hill.
   C++ version by John Burkardt.
   
   Reference:
   
   David Hill,
   Algorithm AS 66:
   The Normal Integral,
   Applied Statistics,
   Volume 22, Number 3, 1973, pages 424-427.
   
   Parameters:
   
   Input, double X, is one endpoint of the semi-infinite interval
   over which the integration takes place.
   
   Input, int UPPER, determines whether the upper or lower
   interval is to be integrated:
   1  => integrate from X to + Infinity;
   0 => integrate from - Infinity to X.
   
   Output, double ALNORM, the integral of the standard normal
   distribution over the desired interval.
   */
{
  double a1 = 5.75885480458;
  double a2 = 2.62433121679;
  double a3 = 5.92885724438;
  double b1 = -29.8213557807;
  double b2 = 48.6959930692;
  double c1 = -0.000000038052;
  double c2 = 0.000398064794;
  double c3 = -0.151679116635;
  double c4 = 4.8385912808;
  double c5 = 0.742380924027;
  double c6 = 3.99019417011;
  double con = 1.28;
  double d1 = 1.00000615302;
  double d2 = 1.98615381364;
  double d3 = 5.29330324926;
  double d4 = -15.1508972451;
  double d5 = 30.789933034;
  double ltone = 7.0;
  double p = 0.398942280444;
  double q = 0.39990348504;
  double r = 0.398942280385;
  int up;
  double utzero = 18.66;
  double value;
  double y;
  double z;
  
  up = upper;
  z = x;
  
  if ( z < 0.0 )
  {
    up = !up;
    z = - z;
  }
  
  if ( ltone < z && ( ( !up ) || utzero < z ) )
  {
    if ( up )
    {
      value = 0.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }
  
  y = 0.5 * z * z;
  
  if ( z <= con )
  {
    value = 0.5 - z * ( p - q * y 
                          / ( y + a1 + b1 
                                / ( y + a2 + b2 
                                / ( y + a3 ))));
  }
  else
  {
    value = r * exp ( - y ) 
    / ( z + c1 + d1 
          / ( z + c2 + d2 
          / ( z + c3 + d3 
          / ( z + c4 + d4 
          / ( z + c5 + d5 
          / ( z + c6 ))))));
  }
  
  if ( !up )
  {
    value = 1.0 - value;
  }
  
  return value;
}
/******************************************************************************/

void owen_values ( int *n_data, double *h, double *a, double *t )
  
  /******************************************************************************/
  /*
   Purpose:
   
   OWEN_VALUES returns some values of Owen's T function.
   
   Discussion:
   
   Owen's T function is useful for computation of the bivariate normal
   distribution and the distribution of a skewed normal distribution.
   
   Although it was originally formulated in terms of the bivariate
   normal function, the function can be defined more directly as
   
   T(H,A) = 1 / ( 2 * pi ) *
   Integral ( 0 <= X <= A ) e^(H^2*(1+X^2)/2) / (1+X^2) dX
   
   In Mathematica, the function can be evaluated by:
   
   fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
   
   Licensing:
   
   This code is distributed under the GNU LGPL license.
   
   Modified:
   
   15 December 2011
   
   Author:
   
   John Burkardt
   
   Reference:
   
   Stephen Wolfram,
   The Mathematica Book,
   Fourth Edition,
   Cambridge University Press, 1999,
   ISBN: 0-521-64314-7,
   LC: QA76.95.W65.
   
   Parameters:
   
   Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
   first call.  On each call, the routine increments N_DATA by 1, and
   returns the corresponding data; when there is no more data, the
   output value of N_DATA will be 0 again.
   
   Output, double *H, a parameter.
   
   Output, double *A, the upper limit of the integral.
   
   Output, double *T, the value of the function.
   */
{
# define N_MAX 28
  
  static double a_vec[N_MAX] = {
    0.2500000000000000E+00,
    0.4375000000000000E+00,
    0.9687500000000000E+00,
    0.0625000000000000E+00,
    0.5000000000000000E+00,
    0.9999975000000000E+00,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.1000000000000000E+02,
    0.1000000000000000E+03 };
  
  static double h_vec[N_MAX] = {
    0.0625000000000000E+00,
    6.5000000000000000E+00,
    7.0000000000000000E+00,
    4.7812500000000000E+00,
    2.0000000000000000E+00,
    1.0000000000000000E+00,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.2500000000000000E+00,
    0.2500000000000000E+00,
    0.2500000000000000E+00,
    0.2500000000000000E+00,
    0.1250000000000000E+00,
    0.1250000000000000E+00,
    0.1250000000000000E+00,
    0.1250000000000000E+00,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02 };
  
  static double t_vec[N_MAX] = {
    3.8911930234701366E-02,
    2.0005773048508315E-11,
    6.3990627193898685E-13,
    1.0632974804687463E-07,
    8.6250779855215071E-03,
    6.6741808978228592E-02,
    0.4306469112078537E-01,
    0.6674188216570097E-01,
    0.7846818699308410E-01,
    0.7929950474887259E-01,
    0.6448860284750376E-01,
    0.1066710629614485E+00,
    0.1415806036539784E+00,
    0.1510840430760184E+00,
    0.7134663382271778E-01,
    0.1201285306350883E+00,
    0.1666128410939293E+00,
    0.1847501847929859E+00,
    0.7317273327500385E-01,
    0.1237630544953746E+00,
    0.1737438887583106E+00,
    0.1951190307092811E+00,
    0.7378938035365546E-01,
    0.1249951430754052E+00,
    0.1761984774738108E+00,
    0.1987772386442824E+00,
    0.2340886964802671E+00,
    0.2479460829231492E+00 };
  
  if ( *n_data < 0 )
  {
    *n_data = 0;
  }
  
  *n_data = *n_data + 1;
  
  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *h = 0.0;
    *a = 0.0;
    *t = 0.0;
  }
  else
  {
    *h = h_vec[*n_data-1];
    *a = a_vec[*n_data-1];
    *t = t_vec[*n_data-1];
  }
  
  return;
# undef N_MAX
}
/******************************************************************************/

double tfn ( double x, double fx )
  
  /******************************************************************************/
  /*
   Purpose:
   
   TFN calculates the T-function of Owen.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license. 
   
   Modified:
   
   02 November 2010
   
   Author:
   
   Original FORTRAN77 version by JC Young, Christoph Minder.
   C version by John Burkardt.
   
   Reference:
   
   MA Porter, DJ Winstanley,
   Remark AS R30:
   A Remark on Algorithm AS76:
   An Integral Useful in Calculating Noncentral T and Bivariate
   Normal Probabilities,
   Applied Statistics,
   Volume 28, Number 1, 1979, page 113.
   
   JC Young, Christoph Minder,
   Algorithm AS 76: 
   An Algorithm Useful in Calculating Non-Central T and 
   Bivariate Normal Distributions,
   Applied Statistics,
   Volume 23, Number 3, 1974, pages 455-457.
   
   Parameters:
   
   Input, double X, FX, the parameters of the function.
   
   Output, double TFN, the value of the T-function.
   */
{
# define NG 5
  
  double fxs;
  int i;
  double r[NG] = {
    0.1477621, 
    0.1346334, 
    0.1095432, 
    0.0747257, 
    0.0333357 };
  double r1;
  double r2;
  double rt;
  double tp = 0.159155;
  double tv1 = 1.0E-35;
  double tv2 = 15.0;
  double tv3 = 15.0;
  double tv4 = 1.0E-05;
  double u[NG] = {
    0.0744372, 
    0.2166977, 
    0.3397048,
    0.4325317, 
    0.4869533 };
  double value;
  double x1;
  double x2;
  double xs;
  /*
   Test for X near zero.
   */
  if ( fabs ( x ) < tv1 )
  {
    value = tp * atan ( fx );
    return value;
  }
  /*
   Test for large values of abs(X).
   */
  if ( tv2 < fabs ( x ) )
  {
    value = 0.0;
    return value;
  }
  /*
   Test for FX near zero.
   */
  if ( fabs ( fx ) < tv1 )
  {
    value = 0.0;
    return value;
  }
  /*
   Test whether abs ( FX ) is so large that it must be truncated.
   */
  xs = - 0.5 * x * x;
  x2 = fx;
  fxs = fx * fx;
  /*
   Computation of truncation point by Newton iteration.
   */
  if ( tv3 <= log ( 1.0 + fxs ) - xs * fxs )
  {
    x1 = 0.5 * fx;
    fxs = 0.25 * fxs;
    
    for ( ; ; )
    {
      rt = fxs + 1.0;
      
      x2 = x1 + ( xs * fxs + tv3 - log ( rt ) ) 
        / ( 2.0 * x1 * ( 1.0 / rt - xs ) );
      
      fxs = x2 * x2;
      
      if ( fabs ( x2 - x1 ) < tv4 )
      {
        break;
      }
      x1 = x2;
    }
  }
  /*
   Gaussian quadrature.
   */
  rt = 0.0;
  for ( i = 0; i < NG; i++ )
  {
    r1 = 1.0 + fxs * pow ( 0.5 + u[i], 2.0 );
    r2 = 1.0 + fxs * pow ( 0.5 - u[i], 2.0 );
    
    rt = rt + r[i] * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 );
  }
  
  value = rt * x2 * tp;
  
  return value;
# undef NG
}
/******************************************************************************/

double tha ( double h1, double h2, double a1, double a2 )
  
  /******************************************************************************/
  /*
   Purpose:
   
   THA computes Owen's T function.
   
   Discussion:
   
   This function computes T(H1/H2, A1/A2) for any real numbers H1, H2, 
   A1 and A2.
   
   Licensing:
   
   This code is distributed under the GNU LGPL license. 
   
   Modified:
   
   02 November 2010
   
   Author:
   
   Original FORTRAN77 version by JC Young, Christoph Minder.
   C version by John Burkardt.
   
   Reference:
   
   Richard Boys,
   Remark AS R80:
   A Remark on Algorithm AS76:
   An Integral Useful in Calculating Noncentral T and Bivariate
   Normal Probabilities,
   Applied Statistics,
   Volume 38, Number 3, 1989, pages 580-582.
   
   Youn-Min Chou,
   Remark AS R55:
   A Remark on Algorithm AS76:
   An Integral Useful in Calculating Noncentral T and Bivariate
   Normal Probabilities,
   Applied Statistics,
   Volume 34, Number 1, 1985, pages 100-101.
   
   PW Goedhart, MJW Jansen,
   Remark AS R89:
   A Remark on Algorithm AS76:
   An Integral Useful in Calculating Noncentral T and Bivariate
   Normal Probabilities,
   Applied Statistics,
   Volume 41, Number 2, 1992, pages 496-497.
   
   JC Young, Christoph Minder,
   Algorithm AS 76: 
   An Algorithm Useful in Calculating Noncentral T and 
   Bivariate Normal Distributions,
   Applied Statistics,
   Volume 23, Number 3, 1974, pages 455-457.
   
   Parameters:
   
   Input, double H1, H2, A1, A2, define the arguments
   of the T function.
   
   Output, double THA, the value of Owen's T function.
   */
{
  double a;
  double absa;
  double ah;
  double c1;
  double c2;
  double ex;
  double g;
  double gah;
  double gh;
  double h;
  double lam;
  double twopi = 6.2831853071795864769;
  double value;
  
  if ( h2 == 0.0 )
  {
    value = 0.0;
    return value;
  }
  
  h = h1 / h2;
  
  if ( a2 == 0.0 )
  {
    g = alnorm ( h, 0 );
    
    if ( h < 0.0 )
    {
      value = g / 2.0;
    }
    else
    {
      value = ( 1.0 - g ) / 2.0;
    }
    
    if ( a1 < 0.0 )
    {
      value = - value;
    }
    return value;
  }
  
  a = a1 / a2;
  
  if ( fabs ( h ) < 0.3 && 7.0 < fabs ( a ) )
  {
    lam = fabs ( a * h );
    ex = exp ( - lam * lam / 2.0 );
    g = alnorm ( lam, 0 );
    c1 = ( ex / lam + sqrt ( twopi ) * ( g - 0.5 ) ) / twopi;
    c2 = ( ( lam * lam + 2.0 ) * ex / lam / lam / lam
             + sqrt ( twopi ) * ( g - 0.5 ) ) / ( 6.0 * twopi );
    ah = fabs ( h );
    value = 0.25 - c1 * ah + c2 * ah * ah * ah;
    if ( a < 0.0 )
    {
      value = - fabs ( value );
    }
    else
    {
      value = fabs ( value );
    }
  }
  else
  {
    absa = fabs ( a );
    
    if ( absa <= 1.0 )
    {
      value = tfn ( h, a );
      return value;
    }
    
    ah = absa * h;
    gh = alnorm ( h, 0 );
    gah = alnorm ( ah, 0 );
    value = 0.5 * ( gh + gah ) - gh * gah - tfn ( ah, 1.0 / absa );
    
    if ( a < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
