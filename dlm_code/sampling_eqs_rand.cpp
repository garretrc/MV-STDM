#define ARMA_64BIT_WORD 

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <cmath>

//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>
#include <thread>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

// Clocks
static Rcpp::Clock clock_mcmc; //time each sampling step
static Rcpp::Clock clock_ffbs; //time ffbs sampler with more detail

// Global Variables
static double b_sigma_0,a_sigma_0, a_tau_0, b_tau_0;
static arma::vec tau_curr;
static arma::uvec ind_X_alpha, ind_Y_alpha, ind_response, ind_beta_A, ind_A;
static arma::mat ytilde, ytilde_hat, sig_curr, A_curr, BTB, BTBinv, alpha_curr, Phi;
static arma::sp_mat LinvBTBinv;
static int k, nt, n, m, k_rand;

// we use Rcpp::(n,0,1) in most palces, but keep this for the alternate A sampler
// rng seed is set randomly every time for now
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::arma_rng::set_seed_random();
  arma::mat Y = arma::randn(n, ncols);
  
  arma::mat mvsamp = arma::repmat(mu, 1, n).t() + Y*arma::chol(sigma);
  return mvsamp;
}

void svd_rand(const arma::mat X, arma::mat &U, arma::vec &d, arma::mat &V){
  arma::mat Q;
  arma::mat R;
  //int m = X.n_rows;
  //int n = X.n_cols;
  
  arma::mat Omega = arma::randn(X.n_cols, k_rand+10); //random normal matrix
  arma::mat A = X*Omega; //sample matrix Y
  
  arma::qr(Q,R,A);
  arma::mat B = Q.t() * A;
  arma::svd(U,d,V,B); //arma::svd_econ(U,d,V,B,"left"); //svd econ seems slower
  U = Q*U;
  
  //now U is m x k_rand, s is k_rand x k_rand, and V is k_rand x krand
  U = U.cols(0,k_rand-1);
  d = d.subvec(0,k_rand-1);
  V = V.cols(0,k_rand-1);
  V = V.rows(0,k_rand-1);
}

static arma::mat sample_sigma() {
  arma::uvec good;
  double rate_sigma;
  arma::mat sigma2(nt,m);
  arma::uvec j_ind0(n);
  for(int i=0; i<n; i++){
    j_ind0(i) = i;
  }
  arma::uvec j_ind;
  arma::vec tmp_vec, tmp_vec2;
  double sig_tmp;
  int n_good;
  
  for(int j=0; j<m; j++){
    j_ind = j_ind0 + j*n;
    for(int t=0; t<nt; ++t){
      tmp_vec = ytilde.row(t).t() - ytilde_hat.row(t).t();
      
      tmp_vec2 = tmp_vec.elem(j_ind);
      good = find_finite(tmp_vec2); //indices of finite, non inf or nan values
      n_good = size(tmp_vec2.elem(good))(0);
      if(n_good>0){
        rate_sigma = b_sigma_0 + 0.5*dot(tmp_vec2.elem(good),tmp_vec2.elem(good));
        sig_tmp = R::rgamma(a_sigma_0 + n_good/2,1/rate_sigma); //R::rgamma uses scale
      }else{ //if entire time period is missing, i don't think we should estimate sigmat using only prior info
        sig_tmp = R::rgamma(a_sigma_0, 1/b_sigma_0); //R::rgamma uses scale
      }
      sigma2(t,j) = 1/sig_tmp;
    }
  }
  
  return sigma2;
}


static arma::vec sample_tau() {
  arma::vec tau2(m);
  double shape_tau = a_tau_0 + k*nt/2;
  double rate_tau;
  double tau_tmp;
  
  arma::mat tmpnu(1,m*k);
  arma::mat nu(1,k);
  arma::vec tmprate(nt);
  arma::mat tmpmat(1,1);
  arma::uvec ind0(k);
  arma::uvec ind;
  
  for(int i=0; i<k; i++){
    ind0(i) = i;
  }
  
  for(int j=0; j<m; j++){
    ind = ind0 + j*k;
    for(int t=0; t<nt; t++){
      //Using cholesky (check speed)
      tmpnu = alpha_curr.row(t+1).t() - A_curr*alpha_curr.row(t).t(); //edit
      nu = LinvBTBinv*tmpnu(ind); 
      tmprate(t) = dot(nu,nu);
      
      //Using full BTB matrix
      //nu = alpha_curr.row(t+1).t() - A_curr*alpha_curr.row(t).t();
      //tmprate(t) = as_scalar(nu(ind).t()*BTB*nu(ind));
    }
    rate_tau = b_tau_0 + 0.5*accu(tmprate);
    tau_tmp = R::rgamma(shape_tau,1/rate_tau);
    tau2(j) = 1/tau_tmp;
  }
  
  return tau2;
}

//This code not working after fixing the sampler for tau2 so I implemented a version without accounting for Sigma_alpha
//Even though we tried this cholesky trick, it doesn't result in proper samples anyways, so I'll rewrite sample_A() with
//the block-wise vector AR structure
// static arma::mat sample_A() { 
//   
//   //Update Innovation Variance
//   arma::mat tau_mat = arma::diagmat(tau_curr);
//   arma::mat Sigma_alpha_A = kron(tau_mat,BTBinv);
//   arma::mat ar1_in = solve(arma::chol(Sigma_alpha_A),alpha_curr.t()); //sigma_alpha defined differently across scruots so i hope this is correct
//   arma::vec ar1_out(k);
//   arma::vec atmp(nt);
//   arma::vec f1(nt-1);
//   arma::vec f2(nt-1);
//   
//   for(int i=0; i<k; i++){
//     atmp = ar1_in.row(i).t();
//     f1 = atmp.subvec(0,nt-1);
//     f2 = atmp.subvec(1,nt);
//     double FF_inv = 1/dot(f1,f1);
//     double a_hat = dot(f1,f2)*FF_inv; //% is element-wise multiplication in c++
//     ar1_out(i) = R::rnorm(a_hat, sqrt(tau_curr(0,0)*FF_inv));
//   }
//   
//   arma::mat A = diagmat(ar1_out);
//   
//   return A;
// }


//univariate case, same as above, without using cholesky. Samples are bad though...
static arma::mat sample_A_ar() { 
  
  //Update Innovation Variance
  arma::vec ar1_out(k);
  arma::vec atmp(nt);
  arma::vec f1(nt-1);
  arma::vec f2(nt-1);
  double FF_inv;
  double a_hat;
  
  for(int i=0; i<k; i++){
    atmp = alpha_curr.col(i);
    f1 = atmp.subvec(0,nt-1);
    f2 = atmp.subvec(1,nt);
    FF_inv = 1/dot(f1,f1);
    a_hat = dot(f1,f2)*FF_inv;
    ar1_out(i) = R::rnorm(a_hat, sqrt(tau_curr(0)*FF_inv)); //tau_curr(i)?
  }
  
  arma::mat A = diagmat(ar1_out);
  
  return A;
}
//implement fixed spatial A matrix

//Vector AR gibbs sampler for A, assuming one m by m block of coefficients for each basis center
//This is using a flat prior, we could implement minnesota as well.
static arma::mat sample_A() { 
  arma::vec Y_alpha(k*m*(nt+1));
  arma::mat X_alpha_tmp(k*m*nt,m*m*k);
  arma::vec z(k*m*m);
  //arma::sp_mat X_alpha(k*m*nt,m*m*k);
  arma::sp_mat L_inv_Sigma_full(k*m*nt,k*m*nt);
  arma::sp_mat I_nt = arma::speye(nt,nt);
  arma::sp_mat sqrt_tau_inv(m,m);
  arma::mat L_inv_X(k*m*nt,m*m*k);
  arma::mat L_inv_Y(k*m*nt,1);
  arma::mat precision_beta_A(m*m*k,m*m*k);
  arma::mat beta_A_hat(m*m*k,1);
  arma::mat beta_A(m*m*k,1);
  arma::mat A_out(k*m,k*m, arma::fill::zeros);
  
  sqrt_tau_inv = inv(diagmat(sqrt(tau_curr)));
  
  //Sample A matrix (Gibbs sample with flat prior)
  //Create vector by reordering the alphas by time, variable, basis
  Y_alpha = vectorise(alpha_curr.t());
  
  //Place the correct elements of Y_alpha into X_alpha (sparse design matrix)
  X_alpha_tmp(ind_X_alpha) = Y_alpha(ind_Y_alpha);
  arma::sp_mat X_alpha(X_alpha_tmp);
  
  //Cholesky decomposition of linear model covariance inverse
  L_inv_Sigma_full = kron(I_nt,kron(sqrt_tau_inv,LinvBTBinv));
  
  //Use cholesky to give diagonal variance
  L_inv_Y = L_inv_Sigma_full*Y_alpha(ind_response);
  L_inv_X = L_inv_Sigma_full*X_alpha;
  
  //posterior precision and mean
  //Lot of big multiplcations and solves here
  //the added identity matrix comes from the prior on the variance
  int c = 9;                                                //c is equivalent to the ridge penalty
  arma::mat precision_prior = c*arma::eye(m*m*k,m*m*k);     //prior precision matrix
  arma::vec mean_prior = arma::vec(m*m*k,arma::fill::zeros);
  int ind = 0;
  for(int i=0; i<m; i++){
    for(int h=0; h<k; h++){
      for(int j=0; j<m; j++){
        if(i==j){
          mean_prior(ind) = 1;
        }
        ind+=1;
        }
     }
  }
  //arma::mat precision_mean_prior = c*arma::diagmat(mean_prior); //prior for precision times prior for the mean
  precision_beta_A = L_inv_X.t()*L_inv_X + precision_prior; //this step takes the longest (big matrix multiplication)
  beta_A_hat = solve(precision_beta_A, L_inv_X.t()*L_inv_Y + c*mean_prior);
  
  //draw posterior sample
  //z = Rcpp::rnorm(m*m*k,0,1);
  z = arma::randn(m*m*k);
  //in R, we call backsolve() here, which is fast for upper triangular
  //there's no backsolve in armadillo, so I wonder how much this will slow us down...
  beta_A = beta_A_hat + solve(chol(precision_beta_A),z);
  //finally, place the beta_A terms in the correct spot
  A_out(ind_A) = beta_A(ind_beta_A);
  return A_out;
}

// MV option (fully dense)
static arma::mat sample_A_var() {
  
  double c0 = 1;//alpha_t prior variance
  //int p = 1; //order 1
  //int i = 0; //no intercept, i=1 is intercept
  arma::mat X = alpha_curr;
  arma::mat XX = X.t() * X;
  arma::mat Vprior_inv = c0 * arma::eye(k*m,k*m); //may want this to be X'X;
  arma::mat Vpost = inv(Vprior_inv + XX);
  
  arma::mat Ahat = solve(XX,  X.t() * alpha_curr);
  arma::mat Aprior=Ahat; //prior mean is OLS estimate?
  arma::mat Apost = Vpost * (Vprior_inv * Aprior + XX * Ahat);
  arma::mat Amat = mvrnormArma(1, Apost.as_col(), kron(kron(arma::diagmat(tau_curr), BTBinv), Vpost));
  arma::mat A(k*m, k*m); //.reshape not working
  
  for(unsigned int i=0; i<A.n_cols; i++){
    A.col(i) = Amat.cols(i*A.n_cols, (i+1)*A.n_cols -1).t();
  }
  return A;
}

//--------------------------------------------------------------------------//
//standard kalman filter
//--------------------------------------------------------------------------//
static arma::mat sample_alpha() {
  clock_ffbs.tick("setup");
  
  //double eps = arma::datum::eps;
  double eps = 0.00002;
  int i, row, col, ind;
  //only implemented for our config but could read in the list of dlm params and implement all the flags to make more general
  arma::vec m0(k*m,arma::fill::zeros);
  arma::mat C0(k*m,k*m, arma::fill::zeros);
  C0.diag().fill(10);
  arma::mat FF = kron(arma::eye(m,m),Phi);
  arma::sp_mat GG = arma::sp_mat(A_curr);
  arma::mat V(n*m,n*m, arma::fill::zeros);
  arma::sp_mat Vinv(n*m,n*m);
  arma::sp_mat sqrtVinv(n*m,n*m);
  arma::mat tau_mat = arma::diagmat(tau_curr);
  arma::mat W = kron(tau_mat,BTBinv); 
  arma::umat JV, tmp_mat, tmp_mat1;
  tmp_mat.eye(m,m);
  for(i=0; i<m; i++){
    tmp_mat(i,i)=i+1;
  }
  tmp_mat1.eye(n,n);
  arma::umat JV_init =  kron(tmp_mat,tmp_mat1);// #The index "i" means that at time t, replace the matching value in V with X[t,i]
  arma::mat X = sig_curr;
  arma::mat y = ytilde; //should remove and rename
  
  //move all objects upfront eventually
  arma::vec tmp_d, Dv, Dv_inv, DR_inv, foo, tmp_vec;
  arma::umat v_index, x_index;
  arma::mat tF_Vinv, sqrtW, tmp_U, tmp_V;
  
  //missing data objects
  arma::sp_mat tF_VinvTMP;
  arma::sp_mat sqrtVinvTMP;
  arma::sp_mat VinvTMP;
  arma::uvec good;
  arma::mat yTMP, fTMP;

  //should build in flag functions for real
  bool tvW = false;
  bool tvFF = false;
  bool tvV = true;
  bool tvGG = false;
  bool tvFV = tvFF || tvV;
  
  if(tvV){
    arma::uvec JV_ind = find(JV_init!=0);
    arma::uvec colvec = floor(JV_ind/JV_init.n_rows);
    arma::uvec rowvec = JV_ind - colvec*JV_init.n_rows;
    arma::uvec jv_tmp = JV_init.elem(JV_ind)-1;
    JV = join_horiz(rowvec,colvec,jv_tmp);
  }
  
  arma::mat mfill;
  mfill.zeros(y.n_rows,m0.n_elem);
  arma::mat mmat = arma::join_cols(m0.t(),mfill);
  arma::mat a(y.n_rows, m0.n_elem, arma::fill::zeros);
  arma::mat f(y.n_rows,y.n_cols, arma::fill::zeros);
  //arma::cube UC(m0.n_elem,m0.n_elem,y.n_rows+1,arma::fill::zeros);
  //arma::mat DC(y.n_rows+1, m0.n_elem,arma::fill::zeros);
  //arma::cube UR(m0.n_elem,m0.n_elem,y.n_rows,arma::fill::zeros);
  //arma::mat DR(y.n_rows,m0.n_elem, arma::fill::zeros);
  arma::cube UC(m0.n_elem, k_rand, y.n_rows+1,arma::fill::zeros);
  arma::mat DC(y.n_rows+1, k_rand, arma::fill::zeros);
  arma::cube UR(m0.n_elem, k_rand, y.n_rows,arma::fill::zeros);
  arma::mat DR(y.n_rows, k_rand, arma::fill::zeros);
  
  if ( !tvW ) {
    tmp_U.clear();
    tmp_d.clear();
    tmp_V.clear();
    arma::svd_econ(tmp_U,tmp_d,tmp_V,W, "right");
    arma::mat sqrtW_tmp(tmp_V.n_rows, tmp_V.n_rows);
    //tmp_d.elem(find(sqrt(tmp_d)<eps)).fill(eps);
    tmp_V = tmp_V.t();
    for(i=0; i<m*k; i++){
      sqrtW_tmp.row(i) = sqrt(tmp_d(i)) * tmp_V.row(i);
    }
    sqrtW = sqrtW_tmp;
    //test if sqrtW is calculated correctly
    //sqrtW = sqrtW_tmp.t()*sqrtW_tmp;
  }
  
  arma::mat tmp_vt, tmp_matU, Uv, tmp_prior, joined_ab;
  
  //UC.slice(0) = arma::eye(C0.n_rows,C0.n_rows);
  //DC.row(0) = sqrt(C0.diag().t());
  UC.slice(0) = arma::eye(m0.n_elem,k_rand);
  arma::vec temp = sqrt(C0.diag());
  DC.row(0) = temp.subvec(0,k_rand-1).t();
  
  //--------------------------------------------------------------------------//
  //Filtering
  //should check this or replace this with C_dlmFilter0
  //--------------------------------------------------------------------------//
  //Rcpp::Rcout << "filtering" << std::endl;
  clock_ffbs.tock("setup");
  clock_ffbs.tick("ff");
  for(int t=0; t<nt; t++){
    clock_ffbs.tick("ff: Vinv, sqrtVinv");
    if (tvV) {
      //setup V
      for(unsigned int i=0; i<JV.n_rows; i++){ //this could be its own function
        row = JV(i,0);
        col = JV(i,1);
        ind = JV(i,2);
        V(row, col) = X(t,ind); //is this correct?
        //Vinv(row, col) = 1/X(t,ind); //is this correct?
      }
      
      //Vinv, sqrtVinv in the diagonal case
      Vinv.diag() = 1/V.diag();
      sqrtVinv.diag() = sqrt(1/V.diag());
      
      // //above should be way faster, but gives a different square root than below code
      // //matrices can have many square roots, I don't think it matters which we use
      // tmp_d = V.diag();
      // tmp_vt.eye(V.n_rows,V.n_rows);
      // Dv = sqrt(tmp_d);
      // Dv_inv = 1/Dv;
      // Dv_inv.replace(arma::datum::inf, 0);
      // sqrtVinv = tmp_vt;
      // for(unsigned int i=0; i<sqrtVinv.n_rows; i++){
      //   sqrtVinv.row(i) *= Dv_inv(i);
      // }

  }
  clock_ffbs.tock("ff: Vinv, sqrtVinv");

    // if ( tvW ) { //should build this in if we want time varying Sigma_alpha?
    //   mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
    //   W <- La.svd(mod$W,nu=0)
    //   sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W
    // }

  clock_ffbs.tick("ff: tF_Vinv");
  if ( tvFV ) { 
    //tF_Vinv = FF.t()/V(0,0);
    tF_Vinv = FF.t()*Vinv;
    //tF_Vinv = FF.t()*(sqrtVinv.t()*sqrtVinv);
  }
  clock_ffbs.tock("ff: tF_Vinv");

  clock_ffbs.tick("ff: find_finite");
  good.clear();
  good = find_finite(y.row(t)); //indices of finite, non inf or nan values
  clock_ffbs.tock("ff: find_finite");

  if ( good.n_elem==y.n_cols ) {//no missing values
    //Rcpp::Rcout <<"t= "<< t << ", complete data" << std::endl;
  
    //prior
    a.row(t) = (GG*mmat.row(t).t()).t();
    tmp_prior = GG*UC.slice(t);
    
    for(unsigned int i=0; i<tmp_prior.n_cols; i++){
      tmp_prior.col(i) *= DC(t,i);
    }
    joined_ab  = arma::join_horiz(tmp_prior, sqrtW);
    tmp_U.clear();
    tmp_d.clear();
    tmp_V.clear();
    
    //clock_ffbs.tick("ff: svd_econ prior");
    //arma::svd_econ(tmp_U,tmp_d,tmp_V,joined_ab, "right");
    //UR.slice(t) = tmp_V; 
    //DR.row(t) = tmp_d.t();
    //clock_ffbs.tock("ff: svd_econ prior");
    
    clock_ffbs.tick("ff: svd prior");
    //Rcpp::Rcout << t << std::endl;
    //TODO: only use randomized for t=0
    svd_rand(joined_ab,tmp_U,tmp_d,tmp_V);
    UR.slice(t) = tmp_U; 
    DR.row(t) = tmp_d.t();
    clock_ffbs.tock("ff: svd prior");
    
    
    //one-step forecast
    clock_ffbs.tick("ff: one step ahead");
    f.row(t) =  a.row(t) * FF.t();
    clock_ffbs.tock("ff: one step ahead");
    
    //posterior
    clock_ffbs.tick("ff: DR_inv");
    DR_inv = 1/DR.row(t).t();
    clock_ffbs.tock("ff: DR_inv");
    
    clock_ffbs.tick("ff: posterior join");
    DR_inv.replace(arma::datum::inf, 0);
    joined_ab = arma::join_vert(sqrtVinv * FF * UR.slice(t), arma::diagmat(DR_inv));
    clock_ffbs.tock("ff: posterior join");
    
    clock_ffbs.tick("ff: clear posterior");
    tmp_U.clear();
    tmp_d.clear();
    tmp_V.clear();
    clock_ffbs.tock("ff: clear posterior");
    clock_ffbs.tick("ff: svd_econ posterior");
    //arma::svd(tmp_U, tmp_d, tmp_V, joined_ab); //this svd is the most expensive thing
    arma::svd_econ(tmp_U, tmp_d, tmp_V, joined_ab, "right", "std");
    clock_ffbs.tock("ff: svd_econ posterior");
    
    clock_ffbs.tick("ff: U*tmp_V");
    UC.slice(t+1) = UR.slice(t) * tmp_V;
    clock_ffbs.tock("ff: U*tmp_V");
    
    clock_ffbs.tick("ff: DC");
    foo = 1/tmp_d;
    foo.replace(arma::datum::inf, 0);
    DC.row(t+1) = foo.t();
    clock_ffbs.tock("ff: DC");
    
    clock_ffbs.tick("ff: tmp_matU");
    tmp_matU = UC.slice(t+1).t(); 
    clock_ffbs.tock("ff: tmp_matU");
    
    clock_ffbs.tick("ff: tmp_vec setup");
    tmp_vec.clear();
    tmp_vec = DC.row(t+1).t(); 
    for(unsigned int i=0; i<tmp_matU.n_rows; i++){
      tmp_matU.row(i) *= tmp_vec(i);
    }
    tmp_vec.clear();
    clock_ffbs.tock("ff: tmp_vec setup");
    
    clock_ffbs.tick("ff: calculate mmat");
    tmp_vec = (tmp_matU.t() * tmp_matU) * tF_Vinv * (y.row(t)-f.row(t)).t(); //high values driven by high y values
    mmat.row(t+1) = a.row(t) + tmp_vec.t(); 
    clock_ffbs.tock("ff: calculate mmat");
  } else if(good.n_elem==0){ //all missing values
    //Rcpp::Rcout <<"t= "<< t << ", all missing" << std::endl;
    //prior
    clock_ffbs.tick("ff: svd prior + join all missing");
    a.row(t) = (GG*mmat.row(t).t()).t();
    tmp_prior = GG*UC.slice(t);
    for(unsigned int i=0; i<tmp_prior.n_cols; i++){
      tmp_prior.col(i) *= DC(t,i);
    }
    //joined_ab  = arma::join_cols(tmp_prior.t(), sqrtW); 
    joined_ab  = arma::join_horiz(tmp_prior, sqrtW); 
    tmp_U.clear();
    tmp_d.clear();
    tmp_V.clear();
    //arma::svd_econ(tmp_U,tmp_d,tmp_V,joined_ab, "right");
    svd_rand(joined_ab,tmp_U,tmp_d,tmp_V);
    UR.slice(t) = tmp_U; 
    DR.row(t) = tmp_d.t();
    UC.slice(t+1) = tmp_U; 
    DC.row(t+1) = tmp_d.t();
    
    //posterior
    mmat.row(t+1) = a.row(t);
    
    //one-step forecast
    f.row(t) =  a.row(t) * FF.t();
    clock_ffbs.tock("ff: svd prior + join all missing");
  } else{ //some components missing
    //Rcpp::Rcout <<"t= "<< t << ", some missing" << std::endl;
    
    clock_ffbs.tick("ff: sqrtVinv some missing");
    tmp_d = V.diag();
    
    sqrtVinvTMP.eye(good.n_elem,good.n_elem);
    sqrtVinvTMP.diag() = sqrt(1/tmp_d(good));
    clock_ffbs.tock("ff: sqrtVinv some missing");
    
    clock_ffbs.tick("ff: tF*Vinv some missing");
    VinvTMP.eye(good.n_elem,good.n_elem);
    VinvTMP.diag() = 1/tmp_d(good);
    tF_VinvTMP = FF.rows(good).t()*VinvTMP;
    clock_ffbs.tock("ff: tF*Vinv some missing");
    
    //prior
    clock_ffbs.tick("ff: svd prior some missing");
    a.row(t) = (GG*mmat.row(t).t()).t();
    tmp_prior = GG*UC.slice(t);
    for(unsigned int i=0; i<tmp_prior.n_cols; i++){
      tmp_prior.col(i) *= DC(t,i);
    }
    //joined_ab  = arma::join_cols(tmp_prior.t(), sqrtW); 
    joined_ab  = arma::join_horiz(tmp_prior, sqrtW); 
    tmp_U.clear();
    tmp_d.clear();
    tmp_V.clear();
    //arma::svd_econ(tmp_U,tmp_d,tmp_V,joined_ab, "right");
    svd_rand(joined_ab,tmp_U,tmp_d,tmp_V);
    UR.slice(t) = tmp_U; 
    DR.row(t) = tmp_d.t();
    clock_ffbs.tock("ff: svd prior some missing");
    
    //one-step forecast
    clock_ffbs.tick("ff: one step some missing");
    f.row(t) =  a.row(t) * FF.t();
    clock_ffbs.tock("ff: one step some missing");
    
    //posterior
    clock_ffbs.tick("ff: svd posterior some missing");
    DR_inv = 1/DR.row(t).t();
    DR_inv.replace(arma::datum::inf, 0);
    //joined_ab = arma::join_cols(sqrtVinvTMP * FF.rows(good) * UR.slice(t), arma::diagmat(DR_inv));
    joined_ab = arma::join_vert(sqrtVinvTMP * FF.rows(good) * UR.slice(t), arma::diagmat(DR_inv));
    tmp_U.clear();
    tmp_d.clear();
    tmp_V.clear();
    //arma::svd(tmp_U, tmp_d, tmp_V, joined_ab); //this svd is the most expensive thing, dc not working right
    //arma::svd_econ(tmp_U, tmp_d, tmp_V, joined_ab, "right");
    arma::svd_econ(tmp_U, tmp_d, tmp_V, joined_ab, "right", "std");
    clock_ffbs.tock("ff: svd posterior some missing");
    
    
    clock_ffbs.tick("ff: tmp_matU some missing");
    UC.slice(t+1) = UR.slice(t) * tmp_V;
    foo = 1/tmp_d;
    foo.replace(arma::datum::inf, 0);
    DC.row(t+1) = foo.t();
    tmp_matU = UC.slice(t+1).t(); //check if this should be transpose
    tmp_vec.clear();
    tmp_vec = DC.row(t+1).t(); 
    for(unsigned int i=0; i<tmp_matU.n_rows; i++){
      tmp_matU.row(i) *= tmp_vec(i);
    }
    clock_ffbs.tock("ff: tmp_matU some missing");
    
    clock_ffbs.tick("ff: mmat some missing");
    tmp_vec.clear();
    yTMP.clear();
    fTMP.clear();
    yTMP = y.cols(good);
    fTMP = f.cols(good);
    tmp_vec = (tmp_matU.t() * tmp_matU) * tF_VinvTMP * (yTMP.row(t)-fTMP.row(t)).t(); //high values driven by high y values
    mmat.row(t+1) = a.row(t) + tmp_vec.t(); 
    clock_ffbs.tock("ff: mmat some missing");
  }
}
clock_ffbs.tock("ff");
//--------------------------------------------------------------------------//
  //Backwards Sampling
//--------------------------------------------------------------------------//
  
  clock_ffbs.tick("bs");
//Rcpp::Rcout << "back sampling" << std::endl;
arma::mat tG_Winv, sqrtWinv, Dinv, UH, tmp_math;
arma::vec Dw, Dw_inv, DH, h;
int p = mmat.n_cols;

arma::mat theta(nt + 1, p, arma::fill::zeros);
if (!tvW) {
  tmp_U.clear();
  tmp_d.clear();
  tmp_V.clear();
  arma::svd_econ(tmp_U, tmp_d, tmp_V, W, "right");
  tmp_V = tmp_V.t();
  Dw = sqrt(tmp_d);
  Dw.elem( find(Dw<eps) ).fill(eps);
  Dw_inv = 1/Dw;
  sqrtWinv.zeros(W.n_rows, W.n_cols);
  for(unsigned int i=0; i<sqrtWinv.n_rows; i++){
    sqrtWinv.row(i) = Dw_inv(i) * tmp_V.row(i);
  }
  if (!tvGG){
    tG_Winv = GG.t() * sqrtWinv.t() * sqrtWinv;
  }
}

tmp_vec.clear();
//tmp_vec = Rcpp::rnorm(p,0,1);
// for(i=0; i<p; i++){
  //   tmp_vec(i) *= DC(nt,i);
  // }
//tmp_vec = Rcpp::rnorm(k_rand,0,1);
tmp_vec = arma::randn(k_rand);
for(i=0; i<k_rand; i++){
  tmp_vec(i) *= DC(nt,i);
}
theta.row(nt) = mmat.row(nt) + (UC.slice(nt) * tmp_vec).t();

//Rcpp::Rcout << "back sampling loop start" << std::endl;
for(int t=nt-1; t>=0; t--){
  Dinv = 1/DC.row(t);
  Dinv.replace(arma::datum::inf, 0);
  
  tmp_U.clear();
  tmp_d.clear();
  tmp_V.clear();
  arma::mat tmp_matW = sqrtWinv * GG * UC.slice(t);
  arma::mat tmp_matW2 = join_cols(tmp_matW, arma::diagmat(Dinv));
  arma::svd_econ(tmp_U, tmp_d, tmp_V, tmp_matW2, "right");
  //UH = UC.slice(t) * tmp_V.t(); //is this the same as tmp_U?
    UH = UC.slice(t) * tmp_V; //is this the same as tmp_U?
      DH = 1/tmp_d;
      DH.replace(arma::datum::inf,0);
      
      
      tmp_math=UH.t();
      for(unsigned int i=0; i<UH.n_cols; i++){
        tmp_math.row(i) *=DH(i);
      }
      
      tmp_vec.clear();
      tmp_vec = (tmp_math.t() * tmp_math) * tG_Winv * (theta.row(t+1) - a.row(t)).t();
      
      h = mmat.row(t).t() + tmp_vec; //mmat.row(0) = 0
      tmp_vec.clear();
      // tmp_vec = Rcpp::rnorm(p,0,1);
      // for(i=0; i<p; i++){
        //   tmp_vec(i) *= DH(i);
        // }
      //tmp_vec = Rcpp::rnorm(k_rand,0,1);
      tmp_vec = arma::randn(k_rand);
      for(i=0; i<k_rand; i++){
        tmp_vec(i) *= DH(i);
      }
      
      //Rcpp::Rcout << "theta(nt), " <<min(theta.row(nt)) <<", "<< mean(theta.row(nt))<<", " << max(theta.row(nt)) << std::endl;
      theta.row(t) = (h + UH * tmp_vec).t();
}
clock_ffbs.tock("bs");

return(theta);
}

//sig_init should be nt x m matrix
//added a bunch of static index vectors in the arguments to avoid rewriting the code in c++
//Eventually, we can generate the ind_ vectors in c++ or in an R wrapper
// [[Rcpp::export]]
List mcmc_dlm(arma::mat y, arma::mat B, arma::mat Phi_basis, int n_mcmc, int n_burn, int k_svd,
              arma::mat alpha_init, arma::mat A_init,
              arma::mat sig_init, arma::vec tau_init,
              double a_sigma, double b_sigma, double a_tau, double b_tau) {
  //make blank clocks to overwrite the global clocks
  //This way, the clocks will refresh each time you call mcmc_dlm
  Rcpp::Clock c1;
  Rcpp::Clock c2;
  clock_mcmc = c1;
  clock_ffbs = c2;
  
  Phi = Phi_basis;
  ytilde = y;
  nt = ytilde.n_rows;
  k = Phi.n_cols;
  n = Phi.n_rows;
  m = y.n_cols/n;
  k_rand = k_svd;
  
  //prior params
  a_sigma_0 = a_sigma;
  b_sigma_0 = b_sigma;
  a_tau_0 = a_tau;
  b_tau_0 = b_tau;
  
  // initialize parameters
  sig_curr = sig_init;
  tau_curr = tau_init;
  A_curr = A_init;
  alpha_curr = alpha_init;
  
  
  ytilde_hat = alpha_curr.submat(1,0,nt,k*m-1)*kron(arma::eye(m,m),Phi).t(); //nt x n
  //upfront computations
  //Binv = solve(B);
  BTB = B.t() * B;
  BTBinv = inv(BTB);
  LinvBTBinv = inv(chol(BTBinv).t());
  
  arma::mat tau_samps(m, n_mcmc);
  arma::cube sig_samps(nt,m, n_mcmc);
  arma::cube A_samps(k*m,k*m,n_mcmc);
  arma::cube alpha_samps(nt+1,k*m,n_mcmc); //keep alpha0
  
  //sample_A setup, indices for estimating the nonzero terms of the A matrix as a vector
  arma::uvec beta_A(k*m*m); 
  for(unsigned int i=0; i < beta_A.n_elem ; i++){
    beta_A(i) = i;
  }
  arma::umat Aind(m*k,m*k);
  arma::umat Atmp_mat1(m,m);
  Atmp_mat1.fill(1);
  arma::umat Atmp_mat2;
  Atmp_mat2.eye(k,k);
  Aind = kron(Atmp_mat1,Atmp_mat2);
  ind_A = find(Aind); //find nonzero elements
  Aind(ind_A) = beta_A;
  Aind = Aind.t();
  ind_beta_A = Aind(ind_A);
  
  arma::uvec ind_response_init(nt*k*m); 
  for(int i=0; i < (nt*k*m); i++){
    ind_response_init(i) = i+k*m;
  }
  ind_response = ind_response_init;
  
  arma::umat X_alpha(nt*k*m, m*m*k, arma::fill::zeros);
  arma::umat Xalpha_tmp(k,k*m, arma::fill::zeros);
  arma::umat meye(m,m);
  meye.eye();
  for(int t=0; t<nt; t++){
    for(int j=0;j<m; j++){
      for(int i=0; i<k; i++){
        Xalpha_tmp(i,i*m+j) = j*k+i + t*k*m + 1;
      }
    }
    X_alpha.rows(t*k*m,(t+1)*k*m-1) = kron(meye, Xalpha_tmp);
  }
  
  ind_X_alpha = find(X_alpha);
  ind_Y_alpha = X_alpha(ind_X_alpha)-1;
  
  // begin mcmc
  
  Progress p(n_mcmc, print);
  for (int iter = 0; iter < n_mcmc; iter++) {
    //if (Progress::check_abort())
      //  return -1.0;
    //Rcpp::Rcout << "iter = " << iter << std::endl;
    
    clock_mcmc.tick("sigma");
    sig_curr = sample_sigma();
    clock_mcmc.tock("sigma");
    
    clock_mcmc.tick("tau");
    tau_curr = sample_tau();
    clock_mcmc.tock("tau");
    
    clock_mcmc.tick("A");
    A_curr = sample_A();
    clock_mcmc.tock("A");
    
    clock_mcmc.tick("alpha");
    alpha_curr = sample_alpha();
    clock_mcmc.tock("alpha");
    
    tau_samps.col(iter) = tau_curr;
    sig_samps.slice(iter) = sig_curr;
    A_samps.slice(iter) = A_curr;
    alpha_samps.slice(iter) = alpha_curr;
    
    ytilde_hat = alpha_curr.submat(1,0,nt,k*m-1)*kron(arma::eye(m,m),Phi).t(); //nt x n
    
    p.increment();  // update progress
  }
  clock_mcmc.stop("clock_mcmc");
  clock_ffbs.stop("clock_ffbs");
  
  //subetting cubes with .slices to implement burnin creating vectors instead of cubes...
  List samps_list = List::create(Rcpp::Named("sigma") = sig_samps, Rcpp::Named("tau") = tau_samps,
                                 Rcpp::Named("A") = A_samps, Rcpp::Named("alpha") = alpha_samps,
                                 Rcpp::Named("Phi") = Phi);
  
  return (samps_list);
}


