data {
  int<lower=1> N;                   // number of reps
  int<lower=1> S;                   // number of species 
  int<lower=1, upper=S> focal[N];   // focal species index
  vector<lower=0>[S] x[N];          // density of competitors 
  vector<lower=0>[N] y;             // response (seeds produced)
  vector<lower=0>[N] h;             // product of competitor density for HOI 
}
transformed data { 
  vector<lower=0>[N] u; 
  
  u = 1 ./ y; 
   
}
parameters {
  vector<lower=0>[S] lambda;        // gain for each species
  vector<lower=0>[S*S] alpha_vec;   // comp 
  vector<lower=0>[S*S] tau_vec;     // curvature on species effects
  vector[S] beta;          // HOI effects  
  vector<lower=0>[S] sigma;         // std dev for each species 
}
transformed parameters{
  vector[N] mu;                     // linear predictor
  matrix<lower=0>[S,S] alpha;
  matrix<lower=0>[S,S] tau; 
  
  alpha = to_matrix(alpha_vec, S, S ); 
  tau = to_matrix(tau_vec, S,S); 
  
  for ( i in 1:N) {
    vector[S] C;        // sum of single species effects 
    real HOI; 
    
    for( j in 1:S){   
      C[j] = (alpha[focal[i], j]*x[i,j])^tau[focal[i], j]; 
    }
    HOI = beta[focal[i]]*(h[i]); // HOI term 
    
    mu[i] = (1 + sum(C) + HOI)/lambda[focal[i]];  // matrix mult alpha by x vector 
  }
}
model {
  // priors
  lambda ~ cauchy(0, 5); 
  sigma ~ cauchy(0,2);
  alpha_vec ~ cauchy(0,1); 
  tau_vec ~ cauchy(0,1);
  beta ~ cauchy(0,1); 

  // likelihood
  for( i in 1:N)
    u[i] ~ normal( mu[i] , sigma[focal[i]] );
}
generated quantities{ 
  vector[N] y_hat; 
  vector[N] u_hat; 
  
  for(i in 1:N)
    u_hat[i] = normal_rng( mu[i], sigma[focal[i]] );
    
  y_hat = 1 ./ u_hat; 
}
