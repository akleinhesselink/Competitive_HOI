data {
  int<lower=1> N;                   // number of reps
  int<lower=1> S;                   // number of species 
  int<lower=1, upper=S> focal[N];   // focal species index
  vector<lower=0>[S] x[N];          // density of competitors 
  vector<lower=0>[N] y;             // response (seeds produced)
}
parameters {
  vector<lower=0>[S] lambda;        // gain for each species
  vector<lower=0>[S*S] alpha_vec;   // comp 
  vector<lower=0>[S*S] tau_vec;     // curvature on species effects 
  real<lower=0> sigma;              // std dev 
}
transformed parameters{
  vector[N] mu;                     // linear predictor
  matrix<lower=0>[S,S] alpha;
  matrix<lower=0>[S,S] tau; 
  
  alpha = to_matrix(alpha_vec, S, S ); 
  tau = to_matrix(tau_vec, S,S); 
  
  for ( i in 1:N) {
    vector[S] C;         // transformed competition 
    for( j in 1:S){   
      C[j] = (alpha[focal[i], j]*x[i,j])^tau[focal[i], j]; 
    }
    mu[i] = lambda[focal[i]]/(1 + sum(C));  // matrix mult alpha by x vector 
  }
}
model {
  // priors
  lambda ~ cauchy(0, 5); 
  sigma ~ cauchy(0,2);
  alpha_vec ~ cauchy(0,1); 
  tau_vec ~ cauchy(0,1);
  
  // likelihood
  y ~ normal( mu , sigma);
}
generated quantities{ 
  vector[N] y_hat; 
  
  for(i in 1:N)
    y_hat[i] = normal_rng( mu[i], sigma);
}
