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
  real<lower=0> sigma;              // std dev 
  vector<lower=0>[S] tau;
}
transformed parameters{
  vector[N] mu;                     // linear predictor
  matrix<lower=0>[S,S] alpha;
  
  alpha = to_matrix(alpha_vec, S, S ); 
  
  for ( i in 1:N) {
    mu[i] = lambda[focal[i]]/(1 + alpha[focal[i],1:S]*x[i] )^tau[focal[i]];  // matrix mult alpha by x vector 
  }
}
model {
  // priors
  lambda ~ cauchy(0, 5); 
  sigma ~ cauchy(0,2);
  alpha_vec ~ cauchy(0,1); 
  tau ~ cauchy(0,1);
  // likelihood
  y ~ normal( mu , sigma);
}
generated quantities{ 
  vector[N] y_hat; 
  
  for(i in 1:N)
    y_hat[i] = normal_rng( mu[i], sigma);
}
