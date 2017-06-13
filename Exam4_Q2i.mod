// Matthew Dennahower
// 103728472
// Monetary Theory - Exam 4
// Question 2 - Part I

// Endogenous Variables
var pi ygap y yn a rn u i;
 
// Exogenous Variables
varexo ea eu;
 
// Parameters
parameters alpha beta sigma phi eps theta rhoa rhou sigmaa sigmau bigt lambda psi K alphax;
 
// Calibration
alpha = 0.33;
beta = 0.99;
sigma = 1;
phi = 1;
eps = 6;
theta = 0.75;
rhoa = 0.95;
rhou = 0.75;
sigmaa = 0.75;
sigmau = 0.5;
bigt = (1-alpha)/(1-alpha+alpha*eps);
lambda = ((1-theta)*(1-beta*theta)/theta)*bigt;
psi = (1+phi)/(sigma*(1-alpha)+alpha+phi);
K = lambda*(sigma + (alpha+phi)/(1-alpha));
alphax = K/eps;
 
// Equations
model;
 
pi = beta*pi(+1) + K*ygap + u;
ygap = y - yn;
yn = psi*a;
ygap = ygap(+1)-(1/sigma)*(i-pi(+1)-rn);
rn = -sigma*psi*(1-rhoa)*a;
a = rhoa*a(-1) + ea;
u = rhou*u(-1) + eu; 

// CB follows optimal commitment policy
pi = -(alphax/K)*ygap + (alphax/K)*ygap(-1);

// CB follows optimal discretion policy
//i = gammapi*alphax*q*rhou*u + sigma*rn;
end;

 
// Steady State
steady;
 
// Blanchard-Kahn conditions
check;
 
// Perturbation analysis
shocks;
var eu; stderr sigmau;
end;
 
// Stochastic simulation
stoch_simul (order=1, irf=100, periods = 250);
