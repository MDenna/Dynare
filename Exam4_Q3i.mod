// Matthew Dennahower
// Monetary Theory - Exam 4
// Question 3 - Part I

// Endogenous Variables
var pi ygap y yn a rn i;

// Exogenous Variables
varexo ea;

// Parameters
parameters alpha beta sigma phi eps theta rhoa fipi fiy sigmaa h X omega bigt lambda gammaf gammab;

alpha = 0.33;
beta = 0.99;
sigma = 1;
phi = 1;
eps = 6;
theta = 0.75;
rhoa = 0.95;
//Central Bank follows Taylor interest rule
fipi = 1.5;
fiy = 0.125; 
sigmaa = 0.5;
h = 0.5;
omega = 0.5;
X = 1/(1+h);
//Internal Parameters
bigt = (1-alpha)/(1-alpha+alpha*eps);
lambda = ((1-theta)*(1-beta*theta))/(theta*(1+theta*beta*omega))*bigt;
gammaf = beta/(1+theta*beta*omega);
gammab = omega/(1+theta*beta*omega);

// Equations
model;

pi = gammaf*pi(+1) + gammab*pi(-1) + lambda/((1-h)*(1-alpha))*((sigma*(1-alpha)+(alpha+phi)*(1-h))/((1-h)*(1-alpha))*ygap - (sigma*h/(1-h))*ygap(-1));
ygap = y - yn;
yn = (sigma*h*(1-alpha))/(sigma*(1-alpha)+(alpha+phi)*(1-h))*yn(-1) - ((1-h)*(1+phi))/(sigma*(1-alpha)+(alpha+phi)*(1-h))*a;
a = rhoa*a(-1)+ea;
ygap = (1-X)*ygap(-1) + X*ygap(+1)-(1/sigma)*((1-h)/(1+h))*(i-pi(+1)-rn);
rn = (sigma/(1-h))*(yn(+1)-yn)- (sigma*h)/(1-h)*(yn-yn(-1));
// I) Taylor Rule 
i = fipi*pi + fiy*ygap;

end;

// Steady State
steady;

// Blanchard-Kahn conditions
check;

// Perturbation analysis
shocks;
// A Technological shock occurs
var ea; stderr sigmaa;
end;

// Stochastic simulation
stoch_simul (order=1, irf=100, periods = 250);
