function [A, SpecAveFit, kappaDiff, ComAveDeltaJ, ComAveDeltaIG, ComAveDeltaIV]...
    = Ecological_Dynamics(thetag, thetav, gen, burnin, Parameters)

% This function finds the relevant ecological quantities for understanding species coexistence in the model. In effect, it takes parameters and the number of time steps to run dynamics and then calculates exact and approximate species coexistence mechanisms.

%% Inputs
% thetag = a 2d vector holding the values thetag1 and thetag2
% thetav = a 2d vector holding the values thetav1 and thetav2
% gen = the length of year over which to evaluate the long-term growth rate
% burnin = the amount of years to run the model to get rid of transient dynamics

% Parameters: a cell array containing all other model parameters in the order

s = Parameters{1};
y = Parameters{2};
rho = Parameters{3};
muG = Parameters{4};
muV = Parameters{5};
sigmaG = Parameters{6};
sigmaV = Parameters{7};
alpha = Parameters{8};

%% Outputs
% A = exact (within computational error) community average stabilizing mechansims.
% SpecAveFit = \kappa_1\prime - \kappa_2\prime, the exact adjusted species average fitness Differences
% kappaDiff = fitness differences using expressions from the small variance approximation
% ComAveDeltaJ = community average relative nonlinearity using expressions from the small varaince approximation
% ComAveDeltaIG = community average storage effect from germination using expressions from the small variance approximation
% ComAveDeltaIV = community average storage effect from vigor using expressions from the small variance approximation

% Generate the distribution of environmental factors (X1,X2,Z1,Z2) for gen time steps
XZ = mvnrnd(zeros(1,4), [eye(2), rho*eye(2); rho*eye(2), eye(2)], gen);
X = XZ(:, 1:2);
Z = XZ(:, 3:4);

% Generate a sequence of environemntal responses for G and V for each species
EG = muG + sigmaG*[sin(thetag'), cos(thetag')]*X';
G = exp(EG)./(1 + exp(EG));
EV = muV + sigmaV*[sin(thetav'), cos(thetav')]*Z';
V = exp(EV);

% Create vecotrs to hold the values for resident population dynamics and competition created by each species.
Nres = zeros(2, gen);
Cres = zeros(2, gen);
Nres(:, 1) = 10;
Cres(:, 1) = log(1 + alpha*G(:, 1).*V(:, 1).*Nres(:, 1));

% Run the dynamics for each species as resident (i.e., without the other species)
for t = 2:gen
    Nres(:, t) = Nres(:, t-1).*(s*(1 - G(:, t-1)) + y'.*G(:, t-1).*V(:, t-1)./exp(Cres(:, t-1)));
    Cres(:, t) = log(1 + alpha*G(:, t).*V(:, t).*Nres(:, t));
end

% Calculate the invader growth rates for both species
inv1 = mean(log(s*(1 - G(1, burnin:end)) + y(1)*G(1, burnin:end).*V(1, burnin:end)./exp(Cres(2, burnin:end))), 2);
inv2 = mean(log(s*(1 - G(2, burnin:end)) + y(2)*G(2, burnin:end).*V(2, burnin:end)./exp(Cres(1, burnin:end))), 2);

% Calculate sundry values used in the low density approximation for the approximate mechanisms equations
Gbarstar = exp(muG)/(1 + exp(muG));
Vbarstar = exp(muV);
Cbarstar = log(y*Gbarstar*Vbarstar*(1-Gbarstar)/(1 - s*(1 - Gbarstar)));
Cstar = mean(Cbarstar);
Estar = log((1 - s)./(y*exp(muV - Cstar) - 1));
Gstar = exp(Estar)./(1+exp(Estar));
beta = 1 - s*(1 - Gstar);

% Calculate all possible covariances between EGs, EVs, and C.
CovEGEVCres = cov([EG(:,burnin:end); EV(:, burnin:end); Cres(:, burnin:end)]');

% Calculate (approximate) Fitness Differences
kappaDiff = log(y(1)/y(2)) + ...
    0.5*(beta(2)-beta(1))*((1-s)*sigmaG.^2/(beta(1)*beta(2))+sigmaV^2) +...
    sigmaG*sigmaV*rho*((1-beta(1))*cos(thetag(1)-thetav(1)) - (1-beta(2))*cos(thetag(2)-thetav(2)));

% Calculate (approximate) Relative Nonlinearity (Delta J)
DeltaJ(1) = 0.5*(beta(1) - beta(2))*CovEGEVCres(6,6);
DeltaJ(2) = 0.5*(beta(2) - beta(1))*CovEGEVCres(5,5);
ComAveDeltaJ = 0.5*sum(DeltaJ);

% Calculate (approximate) Storage Effect From Germination (Delta IG)
DeltaIG(1) = (1-beta(2))*CovEGEVCres(2,6) - (1-beta(1))*CovEGEVCres(1,6);
DeltaIG(2) = (1-beta(1))*CovEGEVCres(1,5) - (1-beta(2))*CovEGEVCres(2,5);
ComAveDeltaIG = 0.5*sum(DeltaIG);

% Calculate (approximate) Storage Effect From Vigor (Delta IV)
DeltaIV(1) = (1-beta(2))*CovEGEVCres(4,6) - (1-beta(1))*CovEGEVCres(3,6);
DeltaIV(2) = (1-beta(1))*CovEGEVCres(3,5) - (1-beta(2))*CovEGEVCres(4,5);
ComAveDeltaIV = 0.5*sum(DeltaIV);

% Calculate (exact) community average stabilizing mechanisms and adjusted species average fitness difference.
A = mean([inv1/beta(1), inv2/beta(2)]);
SpecAveFit = (inv1/beta(1) - inv2/beta(2))/2;
