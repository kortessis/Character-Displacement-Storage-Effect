function[thetag1, thetag2, MeanN] = AdaptDynamics22_PolarCoord(thetaginit, thetav, gen, mutnum, Parameters)

% OUTPUTS
% thetag1 = vector holding values of theta_G1(t)
% thetag1 = vector holding values of theta_G2(t)
% MeanN = long-term average of N1 and N2 for each trait value

% deltaginit = the initial difference in thetag values between species
% deltav = the initial difference in thetav values between species
% gen = the number of generations beyond the burnin to evaluate rbar.
% mutnum = the total number of mutations introduced into the system.
%    mutnum/2 gives the number per species.

% Inputs are in the form of a cell array with elements
% s = seed survival upon dormancy
% rho = maximum possible level of predictive germination
% y = [y1, y2], a 1x2 vector for seed yield per unit biomass
% sigmaG = stdev(EG) = stdev(ln(G/(1-G)))
% sigmaV = stdev(EV) = stdev(ln(V))

% Load all the parameters
s = Parameters{1};
y = Parameters{2};
rho = Parameters{3};
muG = Parameters{4};
muV = Parameters{5};
sigmaG = Parameters{6};
sigmaV = Parameters{7};
alpha = Parameters{8};
% Assume mutation size maximum is 0.01
mut_size = 0.01;
% Set the time for a resident to reach stationary state as burnin years
burnin = 100;

% Calculate G*, V*, and N* for the initial conditions.
Gbarstar = exp(muG)./(1+exp(muG));
Vbarstar = exp(muV);
Nbarstar = (y*Gbarstar*Vbarstar/(1 - s*(1-Gbarstar)) - 1)/(alpha*Gbarstar*Vbarstar);

% Create two vectors to hold the valuues of species traits. One for species 1 and one for species 2.
trait_state1 = zeros(1, mutnum+1);
trait_state2 = zeros(1, mutnum+1);

% Load the values of
thetav1 = thetav(1);
thetav2 = thetav(2);

% Create the covariance structure for the environmental variables [X1,X2,Z1,Z2]
SIGMA = [eye(2), rho*eye(2); rho*eye(2), eye(2)];

% Need a vector with steps of evolutionary time to hold resident trait values.
trait_state1(1) = thetaginit(1);
trait_state2(1) = thetaginit(2);
meanN1 = zeros(1,mutnum);
meanN2 = zeros(1,mutnum);

% Now simulate over all numbers of mutations
for j = 1:mutnum

% Generate the distribution [X1, X2, Z1, Z2]
    X = mvnrnd(zeros(1,4), SIGMA, gen+burnin);
    Z = X(:, 3:4);
    X = X(:, 1:2);
    
    % Create the vigor and germination processes.
    V1 = exp(muV + sigmaV*([sin(thetav1), cos(thetav1)]*Z'));
    V2 = exp(muV + sigmaV*([sin(thetav2), cos(thetav2)]*Z'));
    
    EG1 = muG + sigmaG*([sin(trait_state1(j)), cos(trait_state1(j))]*X');
    EG2 = muG + sigmaG*([sin(trait_state2(j)), cos(trait_state2(j))]*X');
    G1 = exp(EG1)./(1+exp(EG1));
    G2 = exp(EG2)./(1+exp(EG2));

    % Create vectors to hold population dynamics and competition
    N1 = zeros(1,gen+burnin);
    N2 = zeros(1,gen+burnin);
    C = zeros(1,gen+burnin);
    
    % Set intial conditions
    N1(1) = Nbarstar(1)/2;
    N2(1) = Nbarstar(2)/2;
    C(1) = 1+alpha*G1(1)*V1(1)*N1(1) + alpha*G2(1)*V2(1)*N2(1);
    
    % Run ecological dynamics
    for time = 2:gen+burnin
        N1(time) = N1(time-1)*(s*(1-G1(time-1)) + y(1)*G1(time-1)*V1(time-1)/C(time-1));
        N2(time) = N2(time-1)*(s*(1-G2(time-1)) + y(2)*G2(time-1)*V2(time-1)/C(time-1));
        C(time) = 1 + alpha*G1(time)*V1(time)*N1(time) + alpha*G2(time)*V2(time)*N2(time);
    end
    
    % Randomly draw the size of the mutation
    mutation = 2*mut_size*(rand(1)-0.5);
    
    % Calculate the mutant phenotypes and ask if they can invade
    if mod(j,2)
        
        % Calculate mutant phenotype
        xinv = trait_state1(j) + mutation;
    
        if xinv > pi/2
            xinv = pi/2;
        else
            if xinv < -pi/2
                xinv = -pi/2;
            end
        end
        
        % Create mutant germination responses
        EG1inv = muG + sigmaG*([sin(xinv), cos(xinv)]*X');
        G1inv = exp(EG1inv)./(1+exp(EG1inv));        
        
        % Calculate resident and mutant long-term growth rates
        rbar1_wt = mean(log(s*(1-G1(burnin:end)) + y(1)*G1(burnin:end).*V1(burnin:end)./C(burnin:end)));
        rbar1_mut = mean(log(s*(1-G1inv(burnin:end)) + y(1)*G1inv(burnin:end).*V1(burnin:end)./C(burnin:end)));
        
        % Substitute mutant growth rates if they can invade the resident. Throw them out if not.
        if rbar1_mut > rbar1_wt
            trait_state1(j+1) = xinv;
        else
            trait_state1(j+1) = trait_state1(j);
        end
        
        trait_state2(j+1) = trait_state2(j);
        
    else % Do the same for species 2 if the mutation number is even
        xinv = trait_state2(j) + mutation;
    
        if xinv > pi/2
            xinv = pi/2;
        else
            if xinv < -pi/2
                xinv = -pi/2;
            end
        end
        
        EG2inv = muG + sigmaG*([sin(xinv), cos(xinv)]*X');
        G2inv = exp(EG2inv)./(1+exp(EG2inv));
        
        rbar2_wt = mean(log(s*(1-G2(burnin:end)) + y(2)*G2(burnin:end).*V2(burnin:end)./C(burnin:end)));
        rbar2_mut = mean(log(s*(1-G2inv(burnin:end)) + y(2)*G2inv(burnin:end).*V2(burnin:end)./C(burnin:end)));
        
        if rbar2_mut > rbar2_wt
            trait_state2(j+1) = xinv;
        else
            trait_state2(j+1) = trait_state2(j);
        end
        trait_state1(j+1) = trait_state1(j);
        
    end
    
    % Store average population densities
    meanN1(j) = mean(N1(burnin:end));
    meanN2(j) = mean(N2(burnin:end));
end

% Collect all output
MeanN = [meanN1; meanN2];
thetag1 = trait_state1;
thetag2 = trait_state2;


end
