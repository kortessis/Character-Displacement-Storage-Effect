clear
clc

% This script simulates evolution according to adaptive dynamics and then
% uses the output to calculate components relevant to species coexistence.

% The script has 4 sections that change what parameters are varied.
% Section 1. Vary species average log(y) for different sigmaV
% Section 2. Vary thetav1 - thetav2 for different rho.
% Section 3. Vary ydiff.
% Section 4. Vary rho.

% Which sections you want to run can be changed with the object RUN. RUN is
% a 4 element vector and takes values 1 and zero. Each element with the
% value 1 means that the corresponding section will be run. Thus, RUN =
% ones(1,4) runs the entire script. RUN = [1, 0, 1, 0] means that section 1
% and 3 are run.

% Note that we calculate A = 1/2*(rbar1/beta1 + rbar2/beta2) and SpecAveFit
% = rbar1/beta1 - rbar2/beta2 from the function EcologicalDynamics.m

% To determine coexistence, we use these two parameters to recover
% rbar1/beta1 and rbar2/beta2 because
%   rbar1/beta1 = 1/2*SpeAveFit + A
%   rbar2/beta2 = A - 1/2*SpecAveFit
RUN = [0, 0, 0, 1];

burnin = 1000;
gen = 10000 + burnin; 
mutnum = 3000;

s = 0.9;
alpha = 1;
muG = 0;
muV = 0;


%%
if RUN(1) == 1
    % Specific parameters for the sim of sigmaV and ybar
    ybar = linspace(log(5), log(60), 15);
    ydiff = 0;
    sigmaG = sqrt(1);
    sigmaV = sqrt([0.1, 0.5, 1, 2]);
    rho = 1;
    deltav = 0;
    
    % Values of thetag and thetav
    thetav = 0.5*deltav*[1, -1];
    deltaginit = 0.2*pi;
    thetaginit = 0.5*deltaginit*[1, -1];
    
    % Vectors to hold simulation outputs
    deltagSim1 = zeros(length(sigmaV), length(ybar));
    thetag1Sim1 = zeros(length(sigmaV), length(ybar));
    thetag2Sim1 = zeros(length(sigmaV), length(ybar));
    A1 = zeros(length(sigmaV), length(ybar));
    SpecAveFit1 = zeros(length(sigmaV), length(ybar));
    
    % Plotting infomation
    c = viridis(length(sigmaV) + 2);
    legendstring = {['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(1)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(2)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(3)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(4)^2)]};
    
    
    for i = 1:length(sigmaV)
        for j = 1:length(ybar)
            
            % Setting the value of y and other parameters
            y = exp(ybar(j)+0.5*ydiff*[1, -1]);
            Parameters = {s, y, rho, muG, muV, sigmaG, sigmaV(i), alpha};
            
            % Simulation finding trait evolution
            [thetag1, thetag2, MeanN] = AdaptDynamics22_PolarCoord(...
                thetaginit, thetav, gen, mutnum, Parameters);
            deltagvalues = thetag1 - thetag2;
            
            % Taking average of the trait values over the final 50
            % mutations
            thetag1Sim1(i, j) = mean(thetag1(mutnum - 50:mutnum));
            thetag2Sim1(i, j) = mean(thetag2(mutnum - 50:mutnum));
            deltagSim1(i,j) = mean(deltagvalues(mutnum-50:mutnum));
            
            
            % Calculate community stabilizing effect and species average fitnesses
            thetagSim = [thetag1Sim1(i, j), thetag2Sim1(i, j)];
            [A1(i, j), SpecAveFit1(i, j), ~, ~, ~, ~]...
                = Ecological_Dynamics(thetagSim, thetav, gen, burnin, Parameters);
            
            % Ticker to show how far along the simulation is
            [i/length(sigmaV), j/length(ybar)]
        end
    end

    %% Figure Showing Final Trait Differences
    figure()
    plot([min(ybar), max(ybar)], [0, pi], 'Color', 'none', 'HandleVisibility', 'off')
    h = gca;
    h.YTick = [0,pi/2,pi];
    h.YTickLabel = {0,'\pi/2','\pi'};
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('Per-biomass seed yield, ln$y$', 'FontSize', 30, 'Interpreter', 'Latex')
    ylabel('Evolved $\theta_{G_1} - \theta_{G_2}$', 'FontSize', 30, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(ybar, abs(deltagSim1(i,:)));
        set(p, {'LineWidth', 'Marker', 'Color', 'MarkerFaceColor', 'MarkerSize'}, ...
            {4,           'o',      c(i+1,:), c(i+1,:),          10});
    end
    hold off
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
            
    %% Figure Showing Community Average Mechanisms and Species Average Fitness Differences
    fig = figure();
    set(fig, 'defaultAxesColorOrder', zeros(2,3));
    
    yyaxis left
    plot([min(ybar), max(ybar)], [0, 0], 'Color', 'none', 'HandleVisibility', 'off')
    axis([min(ybar), max(ybar), min(min([A1,SpecAveFit1])), max(max([A1,SpecAveFit1]))])

    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('Per-biomass seed yield, ln$y$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('$A$', 'FontSize', 35, 'Interpreter', 'Latex')
    yyaxis right
    ylabel('$|\kappa_1 - \kappa_2|$', 'FontSize', 35, 'Interpreter', 'Latex')
    axis([min(ybar), max(ybar), min(min([A1,SpecAveFit1])), max(max([A1,SpecAveFit1]))])
    hold on
    
    for i = 1:length(sigmaV)
        p1 = plot(ybar, A1(i,:), '-', 'Color', c(i+1,:));
        p2 = plot(ybar, abs(SpecAveFit1(i, :)), ':', 'Color', c(i+1,:), 'HandleVisibility', 'off');
        set([p1, p2], {'LineWidth', 'MarkerSize', 'Marker'},...
                        {5, 10, 'o'});
        p1.MarkerFaceColor = c(i+1,:);
        p2.MarkerFaceColor = 'none';
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
end







%%
if RUN(2) == 1
    %% NOW SECTION ON EFFECT OF DIFFERENCES IN THETAV
    
    % Specific parameters for the sim of rho and thetav
    ybar = log(4);
    ydiff = 0;
    y = exp(ybar + 0.5*ydiff*[1, -1]);
    sigmaV = sqrt(1);
    sigmaG = sqrt(1);
    rho = [0,0.25,0.5,0.75,1];
    deltav = linspace(0, pi, 15);
    
    % Values of thetag and thetav
    thetav = 0.5*[1, -1]'*deltav;
    deltaginit = 0.2*pi;
    thetaginit = 0.5*deltaginit*[1, -1];
    
    % To hold evolutionary outcomes in sympatry
    deltagSim2 = zeros(length(rho), length(deltav));
    thetag1Sim2 = zeros(length(rho), length(deltav));
    thetag2Sim2 = zeros(length(rho), length(deltav));

    % To hold ecological outcomes after evolution in sympatry
    A2 = zeros(length(rho), length(deltav));
    SpecAveFit2 = zeros(length(rho), length(deltav));
    
    % To hold ecological outcomes after evolution in allopatry
    A2_allo = zeros(length(rho), length(deltav));
    SpecAveFit2_allo = zeros(length(rho), length(deltav));
    
    % Plotting code
    c = viridis(length(rho) + 2);
    legendstring = {'Allopatry',...
        ['$\rho =$ ', num2str(rho(1))],...
        ['$\rho =$ ', num2str(rho(2))],...
        ['$\rho =$ ', num2str(rho(3))],...
        ['$\rho =$ ', num2str(rho(4))],...
        ['$\rho =$ ', num2str(rho(5))]};
    legendstringalt = {['$\rho =$ ', num2str(rho(1))],...
        ['$\rho =$ ', num2str(rho(2))],...
        ['$\rho =$ ', num2str(rho(3))],...
        ['$\rho =$ ', num2str(rho(4))],...
        ['$\rho =$ ', num2str(rho(5))]};
    
    
    % Loop over all values of rho
    for i = 1:length(rho)
        
        % Set array of parameters
        Parameters = {s, y, rho(i), muG, muV, sigmaG, sigmaV, alpha};
        
        % Loop over all values of deltav
        for j = 1:length(deltav)
            [thetag1, thetag2, MeanN] = AdaptDynamics22_PolarCoord(...
                thetaginit, thetav(:, j)', gen, mutnum, Parameters);
            
            deltagvalues = thetag1 - thetag2;
            thetag1Sim2(i, j) = mean(thetag1(mutnum - 50:mutnum));
            thetag2Sim2(i, j) = mean(thetag2(mutnum - 50:mutnum));
            deltagSim2(i, j) = mean(deltagvalues(mutnum - 50:mutnum));
            
            % Calculate Abar and species average fitnesses in sympatry
            thetagSim = [thetag1Sim2(i, j), thetag2Sim2(i, j)];
            
            [A2(i, j), SpecAveFit2(i, j), ~,~,~,~]...
                = Ecological_Dynamics(thetagSim, thetav(:, j)', gen, burnin, Parameters);
            
            % Calculate Abar and Delta kappas in allopatry
            [A2_allo(i, j), SpecAveFit2_allo(i, j), ~,~,~,~]...
                = Ecological_Dynamics(thetav(:,j)', thetav(:, j)', gen, burnin, Parameters);
            
            [i/length(rho), j/length(deltav)]
        end
    end
    
    coex2 = A2 > abs(SpecAveFit2);
    mincoex = min(min([A2, SpecAveFit2]));
    maxcoex = max(max([A2, SpecAveFit2]));

    %% Plot the output: Evolved Trait Differences
    fig = figure();
    subplot(1,2,1)
    plot([0, pi], [0, pi], '--', 'Color', 'black')
    h = gca;
    h.YTick = [0, pi/2, pi];
    h.YTickLabel = {'0', '\pi/2', '\pi'};
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [0, pi/2, pi];
    h.XTickLabel = {'0', '\pi/2', '\pi'};
    xlabel('$\theta_{V_1} - \theta_{V_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Sympatric ESC $\theta_{G_1} - \theta_{G_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    legend()
    hold on
    
    for i = 1:length(rho)
        p = plot(deltav, abs(deltagSim2(i,:)), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    %% Plot Stabilizing mechanisms and average fitness differences

    set(fig, 'defaultAxesColorOrder', zeros(2,3));
    subplot(1,2,2)
    yyaxis left
    plot([0, pi], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    axis([0, pi, 0, max(max(A2 - A2_allo))])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [0, pi/2, pi];
    h.XTickLabel = {'0', '\pi/2', '\pi'};
    xlabel('$\theta_{V_1} - \theta_{V_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Sympatric - Allopatric $\bar{A}$', 'FontSize', 35, 'Interpreter', 'Latex')
    
    yyaxis right
    ylabel(["Sympatric - Allopatric"; "$\kappa_1 - \kappa_2$"], 'FontSize', 35, 'Interpreter', 'Latex')
    axis([0, pi, 0, max(max(A2 - A2_allo))])
    hold on
    
    for i = 1:length(rho)
        p1 = plot(deltav, A2(i,:) - A2_allo(i,:), '-', 'Color', c(i+1,:));
        p2 = plot(deltav, abs(SpecAveFit2(i, :) - SpecAveFit2_allo(i, :)), ':', 'Color', c(i+1,:), 'HandleVisibility', 'off');
        set([p1, p2], {'LineWidth', 'MarkerSize', 'Marker'},...
                        {5, 10, 'o'});
        p1.MarkerFaceColor = c(i+1,:);
        p2.MarkerFaceColor = 'none';
    end
    hold off
    
end




%%
if RUN(3) == 1
    %% EFFECT OF DIFFERENCES IN Y
    % Specific parameters for the sim of rho and thetav
    ybar = log(linspace(4, 45, 5)); % Was 6 to 60 with 5.
    ydiff = linspace(0, 0.25, 20); % WAS 0 to 0.2, with 10 values.
    sigmaG = sqrt(1);
    sigmaV = sqrt(1);
    rho = 1;
    deltav = 0;
    
    % Values of thetag and thetav
    thetav = 0.5*[1, -1]*deltav;
    deltaginit = 0.75*pi;
    thetaginit = 0.5*deltaginit*[1, -1];
    
    % Simulation output for sympatric ESC
    deltagSim3 = zeros(length(ybar), length(ydiff));
    thetag1Sim3 = zeros(length(ybar), length(ydiff));
    thetag2Sim3 = zeros(length(ybar), length(ydiff));
    
    % Ecological Output for sympatric ESC
    A3 = zeros(length(ybar), length(ydiff));
    SpecAveFit3 = zeros(length(ybar), length(ydiff));
    
    % Ecological Output for allopatric ESS
    A3_allo = zeros(length(ybar), length(ydiff));
    SpecAveFit3_allo = zeros(length(ybar), length(ydiff));
        
    c = viridis(length(ybar) + 2);
    legendstring = {['$\bar{y} =$ ', num2str(exp(ybar(1)))],...
        ['$\bar{y} =$ ', num2str(exp(ybar(2)))],...
        ['$\bar{y} =$ ', num2str(exp(ybar(3)))],...
        ['$\bar{y} =$ ', num2str(exp(ybar(4)))],...
        ['$\bar{y} =$ ', num2str(exp(ybar(5)))]};
    
    for i = 1:length(ybar)
        y = exp(ybar(i) + 0.5*ydiff'*[1, -1]);
        
        for j = 1:length(ydiff)
             Parameters = {s, y(j, :), rho, muG, muV, sigmaG, sigmaV, alpha};

             [thetag1, thetag2, MeanN] = AdaptDynamics22_PolarCoord(...
                thetaginit, thetav, gen, mutnum, Parameters);
            
            deltagvalues = thetag1 - thetag2;
            thetag1Sim3(i, j) = mean(thetag1(mutnum - 50:mutnum));
            thetag2Sim3(i, j) = mean(thetag2(mutnum - 50:mutnum));
            deltagSim3(i, j) = mean(deltagvalues(mutnum - 50:mutnum));
            
            %Calculate community stabilizing effect and species average fitnesses
            thetagSim3 = [thetag1Sim3(i, j), thetag2Sim3(i, j)];
            
            [A3(i, j), SpecAveFit3(i, j),~,~,~,~]...
                = Ecological_Dynamics(thetagSim3, thetav, gen, burnin, Parameters);
            
            %Calculate A and fitness differences using allopatric ESSs
            [A3_allo(i,j), SpecAveFit3_allo(i,j),~,~,~,~]... 
                = Ecological_Dynamics([0,0], [0,0], gen, burnin, Parameters);
            
            [i/length(ybar), j/length(ydiff)]
        end
    end
    
    
    coex3 = A3 > abs(SpecAveFit3);
    deltagSimCoex3 = NaN(length(ybar), length(ydiff));
    deltagSimCoex3(coex3) = deltagSim3(coex3);
    mincoex = min(min([abs(SpecAveFit3), A3]));
    maxcoex = max(max([abs(SpecAveFit3), A3]));
    
    %% Plot output
    figure()
    subplot(1,3,1)
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'none', 'HandleVisibility', 'off')
    h = gca;
    h.YTick = [0, pi/2, pi];
    h.YTickLabel = {'0', '\pi/2', '\pi'};
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Sympatric ESC $|\theta_{G_1} - \theta_{G_2}|$', 'FontSize', 35, 'Interpreter', 'Latex')
    legend()
    hold on
    
    for i = 1:length(ybar)
        p1 = plot(ydiff, deltagSim3(i,:), 'Color', c(i+1,:));
        p2 = plot(ydiff, deltagSimCoex3(i,:), 'Color', c(i+1,:));
        set([p1,p2], {'LineWidth', 'MarkerSize', 'Marker'}, {5, 20, 'o'})
        p1.HandleVisibility = 'off';
        p1.MarkerFaceColor = 'white';
        p2.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    subplot(1, 3, 2)
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('$\bar{A}_{Sympatric} - \bar{A}_{Allopatric}$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(ybar)
        p = plot(ydiff, A3(i,:) - A3_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    axis([ydiff(1), ydiff(end), 0, 1.1])
    hold off
    
    subplot(1, 3, 3)
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(["Sympatric - Allopatric"; "$\kappa'_1 - \kappa'_2$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(ybar)
        p = plot(ydiff, abs(SpecAveFit3(i,:) - SpecAveFit3_allo(i,:)), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    axis([ydiff(1), ydiff(end), 0, 1.1])
    hold off    
    
    
end


%%
if RUN(4) == 1
    %% EFFECT OF RHO
    % Specific parameters for the sim of rho and thetav
    ybar = log(4); % Was 6 to 60 with 5.
    ydiff = 0; % WAS 0 to 0.2, with 10 values.
    sigmaG = sqrt(1);
    sigmaV = sqrt(linspace(0.5, 2, 5));
    rho = linspace(0, 1, 10);
    deltav = 0;
    y = exp(ybar + 0.5*ydiff'*[1, -1]);

    % Values of thetag and thetav
    thetav = 0.5*[1, -1]*deltav;
    deltaginit = 0.75*pi;
    thetaginit = 0.5*deltaginit*[1, -1];
    
    deltagSim4 = zeros(length(sigmaV), length(rho));
    thetag1Sim4 = zeros(length(sigmaV), length(rho));
    thetag2Sim4 = zeros(length(sigmaV), length(rho));
    A4 = zeros(length(sigmaV), length(rho));
    SpecAveFit4 = zeros(length(sigmaV), length(rho));
    
    c = viridis(length(sigmaV) + 2);
    legendstring = {['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(1)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(2)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(3)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(4)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(5)^2)]};
    
    for i = 1:length(sigmaV)
        
        for j = 1:length(rho)
             Parameters = {s, y, rho(j), muG, muV, sigmaG, sigmaV(i), alpha};

             [thetag1, thetag2, MeanN] = AdaptDynamics22_PolarCoord(...
                thetaginit, thetav, gen, mutnum, Parameters);
            
            deltagvalues = thetag1 - thetag2;
            thetag1Sim4(i, j) = mean(thetag1(mutnum - 50:mutnum));
            thetag2Sim4(i, j) = mean(thetag2(mutnum - 50:mutnum));
            deltagSim4(i, j) = mean(deltagvalues(mutnum - 50:mutnum));
            
            %Calculate community stabilizing effect and species average fitnesses
            thetagSim4 = [thetag1Sim4(i, j), thetag2Sim4(i, j)];
            
            [A4(i, j), SpecAveFit4(i, j), ~, ~, ~, ~]...
                = Ecological_Dynamics(thetagSim4, thetav, gen, burnin, Parameters);
            
            [i/length(sigmaV), j/length(rho)]
        end
    end
    
    
    coex4 = A4 > abs(SpecAveFit4);
    deltagSimCoex4 = NaN(length(sigmaV), length(rho));
    deltagSimCoex4(coex4) = deltagSim4(coex4);
    mincoex = min(min([abs(SpecAveFit4), A4]));
    maxcoex = max(max([abs(SpecAveFit4), A4]));
    
    %% Evolved Trait Differences Figure
    figure()
    plot(rho, zeros(1,length(rho)), '--', 'Color', 'none', 'HandleVisibility', 'off')
    h = gca;
    h.YTick = [0, pi/2, pi];
    h.YTickLabel = {'0', '\pi/2', '\pi'};
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('$\rho$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Evolved $\theta_{G_1} - \theta_{G_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    legend()
    hold on
    
    for i = 1:length(sigmaV)
        p1 = plot(rho, deltagSim4(i,:), 'Color', c(i+1,:));
        p2 = plot(rho, deltagSimCoex4(i,:), 'Color', c(i+1,:));
        set([p1,p2], {'LineWidth', 'MarkerSize', 'Marker'}, {5, 10, 'o'})
        p1.HandleVisibility = 'off';
        p1.MarkerFaceColor = 'white';
        p2.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
       
    %% Full stabilizing mechanisms and average fitness differences
    fig = figure();
    set(fig, 'defaultAxesColorOrder', zeros(2,3));
    
    yyaxis left
    plot([0, max(rho)], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    axis([0, max(rho), mincoex, maxcoex])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('$\rho$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('$A$', 'FontSize', 35, 'Interpreter', 'Latex')
    
    yyaxis right
    ylabel('$|\kappa_1 - \kappa_2|$', 'FontSize', 35, 'Interpreter', 'Latex')
    axis([0, max(rho), mincoex, maxcoex])
    hold on
    
    for i = 1:length(sigmaV)
        p1 = plot(rho, A4(i,:), '-', 'Color', c(i+1,:));
        p2 = plot(rho, abs(SpecAveFit4(i, :)), ':', 'Color', c(i+1,:), 'HandleVisibility', 'off');
        set([p1, p2], {'LineWidth', 'MarkerSize', 'Marker'},...
                        {5, 10, 'o'});
        p1.MarkerFaceColor = c(i+1,:);
        p2.MarkerFaceColor = 'none';
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
end