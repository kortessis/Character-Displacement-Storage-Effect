clear
clc

% This script simulates evolution according to adaptive dynamics and then
% uses the output to calculate components relevant to species coexistence.

% The script has 4 sections that change what parameters are varied.
% Section 1. Vary species average log(y) for different sigmaV
% Section 2. Vary theta1 - theta2 for different rho.
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
RUN = [1, 0, 0, 0];

burnin = 1000;
gen = 50000 + burnin; 
mutnum = 3500; % Was 3000

s = 0.9;
alpha = 1;
muG = 0;
muV = 0;


%%
if RUN(1) == 1
    % Specific parameters for the sim of sigmaV and ybar
    ybar = log(linspace(5, 60, 15));
    ydiff = 0;
    sigmaG = sqrt(1);
    sigmaV = sqrt([0.1, 0.5, 1, 2]);
    rho = 1;
    deltav = 0;
    
    % Values of thetag and thetav
    thetav = 0.5*deltav*[1, -1];
    deltaginit = 0.2*pi;
    thetaginit = 0.5*deltaginit*[1, -1];
    
    % Simulation outputs
    deltagSim1 = zeros(length(sigmaV), length(ybar));
    thetag1Sim1 = zeros(length(sigmaV), length(ybar));
    thetag2Sim1 = zeros(length(sigmaV), length(ybar));
    A1 = zeros(length(sigmaV), length(ybar));
    SpecAveFit1 = zeros(length(sigmaV), length(ybar));
    kappaDiff1 = zeros(length(sigmaV), length(ybar));
    ComAveDeltaN1 = zeros(length(sigmaV), length(ybar));
    ComAveDeltaIG1 = zeros(length(sigmaV), length(ybar));
    ComAveDeltaIV1 = zeros(length(sigmaV), length(ybar));
    
    % Plotting values
    c = viridis(length(sigmaV) + 2);
    legendstring = {['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(1)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(2)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(3)^2)],...
        ['$\sigma_{E_V}^2 =$ ', num2str(sigmaV(4)^2)]};
    
    
    for i = 1:length(sigmaV)
        for j = 1:length(ybar)
            
            y = exp(ybar(j)+0.5*ydiff*[1, -1]);
            Parameters = {s, y, rho, muG, muV, sigmaG, sigmaV(i), alpha};
            
            % Simulation finding trait evolution
            [thetag1, thetag2, MeanN] = AdaptDynamics22_PolarCoord(...
                thetaginit, thetav, gen, mutnum, Parameters);
            deltagvalues = thetag1 - thetag2;
            thetag1Sim1(i, j) = mean(thetag1(mutnum - 50:mutnum));
            thetag2Sim1(i, j) = mean(thetag2(mutnum - 50:mutnum));
            deltagSim1(i,j) = mean(deltagvalues(mutnum-50:mutnum));
            
            
            % Calculate community stabilizing effect and species average fitnesses
            thetagSim = [thetag1Sim1(i, j), thetag2Sim1(i, j)];
            [A1(i, j), SpecAveFit1(i, j), kappaDiff1(i, j), ComAveDeltaN1(i, j), ComAveDeltaIG1(i, j), ComAveDeltaIV1(i, j)]...
                = EcologicalDynamics(thetagSim, thetav, gen, burnin, Parameters);
            
            % Ticker to show how far along the simulation is
            [i/length(sigmaV), j/length(ybar)]
        end
    end
    
    % Use simulation outputs to find the
    Aapprox1 = -ComAveDeltaN1 + ComAveDeltaIG1 + ComAveDeltaIV1;
    rbarapprox1 = kappaDiff1 + Aapprox1;
    kappadiffPrime1 = rbarapprox1 - Aapprox1;
    r11 = SpecAveFit1 + A1;
    r21 = A1 - SpecAveFit1;
    coex1 = r11 > 0 & r21 > 0;
    deltagSim1coex = NaN(length(sigmaV), length(ybar));
    deltagSim1coex(coex1) = deltagSim1(coex1);
    
    mincoex = min(min([Aapprox1, A1, SpecAveFit1, kappadiffPrime1]));
    maxcoex = max(max([Aapprox1, A1, SpecAveFit1, kappadiffPrime1]));
    
    %% Figure Showing Final Trait Differences
    figure()
    plot([min(exp(ybar)), max(exp(ybar))], [0, pi], 'Color', 'none', 'HandleVisibility', 'off')
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
    
    
    %% Stabilizing and Fitness Difference Figures with approx and full
    figure();
    subplot(2, 2, 1)
    
    plot([min(ybar), max(ybar)], [0, 0], 'Color', 'none', 'HandleVisibility', 'off')
    axis([min(exp(ybar)), max(exp(ybar)), mincoex, maxcoex])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('(ln$y_1 + $ln$y_2)/2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('$A$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p1 = plot(exp(ybar), A1(i,:), 'Color', c(i+1,:));
        set(p1, {'LineWidth', 'MarkerSize', 'Marker'},...
                        {5, 10, 'o'});
        p1.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(2, 2, 3)
    plot([min(ybar), max(ybar)], [0, 0], 'Color', 'none', 'HandleVisibility', 'off')
    axis([min(exp(ybar)), max(exp(ybar)), mincoex, maxcoex])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('(ln$y_1 + $ln$y_2)/2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('$\kappa_1 - \kappa_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(exp(ybar), SpecAveFit1(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    
    subplot(2, 2, 2)
    plot([min(ybar), max(ybar)], [0, 0], 'Color', 'none', 'HandleVisibility', 'off')
    axis([min(exp(ybar)), max(exp(ybar)), mincoex, maxcoex])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('(ln$y_1 + $ln$y_2)/2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Approx $A$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(exp(ybar), Aapprox1(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(2, 2, 4)
    plot([min(ybar), max(ybar)], [0, 0], 'Color', 'none', 'HandleVisibility', 'off')
    axis([min(exp(ybar)), max(exp(ybar)), mincoex, maxcoex])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('(ln$y_1 + $ln$y_2)/2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Approx $\kappa_1 - \kappa_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(exp(ybar), kappadiffPrime1(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    
    %% Figure Showing Community Average Mechanisms and Species Average Fitness Differences
    fig = figure();
    set(fig, 'defaultAxesColorOrder', zeros(2,3));
    
    yyaxis left
    plot([min(ybar), max(ybar)], [0, 0], 'Color', 'none', 'HandleVisibility', 'off')
    axis([min(ybar), max(ybar), min(min([A1,SpecAveFit1])), max(max([A1,SpecAveFit1]))])

    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('(ln$y_1$ + ln$y_2$)/2', 'FontSize', 35, 'Interpreter', 'Latex')
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
    
    
    %% Figure Showing Final Trait Values
    figure()
    surf(thetag1Sim1)
    hold on
    surf(thetag2Sim1)
    hold off
    xlabel('$\bar{y}$', 'Interpreter', 'Latex')
    ylabel('$\sigma_{E_V}$', 'Interpreter', 'Latex')
    zlabel('Evolved $\theta_G$', 'Interpreter', 'Latex')
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
    kappaDiff2 = zeros(length(rho), length(deltav));
    ComAveDeltaN2 = zeros(length(rho), length(deltav));
    ComAveDeltaIG2 = zeros(length(rho), length(deltav));
    ComAveDeltaIV2 = zeros(length(rho), length(deltav));

    % To hold ecological outcomes after evolution in allopatry
    A2_allo = zeros(length(rho), length(deltav));
    SpecAveFit2_allo = zeros(length(rho), length(deltav));
    kappaDiff2_allo = zeros(length(rho), length(deltav));
    ComAveDeltaN2_allo = zeros(length(rho), length(deltav));
    ComAveDeltaIG2_allo = zeros(length(rho), length(deltav));
    ComAveDeltaIV2_allo = zeros(length(rho), length(deltav));
    
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
    
    for i = 1:length(rho)
        Parameters = {s, y, rho(i), muG, muV, sigmaG, sigmaV, alpha};
        
        for j = 1:length(deltav)
            [thetag1, thetag2, MeanN] = AdaptDynamics22_PolarCoord(...
                thetaginit, thetav(:, j)', gen, mutnum, Parameters);
            
            deltagvalues = thetag1 - thetag2;
            thetag1Sim2(i, j) = mean(thetag1(mutnum - 50:mutnum));
            thetag2Sim2(i, j) = mean(thetag2(mutnum - 50:mutnum));
            deltagSim2(i, j) = mean(deltagvalues(mutnum - 50:mutnum));
            
            %Calculate Abar and species average fitnesses in sympatry
            thetagSim = [thetag1Sim2(i, j), thetag2Sim2(i, j)];
            
            [A2(i, j), SpecAveFit2(i, j), kappaDiff2(i, j), ComAveDeltaN2(i, j), ComAveDeltaIG2(i, j), ComAveDeltaIV2(i, j)]...
                = EcologicalDynamics(thetagSim, thetav(:, j)', gen, burnin, Parameters);
            
            % Calculate Abar and Delta kappas in allopatry
            [A2_allo(i, j), SpecAveFit2_allo(i, j), kappaDiff2_allo(i, j), ComAveDeltaN2_allo(i, j), ComAveDeltaIG2_allo(i, j), ComAveDeltaIV2_allo(i, j)]...
                = EcologicalDynamics(thetav(:,j)', thetav(:, j)', gen, burnin, Parameters);
            
            [i/length(rho), j/length(deltav)]
        end
    end
    
    Aapprox2 = -ComAveDeltaN2 + ComAveDeltaIG2 + ComAveDeltaIV2;
    rbarapprox2 = kappaDiff2 + Aapprox2;
    kappadiffPrime2 = rbarapprox2 - Aapprox2;
    
    Aapprox2_allo = -ComAveDeltaN2_allo + ComAveDeltaIG2_allo + ComAveDeltaIV2_allo;
    rbarapprox2_allo = kappaDiff2_allo + Aapprox2_allo;
    kappadiffPrime2_allo = rbarapprox2_allo - Aapprox2_allo;
    
    r12 = SpecAveFit2 + A2;
    r22 = A2 - SpecAveFit2;
    coex2 = r12 > 0 & r22 > 0;
    mincoex = min(min([Aapprox2, A2, kappadiffPrime2, SpecAveFit2]));
    maxcoex = max(max([Aapprox2, A2, kappadiffPrime2, SpecAveFit2]));

    %% Evolved Trait Differences Figure
    corrV = cos(deltav);
    corrG = cos(deltagSim2);
    figure()
    subplot(1, 2, 2)
    plot([-1, 1], [-1, 1], '--', 'Color', 'black')
    h = gca;
    h.YTick = [-1,0,1];
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [-1,0,1];
    xlabel('corr$(E_{V_1}, E_{V_2})$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel("Sympatric ESC corr$(E_{G_1}, E_{G_2})$", 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(rho)
        p = plot(corrV, corrG(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(1, 2, 1)
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
        p = plot(deltav, deltagSim2(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    %% Stabilizing and Fitness Difference Figures with approx and full
    figure()
    subplot(2, 2, 1)
    plot([0, pi], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [0, pi/2, pi];
    h.XTickLabel = {'0', '\pi/2', '\pi'};
    xlabel('$\theta_{V_1} - \theta_{V_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('$Sympatric - Allopatric \bar{A}$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(rho)
        p = plot(deltav, A2(i,:) - A2_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(2, 2, 3)
    plot([0, pi], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [0, pi/2, pi];
    h.XTickLabel = {'0', '\pi/2', '\pi'};
    xlabel('$\theta_{V_1} - \theta_{V_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(["Sympatric - Allopatric"; "$\kappa_1 - \kappa_2$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(rho)
        p = plot(deltav, SpecAveFit2(i,:) - SpecAveFit2_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstringalt, 'off'});
    
    subplot(2, 2, 2)
    plot([0, pi], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [0, pi/2, pi];
    h.XTickLabel = {'0', '\pi/2', '\pi'};
    xlabel('$\theta_{V_1} - \theta_{V_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(["Sympatric - Allopatric"; "$-\Delta N + \Delta I_G + \Delta I_V$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(rho)
        p = plot(deltav, Aapprox2(i,:) - Aapprox2_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(2, 2, 4)
    plot([0, pi], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [0, pi/2, pi];
    h.XTickLabel = {'0', '\pi/2', '\pi'};
    xlabel('$\theta_{V_1} - \theta_{V_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(["Sympatric - Allopatric"; "Approximte $\kappa_1 - \kappa_2$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(rho)
        p = plot(deltav, kappadiffPrime2(i,:) - kappadiffPrime2_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    %% Full stabilizing mechanisms and average fitness differences
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
    
    subplot(1,2,2)
    plot([0, pi], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    axis([0, pi, 0, max(max(A2 - A2_allo))])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    h.XTick = [0, pi/2, pi];
    h.XTickLabel = {'0', '\pi/2', '\pi'};
    xlabel('$\theta_{V_1} - \theta_{V_2}$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Sympatric - Allopatric $\bar{A}$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(rho)
        p1 = plot(deltav, A2(i,:) - A2_allo(i,:), '-', 'Color', c(i+1,:));
        set(p1, {'LineWidth', 'MarkerSize', 'Marker'},...
                        {5, 10, 'o'});
        p1.MarkerFaceColor = c(i+1,:);
    end
    hold off
        
    %% Figure showing endvalue of theta g for each species
    figure()
    surf(deltav, rho, thetag1Sim2)
    hold on
    surf(deltav, rho, thetag2Sim2)
    hold off
    xlabel('$\theta_{G_1} - \theta_{G_2}$', 'Interpreter', 'Latex')
    ylabel('$rho$', 'Interpreter', 'Latex')
    zlabel('Evolved $\theta_G$', 'Interpreter', 'Latex')
end




%%
if RUN(3) == 1
    %% NOW SECTION ON EFFECT OF DIFFERENCES IN Y
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
    kappaDiff3 = zeros(length(ybar), length(ydiff));
    ComAveDeltaN3 = zeros(length(ybar), length(ydiff));
    ComAveDeltaIG3 = zeros(length(ybar), length(ydiff));
    ComAveDeltaIV3 = zeros(length(ybar), length(ydiff));
    
    % Ecological Output for allopatric ESS
    A3_allo = zeros(length(ybar), length(ydiff));
    SpecAveFit3_allo = zeros(length(ybar), length(ydiff));
    kappaDiff3_allo = zeros(length(ybar), length(ydiff));
    ComAveDeltaN3_allo = zeros(length(ybar), length(ydiff));
    ComAveDeltaIG3_allo = zeros(length(ybar), length(ydiff));
    ComAveDeltaIV3_allo = zeros(length(ybar), length(ydiff));

    
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
            
            [A3(i, j), SpecAveFit3(i, j), kappaDiff3(i, j), ComAveDeltaN3(i, j), ComAveDeltaIG3(i, j), ComAveDeltaIV3(i, j)]...
                = EcologicalDynamics(thetagSim3, thetav, gen, burnin, Parameters);
            
            %Calculate A and fitness differences using allopatric ESSs
            [A3_allo(i,j), SpecAveFit3_allo(i,j), kappaDiff3_allo(i,j), ComAveDeltaN3_allo(i,j), ComAveDeltaIG3_allo(i,j), ComAveDeltaIV3_allo(i,j)]... 
                = EcologicalDynamics([0,0], [0,0], gen, burnin, Parameters);
            
            [i/length(ybar), j/length(ydiff)]
        end
    end
    
    Aapprox3 = -ComAveDeltaN3 + ComAveDeltaIG3 + ComAveDeltaIV3;
    rbarapprox3 = kappaDiff3 + Aapprox3;
    kappadiffPrime3 = rbarapprox3 - Aapprox3;
    
    Aapprox3_allo = -ComAveDeltaN3_allo + ComAveDeltaIG3_allo + ComAveDeltaIV3_allo;
    rbarapprox3_allo = kappaDiff3_allo - Aapprox3_allo;
    kappadiffPrime3_allo = rbarapprox3_allo - Aapprox3_allo;
    
    r13 = SpecAveFit3 + A3;
    r23 = A3 - SpecAveFit3;
    coex3 = r13 > 0 & r23 > 0;
    deltagSimCoex3 = NaN(length(ybar), length(ydiff));
    deltagSimCoex3(coex3) = deltagSim3(coex3);
    mincoex = min(min([abs(SpecAveFit3), A3, Aapprox3, kappadiffPrime3]));
    maxcoex = max(max([abs(SpecAveFit3), A3, Aapprox3, kappadiffPrime3]));
    
    %% Paper figure for this case of asymmetric competition
    %figure()
    subplot(2,2,2)
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
        set([p1,p2], {'LineWidth', 'MarkerSize', 'Marker'}, {5, 10, 'o'})
        p1.HandleVisibility = 'off';
        p1.MarkerFaceColor = 'white';
        p2.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    subplot(2, 2, 3)
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
    
    subplot(2, 2, 4)
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
    
    %% Evolved Trait Differences Figure
    figure()
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
        set([p1,p2], {'LineWidth', 'MarkerSize', 'Marker'}, {5, 10, 'o'})
        p1.HandleVisibility = 'off';
        p1.MarkerFaceColor = 'white';
        p2.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    %% Stabilizing and Fitness Difference Figures with approx and full
    figure()
    subplot(2, 2, 1)
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
    hold off
    
    subplot(2, 2, 3)
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(["Sympatric - Allopatric"; "$\kappa'_1 - \kappa'_2$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(ybar)
        p = plot(ydiff, SpecAveFit3(i,:) - SpecAveFit3_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    subplot(2, 2, 2)
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(["Sympatric - Allopatric"; "$-\Delta N + \Delta I_G + \Delta I_V$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(ybar)
        p = plot(ydiff, Aapprox3(i,:) - Aapprox3_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(2, 2, 4)
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(['Sympatric - Allopatric';...
        "$\kappa'_1 - \kappa'_2$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(ybar)
        p = plot(ydiff, kappadiffPrime3(i,:) - kappadiffPrime3_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    
    % Seperating into different figures
    figure()
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel("Sympatric - Allopatric $\bar{A}$", 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(ybar)
        p = plot(ydiff, A3(i,:) - A3_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    figure()
    plot(ydiff, zeros(1,length(ydiff)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1$ - ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel(["Sympatric - Allopatric";...
        "$\kappa'_1 - \kappa'_2$"], 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(ybar)
        p = plot(ydiff, SpecAveFit3(i,:) - SpecAveFit3_allo(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    %% Full stabilizing mechanisms and average fitness differences
    fig = figure();
    set(fig, 'defaultAxesColorOrder', zeros(2,3));
    
    yyaxis left
    plot([0, max(ydiff)], [-0.02, 0.1], 'Color', 'none', 'HandleVisibility', 'off')
    axis([0, max(ydiff), mincoex, maxcoex])
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('ln$y_1 - $ln$y_2$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('Sympatric - Allopatric $\bar{A}$', 'FontSize', 35, 'Interpreter', 'Latex')
    
    yyaxis right
    ylabel(["Sympatric - Allopatric"; "$|\kappa_1 - \kappa_2|$"], 'FontSize', 35, 'Interpreter', 'Latex')
    axis([0, max(ydiff), mincoex, maxcoex])
    hold on
    
    for i = 1:length(ybar)
        p1 = plot(ydiff, A3(i,:) - A3_allo(i,:), '-', 'Color', c(i+1,:));
        p2 = plot(ydiff, abs(SpecAveFit3(i, :)) - abs(SpecAveFit3_allo(i,:)), ':', 'Color', c(i+1,:), 'HandleVisibility', 'off');
        set([p1, p2], {'LineWidth', 'MarkerSize', 'Marker'},...
                        {5, 10, 'o'});
        p1.MarkerFaceColor = c(i+1,:);
        p2.MarkerFaceColor = 'none';
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    
    %% Figure showing endvalue of theta g for each species
    figure()
    surf(ydiff, ybar, thetag1Sim3)
    hold on
    surf(ydiff, ybar, thetag2Sim3)
    hold off
    xlabel('ln$y_1$ - ln$y_2$', 'Interpreter', 'Latex')
    ylabel('0.5(ln$y_1$ + ln$y_2$)', 'Interpreter', 'Latex')
    zlabel('Sympatric ESC $\theta_G$', 'Interpreter', 'Latex')
end


%%
if RUN(4) == 1
    %% NOW SECTION ON EFFECT OF DIFFERENCES IN Y
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
    kappaDiff4 = zeros(length(sigmaV), length(rho));
    ComAveDeltaN4 = zeros(length(sigmaV), length(rho));
    ComAveDeltaIG4 = zeros(length(sigmaV), length(rho));
    ComAveDeltaIV4 = zeros(length(sigmaV), length(rho));
    
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
            
            [A4(i, j), SpecAveFit4(i, j), kappaDiff4(i, j), ComAveDeltaN4(i, j), ComAveDeltaIG4(i, j), ComAveDeltaIV4(i, j)]...
                = EcologicalDynamics(thetagSim4, thetav, gen, burnin, Parameters);
            
            [i/length(sigmaV), j/length(rho)]
        end
    end
    
    Aapprox4 = -ComAveDeltaN4 + ComAveDeltaIG4 + ComAveDeltaIV4;
    rbarapprox4 = kappaDiff4 + Aapprox4;
    kappadiffPrime4 = rbarapprox4 - Aapprox4;
    
    r14 = SpecAveFit4 + A4;
    r24 = A4 - SpecAveFit4;
    coex4 = r14 > 0 & r24 > 0;
    deltagSimCoex4 = NaN(length(sigmaV), length(rho));
    deltagSimCoex4(coex4) = deltagSim4(coex4);
    mincoex = min(min([abs(SpecAveFit4), A4, Aapprox4, kappadiffPrime4]));
    maxcoex = max(max([abs(SpecAveFit4), A4, Aapprox4, kappadiffPrime4]));
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
    
    %% Stabilizing and Fitness Difference Figures with approx and full
    figure()
    subplot(2, 2, 1)
    plot(rho, zeros(1,length(rho)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('$\rho$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('$A$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(rho, A4(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(2, 2, 3)
    plot(rho, zeros(1,length(rho)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('$\rho$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel({'ESC Species Average Fitness';...
        'Differences, $\kappa_1 - \kappa_2$'}, 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(rho, SpecAveFit4(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    leg = legend();
    set(leg, {'Interpreter', 'Location', 'String', 'Box'},...
        {'Latex', 'east', legendstring, 'off'});
    
    subplot(2, 2, 2)
    plot(rho, zeros(1,length(rho)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('$\rho$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel('ESC Community Stabilization, $A$', 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(rho, Aapprox4(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
    subplot(2, 2, 4)
    plot(rho, zeros(1,length(rho)), '--', 'Color', 'black')
    h = gca;
    h.FontSize = 30;
    h.FontName = 'Times New Roman';
    xlabel('$\rho$', 'FontSize', 35, 'Interpreter', 'Latex')
    ylabel({'ESC Species Average Fitness';...
        'Differences, $\kappa_1 - \kappa_2$'}, 'FontSize', 35, 'Interpreter', 'Latex')
    hold on
    
    for i = 1:length(sigmaV)
        p = plot(rho, kappadiffPrime4(i,:), 'Color', c(i+1,:));
        p.LineWidth = 5;
        p.MarkerSize = 10;
        p.Marker = 'o';
        p.MarkerFaceColor = c(i+1,:);
    end
    hold off
    
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
    
    
    %% Figure showing endvalue of theta g for each species
    figure()
    surf(rho, sigmaV, thetag1Sim4)
    hold on
    surf(rho, sigmaV, thetag2Sim4)
    hold off
    xlabel('$\rho$', 'Interpreter', 'Latex')
    ylabel('0.5(ln$y_1$ + ln$y_2$)', 'Interpreter', 'Latex')
    zlabel('Evolved $\theta_G$', 'Interpreter', 'Latex')
end