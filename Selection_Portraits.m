clear
clc

% Determine Whether to run 1) effect of competition strength, 2) effect of
% competitive asymmetry, and 3) effect of differences in \theta_V
RUN = [0,0,1];

% Set to 1 if you want the initial condition for the simulation to be reps
% number of random initializations.
RAND = 0;
reps = 30;
shortmutnum = 100;
mutnum = 4000;

% Set to 1 if you want all selection vectors present. Set to 0 if you only
% want selection vectors for the regions of coexistence.
ALLSELVECT = 1;

coex_res = 200;
sel_res = 13;
gen = 50000; 
burnin = 1000;
epsilon = 0.05;

muG = 0;
muV = 0;
sigmaG = sqrt(0.5);
sigmaV = sqrt(0.5);
rho = 1;
alpha = 1;
s = 0.9;

Gbarstar = exp(muG)/(1+exp(muG));
Vbarstar = exp(muV);


if RUN(1) == 1
    %% Effect of strength of competition figures in polar coordinates
    
    diffv = 0*pi;
    d = 0;
    ymean = [3, 6, 9];
    ydiff = 0;
    
    figure()
    title_string = {['(a) Weak Competition, $y$ = ', num2str(ymean(1))],...
        ['(b) Intermediate Competition, $y$ = ', num2str(ymean(2))],...
        ['(c) Strong Copmetition, $y$ = ', num2str(ymean(3))]};
    
    thetag = linspace(-pi/2, pi/2, coex_res);
    
    sel_subset = round(linspace(1,coex_res-1, sel_res));
    thetag1sel = thetag(sel_subset);
    thetag2sel = thetag(sel_subset);
    
    thetav1 = diffv/2;
    thetav2 = -diffv/2;
    thetav = [thetav1, thetav2];
    
    deltagv1 = thetag1sel'*ones(1, sel_res) - thetav1*ones(sel_res);
    deltagv2 = ones(sel_res, 1)*thetag2sel - thetav2*ones(sel_res);
    
    deltag = zeros(sel_res);
    
    for k = 1:length(ymean)
        y = ymean(k) + 0.5*ydiff*[1, -1];
        
        Nbarstar = (y*Gbarstar*Vbarstar/(1 - s*(1-Gbarstar)) - 1)/(alpha*Gbarstar*Vbarstar);

        invr1 = zeros(coex_res);
        invr2 = zeros(coex_res);
        
        sel_grad1 = zeros(sel_res);
        sel_grad2 = zeros(sel_res);
        
        deltag_grad = zeros(sel_res);
        deltagv1_grad = zeros(sel_res);
        deltagv2_grad = zeros(sel_res);
        
        Xvec = mvnrnd(zeros(1,4), [eye(2),rho*eye(2);rho*eye(2),eye(2)], gen);
        X = Xvec(:,1:2);
        Y = Xvec(:,3:4);
        
        parfor i = 1:length(thetag)
            
            thetagvec = [thetag(i), thetag];
            EG = muG + sigmaG*([sin(thetagvec'), cos(thetagvec')]*X');
            EV = sigmaV*([sin(thetav'), cos(thetav')]*Y');
            
            G = exp(EG)./(1+exp(EG));
            Gres = G(1, :);
            Ginv = G(2:end, :);
            
            V = exp(EV);
            
            % Evaluate coexistence
            % Resident Dynamics
            N = zeros(2, gen);
            N(:,1) = Nbarstar;
            lambda = zeros(2, gen-1);
            C = zeros(2, gen);
            C(:, 1) = 1 + alpha*Gres(:, 1).*V(:, 1).*N(:, 1);
            
            for t = 2:gen
                lambda(:, t-1) = s*(1-Gres(:, t-1)) + y'.*Gres(:, t-1).*V(:, t-1)./C(:, t-1);
                N(:, t) = N(:, t-1).*lambda(:, t-1);
                C(:, t) = 1 + alpha*Gres(:, t).*V(:, t).*N(:, t);
            end
            
            % Calculate Invader Growth Rates
            invr1(:, i) = mean(log(s*(1-Ginv(:, burnin:end)) + y(1)*Ginv(:, burnin:end).*(ones(length(thetag),1).*V(1, burnin:end))./(ones(length(thetag), 1).*C(2, burnin:end))), 2);
            invr2(i, :) = mean(log(s*(1-Ginv(:, burnin:end)) + y(2)*Ginv(:, burnin:end).*(ones(length(thetag),1).*V(2, burnin:end))./(ones(length(thetag), 1).*C(1, burnin:end))), 2)';
            
            [k, i/length(thetag)]
        end
        
        coex = invr1 > 0 & invr2 > 0;
        A = diag(diag(coex));
        coex = coex - A;
        
        for i = 1:length(thetag1sel)
            deltag(i,:) = thetag1sel(i) - thetag2sel;
            
            for j = 1:length(thetag2sel)
                
                if coex(sel_subset(i), sel_subset(j)) == 1 || ALLSELVECT == 1
                    
                    thetagvec = [thetag1sel(i), thetag2sel(j)];
                    EG = muG + sigmaG*([sin(thetagvec'), cos(thetagvec')]*X');
                    EV = sigmaV*([sin(thetav'), cos(thetav')]*Y');
                    
                    G = exp(EG)./(1+exp(EG));
                    V = exp(EV);
                    
                    % Evaluate Selection
                    % Calculate Resident Dynamics with 2-species
                    N = zeros(2, gen);
                    N(:,1) = Nbarstar/2;
                    lambda = zeros(2, gen-1);
                    C = zeros(1, gen);
                    C(1) = 1 + alpha*(G(:, 1).*V(:, 1))'*N(:, 1);
                    
                    for t = 2:gen
                        lambda(:, t-1) = s*(1-G(:, t-1)) + y'.*G(:, t-1).*V(:, t-1)/C(t-1);
                        N(:, t) = N(:, t-1).*lambda(:, t-1);
                        C(t) = 1 + alpha*(G(:, t).*V(:, t))'*N(:, t);
                    end
                    
                    thetag1i = thetag1sel(i) + epsilon;
                    thetag2i = thetag2sel(j) + epsilon;
                    thetagi = [thetag1i, thetag2i];
                    
                    EGi = muG + sigmaG*([sin(thetagi'), cos(thetagi')]*X');
                    Gi = exp(EGi)./(1+exp(EGi));
                    
                    rbar_wt = mean(log(s*(1-G(:, burnin:end)) + y'*ones(1,length(burnin:gen)).*G(:, burnin:end).*V(:, burnin:end)./(ones(2,1)*C(burnin:end))), 2);
                    rbar_mut = mean(log(s*(1-Gi(:, burnin:end)) + y'*ones(1,length(burnin:gen)).*Gi(:, burnin:end).*V(:, burnin:end)./(ones(2,1)*C(burnin:end))), 2);
                    
                    sel_grad1(i,j) = (rbar_mut(1) - rbar_wt(1))/epsilon;
                    sel_grad2(i,j) = (rbar_mut(2) - rbar_wt(2))/epsilon;
                    [k, i/sel_res, j/sel_res]
                end
            end
        end
        
        k
        'Simulating Evolution'
        
        Parameters = {s,y,rho,muG,muV,sigmaG,sigmaV,alpha};
        
        if RAND == 0
            reps = 1;
            shortmutnum = mutnum;
        end
        
        evothetag1 = zeros(reps,shortmutnum+1);
        evothetag2 = zeros(reps,shortmutnum+1);
        
        for i = 1:reps
            if RAND == 1
                thetaginit = pi*(rand(1,2) - 0.5);
            else
                thetaginit = -pi/2*[1,-1]*0.2;
            end
            [evothetag1(i,:), evothetag2(i,:), MeanN] = AdaptDynamics22_PolarCoord(thetaginit, diffv/2*[1,-1], gen, shortmutnum, Parameters);
            i/reps
        end
        
        subplot(1, length(ymean), k)
        contourf(thetag, thetag, coex', 0.5*[1, 1])
        colormap(gray)
        hold on
        quiver(thetag1sel'*ones(1, sel_res), ones(sel_res, 1)*thetag2sel, sel_grad1, sel_grad2,...
            'color', 'black', 'LineWidth', 2, 'MaxHeadSize', 0.4)
        plot(thetav1*ones(1, 2), pi/2*[-1, 1], '--', 'color', 'black')
        plot(pi/2*[-1, 1], thetav2*ones(1, 2), '--', 'color', 'black')
        for i = 1:reps
            scatter(evothetag1(i,:), evothetag2(i,:), 'filled', 'CData', viridis(shortmutnum+1))
        end
        hold off
        axis(pi/2*[-1, 1, -1, 1])
        ax = gca;
        ax.XTick = [-pi/2, 0, pi/2];
        ax.YTick = [-pi/2, 0, pi/2];
        ax.XTickLabel = {'-\pi/2', 0, '\pi/2'};
        ax.YTickLabel = {'-\pi/2', 0, '\pi/2'};
        set(ax, {'FontName', 'FontSize'}, {'Times New Roman', 25});
        xlab = xlabel('Species Trait 1, $\theta_{G_1}$');
        ylab = ylabel('Species Trait 2, $\theta_{G_2}$');
        tit = title(title_string(k));
        set([xlab, ylab], {'Interpreter', 'FontSize'}, {'Latex', 30});
        tit.Interpreter = 'Latex';
        tit.FontSize = 25;
        
    end
    
end

if RUN(2) == 1
    %% Effect of competitive asymmetry polar coordinates
    
    diffv = 0*pi;
    d = 0;
    ymean = 7;
    ydiff = 0.2;%[0, 0.2];
    
    figure()
    title_string = 'Asymmetric Competition';
    %title_string = {['(a) Symmetric Competition, $y_1 - y_2$ = ', num2str(ydiff(1))],...
    %    ['(b) Asymmetric Competition, $y_1 - y_2$ = ', num2str(ydiff(2))]};
    
    thetag = linspace(-pi/2, pi/2, coex_res);
    
    sel_subset = round(linspace(1,coex_res-1, sel_res));
    thetag1sel = thetag(sel_subset);
    thetag2sel = thetag(sel_subset);
    
    thetav1 = diffv/2;
    thetav2 = -diffv/2;
    thetav = [thetav1, thetav2];
    
    deltagv1 = thetag1sel'*ones(1, sel_res) - thetav1*ones(sel_res);
    deltagv2 = ones(sel_res, 1)*thetag2sel - thetav2*ones(sel_res);
    
    deltag = zeros(sel_res);
    
    for k = 1:length(ydiff)
        y = ymean + ydiff(k)/2*[1, -1];
        
        Nbarstar = (y*Gbarstar*Vbarstar/(1 - s*(1-Gbarstar)) - 1)/(alpha*Gbarstar*Vbarstar);

        invr1 = zeros(coex_res);
        invr2 = zeros(coex_res);
        
        sel_grad1 = zeros(sel_res);
        sel_grad2 = zeros(sel_res);
        
        deltag_grad = zeros(sel_res);
        deltagv1_grad = zeros(sel_res);
        deltagv2_grad = zeros(sel_res);
        
        Xvec = mvnrnd(zeros(1,4), [eye(2),rho*eye(2);rho*eye(2),eye(2)], gen);
        X = Xvec(:,1:2);
        Y = Xvec(:,3:4);
        
        for i = 1:length(thetag)
            
            thetagvec = [thetag(i), thetag];
            EG = muG + sigmaG*([sin(thetagvec'), cos(thetagvec')]*X');
            EV = muV + sigmaV*([sin(thetav'), cos(thetav')]*Y');
            
            G = exp(EG)./(1+exp(EG));
            Gres = G(1, :);
            Ginv = G(2:end, :);
            
            V = exp(EV);
            
            % Evaluate coexistence
            % Resident Dynamics
            N = zeros(2, gen);
            N(:,1) = Nbarstar;
            lambda = zeros(2, gen-1);
            C = zeros(2, gen);
            C(:, 1) = 1 + alpha*Gres(:, 1).*V(:, 1).*N(:, 1);
            
            for t = 2:gen
                lambda(:, t-1) = s*(1-Gres(:, t-1)) + y'.*Gres(:, t-1).*V(:, t-1)./C(:, t-1);
                N(:, t) = N(:, t-1).*lambda(:, t-1);
                C(:, t) = 1 + alpha*Gres(:, t).*V(:, t).*N(:, t);
            end
            
            % Calculate Invader Growth Rates
            invr1(:, i) = mean(log(s*(1-Ginv(:, burnin:end)) + y(1)*Ginv(:, burnin:end).*(ones(length(thetag),1).*V(1, burnin:end))./(ones(length(thetag), 1).*C(2, burnin:end))), 2);
            invr2(i, :) = mean(log(s*(1-Ginv(:, burnin:end)) + y(2)*Ginv(:, burnin:end).*(ones(length(thetag),1).*V(2, burnin:end))./(ones(length(thetag), 1).*C(1, burnin:end))), 2)';
            
            [k, i/length(thetag)]
        end
        
        coex = invr1 > 0 & invr2 > 0;
        A = diag(diag(coex));
        coex = coex - A;
        
        for i = 1:length(thetag1sel)
            deltag(i,:) = thetag1sel(i) - thetag2sel;
            
            for j = 1:length(thetag2sel)
                
                if coex(sel_subset(i),sel_subset(j)) == 1 || ALLSELVECT == 1
                    
                    thetagvec = [thetag1sel(i), thetag2sel(j)];
                    EG = muG + sigmaG*([sin(thetagvec'), cos(thetagvec')]*X');
                    EV = muV + sigmaV*([sin(thetav'), cos(thetav')]*Y');
                    
                    G = exp(EG)./(1+exp(EG));
                    V = exp(EV);
                    
                    % Evaluate Selection
                    % Calculate Resident Dynamics with 2-species
                    N = zeros(2, gen);
                    N(:,1) = Nbarstar/2;
                    lambda = zeros(2, gen-1);
                    C = zeros(1, gen);
                    C(1) = 1 + alpha*(G(:, 1).*V(:, 1))'*N(:, 1);
                    
                    for t = 2:gen
                        lambda(:, t-1) = s*(1-G(:, t-1)) + y'.*G(:, t-1).*V(:, t-1)/C(t-1);
                        N(:, t) = N(:, t-1).*lambda(:, t-1);
                        C(t) = 1 + alpha*(G(:, t).*V(:, t))'*N(:, t);
                    end
                    
                    thetag1i = thetag1sel(i) + epsilon;
                    thetag2i = thetag2sel(j) + epsilon;
                    thetagi = [thetag1i, thetag2i];
                    
                    EGi = muG + sigmaG*([sin(thetagi'), cos(thetagi')]*X');
                    Gi = exp(EGi)./(1+exp(EGi));
                    
                    rbar_wt = mean(log(s*(1-G(:, burnin:end)) + y'*ones(1,length(burnin:gen)).*G(:, burnin:end).*V(:, burnin:end)./(ones(2,1)*C(burnin:end))), 2);
                    rbar_mut = mean(log(s*(1-Gi(:, burnin:end)) + y'*ones(1,length(burnin:gen)).*Gi(:, burnin:end).*V(:, burnin:end)./(ones(2,1)*C(burnin:end))), 2);
                    
                    sel_grad1(i, j) = (rbar_mut(1) - rbar_wt(1))/epsilon;
                    sel_grad2(i, j) = (rbar_mut(2) - rbar_wt(2))/epsilon;
                    [k, i/sel_res, j/sel_res]

                end
            end
        end
        
        k
        ['Simulating Evolution']
        
        Parameters = {s,y,rho,muG,muV,sigmaG,sigmaV,alpha};
        
        if RAND == 0
            reps = 1;
            shortmutnum = mutnum;
        end
        
        evothetag1 = zeros(reps,shortmutnum+1);
        evothetag2 = zeros(reps,shortmutnum+1);
        
        for i = 1:reps
            if RAND == 1
                thetaginit = pi*(rand(1,2) - 0.5);
            else
                if k == 1
                    thetaginit = 0.4*pi/2*[1, -1];
                else
                    thetaginit = [0.85, -0.25]*pi/2;
                end
            end
            [evothetag1(i,:), evothetag2(i,:), MeanN] = AdaptDynamics22_PolarCoord(thetaginit, diffv/2*[1,-1], gen, shortmutnum, Parameters);
            i/reps
        end
        
        subplot(1, length(ydiff), k)
        contourf(thetag, thetag, coex', 0.5*[1, 1])
        colormap(gray)
        hold on
        quiver(thetag1sel'*ones(1, sel_res), ones(sel_res, 1)*thetag2sel, sel_grad1, sel_grad2,...
            'color', 'black', 'LineWidth', 2, 'MaxHeadSize', 0.4)
        plot(thetav1*ones(1, 2), pi/2*[-1, 1], '--', 'color', 'black', 'LineWidth', 2)
        plot(pi/2*[-1, 1], thetav2*ones(1, 2), '--', 'color', 'black', 'LineWidth', 2)
        for i = 1:reps
            scatter(evothetag1(i,:), evothetag2(i,:), 'filled', 'CData', viridis(shortmutnum+1))
        end
        hold off
        axis(pi/2*[-1, 1, -1, 1])
        ax = gca;
        ax.XTick = [-pi/2, 0, pi/2];
        ax.YTick = [-pi/2, 0, pi/2];
        ax.XTickLabel = {'-\pi/2', 0, '\pi/2'};
        ax.YTickLabel = {'-\pi/2', 0, '\pi/2'};
        set(ax, {'FontName', 'FontSize'}, {'Times New Roman', 25});
        xlab = xlabel('Species Trait 1, $\theta_{G_1}$');
        ylab = ylabel('Species Trait 2, $\theta_{G_2}$');
        tit = title(title_string(k));
        set([xlab, ylab], {'Interpreter', 'FontSize'}, {'Latex', 30});
        tit.Interpreter = 'Latex';
        tit.FontSize = 25;
        
    end
end

if RUN(3) == 1
    %% Effect of Displacement in polar coordinates
    
    diffv = [1/4]*pi;
    d = 0;
    ymean = 4;
    ydiff = 0;
    y = ymean + ydiff/2*[1, -1];
    Nbarstar = (y*Gbarstar*Vbarstar/(1 - s*(1-Gbarstar)) - 1)/(alpha*Gbarstar*Vbarstar);

    
    figure()
    subplot(3,1,1)
    title_string = {['$\theta_{V_1} - \theta_{V_2}$ = ', num2str(diffv(1))]};
    
    thetag = linspace(-pi/2, pi/2, coex_res);
    
    sel_subset = round(linspace(1,coex_res-1, sel_res));
    thetag1sel = thetag(sel_subset);
    thetag2sel = thetag(sel_subset);
    
    deltag = zeros(sel_res);
    
    for k = 1:length(diffv)
        thetav1 = diffv(k)/2;
        thetav2 = -diffv(k)/2;
        thetav = [thetav1, thetav2];
        
        invr1 = zeros(coex_res);
        invr2 = zeros(coex_res);
        
        sel_grad1 = zeros(sel_res);
        sel_grad2 = zeros(sel_res);
        
        deltag_grad = zeros(sel_res);
        deltagv1_grad = zeros(sel_res);
        deltagv2_grad = zeros(sel_res);
        
        Xvec = mvnrnd(zeros(1,4), [eye(2),rho*eye(2);rho*eye(2),eye(2)], gen);
        X = Xvec(:,1:2);
        Y = Xvec(:,3:4);
        
        for i = 1:length(thetag)
            
            thetagvec = [thetag(i), thetag];
            EG = muG + sigmaG*([sin(thetagvec'), cos(thetagvec')]*X');
            EV = muV + sigmaV*([sin(thetav'), cos(thetav')]*Y');
            
            G = exp(EG)./(1+exp(EG));
            Gres = G(1, :);
            Ginv = G(2:end, :);
            
            V = exp(EV);
            
            % Evaluate coexistence
            % Resident Dynamics
            N = zeros(2, gen);
            N(:,1) = Nbarstar;
            lambda = zeros(2, gen-1);
            C = zeros(2, gen);
            C(:, 1) = 1 + alpha*Gres(:, 1).*V(:, 1).*N(:, 1);
            
            for t = 2:gen
                lambda(:, t-1) = s*(1-Gres(:, t-1)) + y'.*Gres(:, t-1).*V(:, t-1)./C(:, t-1);
                N(:, t) = N(:, t-1).*lambda(:, t-1);
                C(:, t) = 1 + alpha*Gres(:, t).*V(:, t).*N(:, t);
            end
            
            % Calculate Invader Growth Rates
            invr1(:, i) = mean(log(s*(1-Ginv(:, burnin:end)) + y(1)*Ginv(:, burnin:end).*(ones(length(thetag),1).*V(1, burnin:end))./(ones(length(thetag), 1).*C(2, burnin:end))), 2);
            invr2(i, :) = mean(log(s*(1-Ginv(:, burnin:end)) + y(2)*Ginv(:, burnin:end).*(ones(length(thetag),1).*V(2, burnin:end))./(ones(length(thetag), 1).*C(1, burnin:end))), 2)';
            
            [k, i/length(thetag)]
        end
        
        coex = invr1 > 0 & invr2 > 0;
        if diffv == 0
            A = diag(diag(coex));
            coex = coex - A;
        end
        
        for i = 1:length(thetag1sel)
            deltag(i,:) = thetag1sel(i) - thetag2sel;
            
            for j = 1:length(thetag2sel)
                
                if coex(sel_subset(i), sel_subset(j)) == 1 || ALLSELVECT == 1
                    
                    thetagvec = [thetag1sel(i), thetag2sel(j)];
                    EG = muG + sigmaG*([sin(thetagvec'), cos(thetagvec')]*X');
                    EV = muV + sigmaV*([sin(thetav'), cos(thetav')]*Y');
                    
                    G = exp(EG)./(1+exp(EG));
                    V = exp(EV);
                    
                    % Evaluate Selection
                    % Calculate Resident Dynamics with 2-species
                    N = zeros(2, gen);
                    N(:,1) = Nbarstar/2;
                    lambda = zeros(2, gen-1);
                    C = zeros(1, gen);
                    C(1) = 1 + alpha*(G(:, 1).*V(:, 1))'*N(:, 1);
                    
                    for t = 2:gen
                        lambda(:, t-1) = s*(1-G(:, t-1)) + y'.*G(:, t-1).*V(:, t-1)/C(t-1);
                        N(:, t) = N(:, t-1).*lambda(:, t-1);
                        C(t) = 1 + alpha*(G(:, t).*V(:, t))'*N(:, t);
                    end
                    
                    res_growth = mean(log(s*(1-G(:, burnin:end)) + (y'*ones(1, length(burnin:gen))).*G(:, burnin:end).*V(:, burnin:end)./(ones(2, 1)*C(burnin:end))), 2);
                    thetag1i = thetag1sel(i) + epsilon;
                    thetag2i = thetag2sel(j) + epsilon;
                    thetagi = [thetag1i, thetag2i];
                    
                    EGi = muG + sigmaG*([sin(thetagi'), cos(thetagi')]*X');
                    Gi = exp(EGi)./(1+exp(EGi));
 
                    rbar_wt = mean(log(s*(1-G(:, burnin:end)) + y'*ones(1,length(burnin:gen)).*G(:, burnin:end).*V(:, burnin:end)./(ones(2,1)*C(burnin:end))), 2);
                    rbar_mut = mean(log(s*(1-Gi(:, burnin:end)) + y'*ones(1,length(burnin:gen)).*Gi(:, burnin:end).*V(:, burnin:end)./(ones(2,1)*C(burnin:end))), 2);
                    
                    sel_grad1(i,j) = (rbar_mut(1) - rbar_wt(1))/epsilon;
                    sel_grad2(i,j) = (rbar_mut(2) - rbar_wt(2))/epsilon;
                    
                    [k, i/sel_res, j/sel_res]
                end
            end
        end
        
        k
        ['Simulating Evolution']
        
        Parameters = {s,y,rho,muG,muV,sigmaG,sigmaV,alpha};
        
        if RAND == 0
            reps = 1;
            shortmutnum = mutnum;
        end
        
        evothetag1 = zeros(reps,shortmutnum+1);
        evothetag2 = zeros(reps,shortmutnum+1);
        
        for i = 1:reps
            if RAND == 1
                thetaginit = pi*(rand(1,2) - 0.5);
            else
                thetaginit = -pi/2*[1,-1]*0.2;
            end
            [evothetag1(i,:), evothetag2(i,:), MeanN] = AdaptDynamics22_PolarCoord(thetaginit, diffv(k)/2*[1,-1], gen, shortmutnum, Parameters);
            i/reps
        end
        
        subplot(1, length(diffv), k)
        contourf(thetag, thetag, coex', 0.5*[1, 1])
        colormap(gray)
        hold on
        quiver(thetag1sel'*ones(1, sel_res), ones(sel_res, 1)*thetag2sel, sel_grad1, sel_grad2,...
            'color', 'black', 'LineWidth', 2, 'MaxHeadSize', 0.4)
        plot(thetav1*ones(1, 2), pi/2*[-1, 1], '--', 'color', 'black', 'LineWidth', 2)
        plot(pi/2*[-1, 1], thetav2*ones(1, 2), '--', 'color', 'black', 'LineWidth', 2)
        for i = 1:reps
            scatter(evothetag1(i,:), evothetag2(i,:), 'filled', 'CData', viridis(shortmutnum+1))
        end
        hold off
        axis(pi/2*[-1, 1, -1, 1])
        ax = gca;
        ax.XTick = [-pi/2, 0, pi/2];
        ax.YTick = [-pi/2, 0, pi/2];
        ax.XTickLabel = {'-\pi/2', 0, '\pi/2'};
        ax.YTickLabel = {'-\pi/2', 0, '\pi/2'};
        set(ax, {'FontName', 'FontSize'}, {'Times New Roman', 25});
        xlab = xlabel('Species Trait 1, $\theta_{G_1}$');
        ylab = ylabel('Species Trait 2, $\theta_{G_2}$');
        tit = title(title_string(k));
        set([xlab, ylab], {'Interpreter', 'FontSize'}, {'Latex', 30});
        tit.Interpreter = 'Latex';
        tit.FontSize = 25;
        
    end
end
