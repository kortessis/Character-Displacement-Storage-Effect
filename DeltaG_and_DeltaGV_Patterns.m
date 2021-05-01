d = 0.7*pi;

vdiff = linspace(-(1-d), (1-d), 10);
gdiff = linspace(-pi, pi, 10);


%Assume symmetry, which is to say that vdiff is centered around the
%midpoint of phenotype space, which is pi/2 because phenotypes are in the
%range [0,pi].

corrEG = zeros(length(gdiff), length(vdiff));
corrEV = zeros(length(gdiff), length(vdiff));
corrEGEV1 = zeros(length(gdiff), length(vdiff));
corrEGEV2 = zeros(length(gdiff), length(vdiff));

for i = 1:length(gdiff)
    
    g1 = (pi + gdiff(i))/2;
    g2 = (pi - gdiff(i))/2;
    
    for j = 1:length(vdiff)
        v1 = (pi + vdiff(j))/2 + d;
        v2 = (pi - vdiff(j))/2 + d;
        
        phenotypes = [g1,g2,v1,v2];
        diffmatrix = phenotypes'*ones(1,4) - ones(4,1)*phenotypes;
        corrmatrix = cos(diffmatrix);
        
        corrEG(i,j) = corrmatrix(1,2);
        corrEV(i,j) = corrmatrix(3,4);
        corrEGEV1(i,j) = corrmatrix(1,3);
        corrEGEV2(i,j) = corrmatrix(2,4);
    end
    i
end
cmap = magma;
subplot(2,2,1)
surf(gdiff, vdiff, corrEG')
x1 = xlabel('$g_1 - g_2$');
y1 = ylabel('$v_1 - v_2$');
colormap(magma());

subplot(2,2,2)
surf(gdiff, vdiff, corrEV')
x2 = xlabel('$g_1 - g_2$');
y2 = ylabel('$v_1 - v_2$');
colormap(magma());

subplot(2,2,3)
surf(gdiff, vdiff, corrEGEV1')
x3 = xlabel('$g_1 - g_2$');
y3 = ylabel('$v_1 - v_2$');
colormap(magma());

subplot(2,2,4)
surf(gdiff, vdiff, corrEGEV2')
x4 = xlabel('$g_1 - g_2$');
y4 = ylabel('$v_1 - v_2$');
colormap(magma());

set([x1,x2,x3,x4,y1,y2,y3,y4], {'Interpreter', 'FontSize'}, {'Latex',30})






mutnum = 2000;
[deltag,deltagv1,deltagv2,MeanN] = AdaptDynamics22_PolarCoord(0.5*pi, 0*pi, 0, 10000, mutnum, 0.9, 1, 5*ones(1,2), 1, 2);


figure()
subplot(2,2,1)
plot(0:mutnum, abs(deltag)/pi)
xlabel('Time (# Mutations)')
ylabel('(g_1 - g_2)/\pi')


subplot(2,2,2)
plot(abs(deltag)/pi, 1/(2*pi)*(abs(deltagv1) + abs(deltagv2)))
xlabel('(g_1 - g_2)/\pi')
ylabel('Species Average (g - v)/pi')


subplot(2,2,3)
plot(0:mutnum, abs(deltagv1)/pi)
xlabel('Time (# Mutations)')
ylabel('(g_1 - v_1)/\pi')


subplot(2,2,4)
plot(0:mutnum, abs(deltagv2)/pi)
xlabel('Time (# Mutations)')
ylabel('(g_2 - v_2)/\pi')


figure()
subplot(2,2,1)
plot(0:mutnum, cos(abs(deltag)))
x1 = xlabel('Time (Number of Mutations)');
y1 = ylabel('corr$(E_{G_1}, E_{G_2})$');


cmap = viridis(length(deltagv1));
subplot(2,2,2)
scatter(cos(abs(deltag)), (1/2)*(cos(abs(deltagv1)) + cos(abs(deltagv2))), 200, cmap, 'filled')
x2 = xlabel('corr$(E_{G_1}, E_{G_2})$');
y2 = ylabel('Species Average corr$(E_G, E_V)$');


subplot(2,2,3)
plot(0:mutnum, cos(abs(deltagv1)))
x3 = xlabel('Time (Number of Mutations)');
y3 = ylabel('corr$(E_{G_1}, E_{V_1})$');


subplot(2,2,4)
plot(0:mutnum, cos(abs(deltagv2)))
x4 = xlabel('Time (Number of Mutations)');
y4 = ylabel('corr$(E_{G_1}, E_{V_2})$');

set([x1,x2,x3,x4,y1,y2,y3,y4], {'Interpreter', 'FontSize'}, {'Latex', 30})

