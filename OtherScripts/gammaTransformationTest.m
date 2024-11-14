% %Get nSR data from dataT
% load("Results/dataT_highRes_Sep12th24.mat")
% 
% nSRcounts500 = dataT.nSRcounts500;
% highSRCoresLog = dataT.meanSR > 8;
% 
% %Set up bin edges on inverse nSR
% invXbinEdges = [0:0.02:15, 100];
% 
% %Fit a gamma distribution to inverse nSR counts, and return histogram
% %counts
% [invSRx, gammaSR, gammainvSR, alpha,invSRhist] = fitGamma2invSR(nSRcounts500, highSRCoresLog, invXbinEdges);

%% Load gamma distribution from python (alpha = 3.82)
shape = 3.82;

invnSR_gamma_py = readmatrix("Results/invnSR_gamma.csv");

x = invnSR_gamma_py(:,1);

%Fit gamma pdf to same x values
invnSR_gamma = gampdf(x, shape, 1/shape);

linewidth = 1;
gamFig = figure;
subplot(2,1,1)
plot(x, invnSR_gamma_py(:,2), 'k-', "DisplayName", "Python Gamma", "LineWidth", linewidth)
hold on
%plot(x, invnSR_gamma, 'r--', "DisplayName", "Matlab Gamma", "LineWidth", linewidth)
xlabel("inverse nSR")
xlim([0 5])
legend()

%% Load inverse gamma distribution from python (alpha = 3.82)

nSR_inversegamma_py = readmatrix("Results/nSR_inversegamma.csv");

%Convert gamma pdf to nSR using my own function

[nSRx, nSR_gammaTransformed] = invSRtoSR(x, invnSR_gamma_py(:,2));

figure(gamFig)
subplot(2,1,2)
plot(nSR_inversegamma_py(:,1), nSR_inversegamma_py(:,2), 'k-',"DisplayName", "Python Inverse Gamma Distribution", "LineWidth", linewidth)
hold on
plot(nSRx, nSR_gammaTransformed, 'r--', "DisplayName", "Matlab Transformation of Python Gamma", "LineWidth", linewidth)
xlabel("nSR")
xlim([0 5])
legend()
