%% README
%Compare many bacon gamma acc rate distributions to BIGMACS nSR dist

addpath('Functions')

%% Load BIGMACS nSR distribution
BMLogNorm = load("lognormal_BIGMACS.txt");

%% Interp values onto same scale as gamma values will be plotted on
[convnSRExample, convnSRprobExample] = gammaAccRate2nSR(3/2);
BMLogNormProbInterp = interp1(BMLogNorm(:,1), BMLogNorm(:,2), convnSRExample);

%% Convert acc rate dists to nSR
% The value of alpha used is described in halves, i.e. alpha = 1.5 = 3/2,
% is called alpha3h

% [alph1h.nSR, alph1h.nSRprob] = gammaAccRate2nSR(1/2);
% [alph2h.nSR, alph2h.nSRprob] = gammaAccRate2nSR(2/2);
% [alph3h.nSR, alph3h.nSRprob] = gammaAccRate2nSR(3/2);
% [alph4h.nSR, alph4h.nSRprob] = gammaAccRate2nSR(4/2);
% [alph5h.nSR, alph5h.nSRprob] = gammaAccRate2nSR(5/2);
% [alph6h.nSR, alph6h.nSRprob] = gammaAccRate2nSR(6/2);
% [alph7h.nSR, alph7h.nSRprob] = gammaAccRate2nSR(7/2);

figure;
subplot(2,1,1)
hold on
plot(BMLogNorm(:,1), BMLogNorm(:,2), 'k-', 'DisplayName', "BIGMACS MLN", 'LineWidth', 2)

for i = [1,2,3,4,5,6,7]
    [convnSR, convnSRprob] = gammaAccRate2nSR(i/2);
    plot(convnSR, convnSRprob, '-', 'DisplayName', "alpha = " + num2str(i)+ "/2", 'LineWidth', 1)
end

xlabel("nSR")
ylabel("Probability Density")
xlim([0 6])
legend('location', 'northeast')

subplot(2,1,2)
hold on

for i = [1,2,3,4,5,6,7]
    [convnSR, convnSRprob] = gammaAccRate2nSR(i/2);
    plot(convnSR, BMLogNormProbInterp - convnSRprob, '-', 'DisplayName', "alpha = " + num2str(i)+ "/2", 'LineWidth', 1)
    % NEED TO FIX THE FACT THAT COLOURS ON SUBPLOT 1 and 2 are different!
end

xlabel("nSR")
ylabel("Probability Density Anomaly")
xlim([0 6])
legend('location', 'northeast')