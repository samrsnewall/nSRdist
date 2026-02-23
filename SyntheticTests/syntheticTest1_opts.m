%% Synthetic Test for recovering nSR distribution
%Synthetic test to evaluate whether we can recover some "true" nSR
%distribution from a number of cores with radiocarbon dates. We want to
%know if there is any bias is introduced or if there are consistency
%problems

%% Set up useful settings
S.evenSpacing   = 'depth';    %Set up whether the spacing is even in time ('time') or depth ('depth')
S.weighting     = 'depth';    %Set up whether the weighting uses time ('time') or depth ('depth')
S.X_SR_or_DT    = 'DT';       %Use a distribution of sedimentation rate (SR) or deposition time (DT; inverse of SR);
S.pdfMeanIsOne  = 1;          %Set up whether the mean of the lognormal pdf should be exactly 1;
S.useTrueSR     = 0;          %Use the "true" SR (1) when normalizing or the "observed" SR (0)
S.number_of_cores = 100;      %Number of cores to simulate
S.dates_per_core = 10;       %Number of dates per core
S.true_lognormX_variance = 0.3;% Variance 
S.fit_pdf_mean_is_one = true;

%% Synthetic Set Up
%Set up the "true" nSR distribution
t_lognX_var = S.true_lognormX_variance; %variance (in log space)
if S.pdfMeanIsOne
    t_lognX_mu = -0.5*t_lognX_var; %mean - This formula sets the mean of the lognormal to 1 in linear space
else
    t_lognX_mu = 0.159458;
end
t_lognX_sigma = sqrt(t_lognX_var);
t_nX_pd     = makedist("Lognormal", 'mu', t_lognX_mu, 'sigma', t_lognX_sigma);
t_lnX_pd    = makedist("Normal", "mu", t_lognX_mu, 'sigma', t_lognX_sigma);


%How many times to repeat synthetic test
J = 1000;

%Overall settings for synthetic cores
M = S.number_of_cores;
N = S.dates_per_core;
aveSR = 20;
if contains(S.evenSpacing,'depth')
    t_depthSpacing = 15;
elseif contains(S.evenSpacing, 'time')
    t_timeSpacing = 10/12;
else
    error("Spacing not defined")
end

%Initialize vector to hold reconstructed values
r_mu = NaN(J,1);
r_sigma = NaN(J,1);

% Perform simulations
for j = 1:J

%Initialize some cells
t_nSRsamps = cell(M,1);
t_coreDepths = cell(M,1);
t_aveSR = NaN(M,1);
t_SRsamps = cell(M,1);
t_coreAgeIncrements = cell(M,1);
t_coreAges = cell(M,1);
t_coreDepthIncrements = cell(M,1);

%Set up the cores
for m = 1:M
    %Sample some nSR values
    if contains(S.X_SR_or_DT, 'SR')
        t_nSRsamps{m} = random(t_nX_pd, N, 1);
    elseif contains(S.X_SR_or_DT, 'DT') %If the distribution is a distribution of DT, convert it to SR
        t_nSRsamps{m} = 1./random(t_nX_pd, N, 1);
    end

    %Choose an average SR for the core
    t_aveSR(m) = aveSR;

    %Set up the SR values
    t_SRsamps{m} = t_nSRsamps{m}.*t_aveSR(m);

    if contains(S.evenSpacing,'depth')
        %Set up the depth spacing of each core
        t_coreDepths{m} = (0:t_depthSpacing:t_depthSpacing*(N))';

        %Using nSR values and average SR, set up age values at each interval
        t_coreAgeIncrements{m} = diff(t_coreDepths{m})./t_SRsamps{m};
        t_coreAges{m} = [0; cumsum(t_coreAgeIncrements{m})];
    elseif contains(S.evenSpacing, 'time')
        t_coreAges{m} = (0:t_timeSpacing:t_timeSpacing*(N))';

        t_coreDepthIncrements{m} = diff(t_coreAges{m}).*t_SRsamps{m};
        t_coreDepths{m} = [0; cumsum(t_coreDepthIncrements{m})];
    end
end

if j == 1
%Plot some synthetic cores
figure;
for m = 1:9
    subplot(3,3,m)
    plot(t_coreDepths{m}, t_coreAges{m}, 'k*')
    hold on
    plot(t_coreDepths{m}, t_coreDepths{m}./t_aveSR(m), '--k')
    if m ==8
        xlabel('Depth (cm)')
    elseif m == 4
        ylabel("Age (kyr)")
    end
end
sgtitle("Some Synthetic Cores")
end

%% Perform Evaluation using our method
%Now, evaluate using methodology in our paper
allnSRs = [];
weights_age = [];
weights_depth = [];
r_aveSR = NaN(M,1);
r_meanSR = NaN(M,1);
r_SRs = cell(M,1);
r_nSRs = cell(M,1);
for m = 1:M
    r_aveSR(m) = (t_coreDepths{m}(end)-t_coreDepths{m}(1))/(t_coreAges{m}(end) - t_coreAges{m}(1));
    r_SRs{m} = diff(t_coreDepths{m})./diff(t_coreAges{m});
    r_meanSR(m) = mean(r_SRs{m});
    if S.useTrueSR == 1
        r_nSRs{m} = r_SRs{m}./aveSR;
    else
        %r_nSRs{m} = r_SRs{m}./r_aveSR(m);
        r_nSRs{m} = r_SRs{m}./r_meanSR(m);
    end
    allnSRs = [allnSRs; r_nSRs{m}];
    weights_age = [weights_age; diff(t_coreAges{m})];
    weights_depth = [weights_depth; diff(t_coreDepths{m})];
end

%Test whether the r_aveSR are properly representing the true aveSR
if j ==1
    figure
    subplot(4,3,[1,2,3])
    histogram(r_aveSR, 50)
    hold on
    xline(aveSR, '--r', 'LineWidth', 1)
    xlabel("Average SR (cm/kyr)")
    ylims = ylim;
    xlims = xlim;
    text(0.8*max(xlims), 0.8*max(ylims),"median = " + num2str(median(r_aveSR)) + ",\newline mean = " + num2str(mean(r_aveSR)))
    title("Whole Core SR of Synthetic Cores")
end

%Perform weighting
if contains(S.weighting, 'depth')
    weights = weights_depth;
elseif contains(S.weighting, 'time')
    weights = weights_age;
end
weighted_nSRs = makeWeightedReplicates(allnSRs, weights, 3, 1); %Weight by replicating data according to weighting

%Fit distribution
if S.fit_pdf_mean_is_one
    [mu_hat, sigma_hat, ~] = fitLognormalMean1(weighted_nSRs);

    r_mu(j) = mu_hat;
    r_sigma(j) = sigma_hat;
else
nSRdist = fitdist(weighted_nSRs, 'Lognormal');
lnSRdist = fitdist(log(weighted_nSRs), "Normal");

r_mu(j) = nSRdist.mu;
r_sigma(j) = nSRdist.sigma;
end

if j == 1
%% Compare true with fitted
%Compare nSRs against true distribution

subplot(4,3,[4,5,6])
yyaxis left
histogram(log(allnSRs), 'DisplayName', "Raw nSR")

yyaxis right
xtoplot = -4:0.01:4;
y_trueDist = pdf(t_lnX_pd, xtoplot);
plot(xtoplot, y_trueDist, 'Color', 'blue', 'DisplayName', "True Distribution")
xlabel("log nSR")
legend()

subplot(4,3,[7,8,9])
yyaxis left
histogram(log(weighted_nSRs), 'DisplayName', 'Weighted nSRs')

yyaxis right
plot(xtoplot, y_trueDist, 'Color', 'blue', 'DisplayName', "True Distribution")
xlabel("log nSR")
legend()

% subplot(4,3,[10,11,12])
% %y_recDist = pdf(lnSRdist, xtoplot);
% plot(xtoplot, y_recDist, 'Color', 'r', "DisplayName", "Reconstructed Distribution")
% hold on
% plot(xtoplot, y_trueDist, 'Color', 'blue',"DisplayName", "True Distribution")
% legend()
% xlabel("log nSR")
end
end

%% Compare all reconstructed parameters with true parameters

%Calculate bias of each parameter
bias_sigma = mean(r_sigma)-t_lognX_sigma;
bias_mu = mean(r_mu) - t_lognX_mu;

%Calculate variance of each parameter
var_sigma = var(r_sigma);
var_mu = var(r_mu);

%Calculate standard error of each parameter
mse_sigma = bias_sigma^2 + var_sigma;
mse_mu = bias_mu^2 + var_mu;

figure;
plot(r_mu, r_sigma, 'kx', 'DisplayName', 'Reconstructed Parameters')
hold on
plot(mean(r_mu), mean(r_sigma), 'Marker','o','MarkerEdgeColor', 'k',  'MarkerFaceColor', 'g', 'DisplayName', 'Expected Value (reconstructed)', 'LineStyle','none')
plot(t_lognX_mu, t_lognX_sigma, 'Marker','o','MarkerEdgeColor', 'k',  'MarkerFaceColor', 'b', 'DisplayName', 'True Parameters', 'LineStyle','none')
xlabel("Mean")
ylabel("Sigma")
formatSpec = '%.2g';
title({ ...
    "Cores = " + num2str(M) + "; AgesPerCore = " + num2str(N), ...
    "Bias in mu = " + num2str(bias_mu, formatSpec) + "; bias in sigma = " + num2str(bias_sigma, formatSpec) ...
});
legend()

%Calculate mean of pdf

mean_estimated_pdf = exp(r_mu+(r_sigma.^2)/2);

figure;
plot(mean_estimated_pdf)
title("Mean of estimated pdfs")
xlabel("Estimated pdf #")
ylabel("Mean of pdf")


%% Specify outcomes and write results
outputs.bias_mu = bias_mu;
outputs.bias_sigma = bias_sigma;
outputs.var_mu = var_mu;
outputs.var_sigma = var_sigma;
outputs.mse_mu = mse_mu;
outputs.mse_sigma = mse_sigma;
outputs.timestamp = datetime('now');


log_results(S, outputs, 'simulation_results_1.csv');


%%

function log_results(S, outputs, filename)

    if nargin < 3
        filename = 'simulation_results.csv';   % default
    end

    % Convert S struct to a table row
    S_table = struct2table(S, 'AsArray', true);

    % Convert outputs struct to table row
    outputs_table = struct2table(outputs, 'AsArray', true);

    % Combine them horizontally
    new_row = [S_table outputs_table];

    % If file does not exist, create it
    if ~isfile(filename)
        writetable(new_row, filename);
        return;
    end

    % Otherwise append without writing header
    writetable(new_row, filename, 'WriteMode', 'append', 'WriteVariableNames', false);
end


%% Functions
function[repData] = makeWeightedReplicates(data, weight, dataRoundingDP, weightingInflationMultiplier)
%%% This function takes some data and some weighting and creates
%%% the variable repData, which is a single variable that approximates the
%%% weighting of the data by replicating each x by approximately it's
%%% weighting. For example, if x=3 has a weighting y=4, then the value 3
%%% will be repeated 4 times in repData.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 data           = round(data,dataRoundingDP);                   %The data are rounded so that there are less unique values to help count weightings
 data           = data(data ~=0);                               %Zeros are removed (choice for use in gamma fitting);
 weight         = weight(data ~=0);
[data_u,~, IC]  = unique(data);                                 %The unique values of the data are found, with their indices
 weight_u       = accumarray(IC,weight);                        %The weighting of each unique value is combined
 weight_uR      = round(weight_u.*weightingInflationMultiplier);%The weighting of each unique value is converted to an integer
 repData        = repelem(data_u,weight_uR);                    %The unique values are replicated according to their integer weighting
end

function [mu_hat, sigma_hat, negloglik] = fitLognormalMean1(x)
%FITLOGNORMALMEAN1  MLE for a lognormal distribution with mean fixed to 1.
%
%   [mu_hat, sigma_hat, negloglik] = fitLognormalMean1(x)
%
%   Assumes X ~ Lognormal(mu, sigma^2) with constraint:
%       E[X] = exp(mu + 0.5 * sigma^2) = 1
%   =>  mu = -0.5 * sigma^2
%
%   Inputs
%   ------
%   x : vector of positive data
%
%   Outputs
%   -------
%   mu_hat      : MLE of mu under the mean=1 constraint
%   sigma_hat   : MLE of sigma under the mean=1 constraint
%   negloglik   : value of the negative log-likelihood at the optimum

    % Ensure column vector
    x = x(:);

    if any(x <= 0)
        error('All data must be positive for a lognormal model.');
    end

    % Negative log-likelihood as a function of sigma only (mu is implied)
    negloglik_fun = @(sigma) negloglik_lognormal_mean1(sigma, x);

    % Choose a reasonable search interval for sigma
    % (you can tweak these bounds if needed)
    sigma_lower = 1e-6;
    % Use variability in log(x) to set an upper bound that is not absurd
    sigma_upper = max(10 * std(log(x)), 1e-2);  

    % Constrain sigma > 0 via fminbnd over [sigma_lower, sigma_upper]
    [sigma_hat, negloglik] = fminbnd(negloglik_fun, sigma_lower, sigma_upper);

    % Recover mu from the constraint E[X] = 1
    mu_hat = -0.5 * sigma_hat.^2;
end


function nll = negloglik_lognormal_mean1(sigma, x)
%NEGLOGLIK_LOGNORMAL_MEAN1  Negative log-likelihood for lognormal with mean=1.
%
%   sigma > 0. The mean constraint implies mu = -0.5 * sigma^2.

    if sigma <= 0
        nll = Inf;
        return
    end

    n = numel(x);
    mu = -0.5 * sigma.^2;

    % Lognormal pdf:
    % f(x|mu,sigma) = (1 ./ (x * sigma * sqrt(2*pi))) ...
    %                 .* exp( -0.5 * ((log(x)-mu)./sigma).^2 )
    %
    % So the negative log-likelihood is:
    % nll = sum( log(sigma) + log(x) + 0.5*log(2*pi) + 0.5*((log(x)-mu)/sigma).^2 )

    z = (log(x) - mu) ./ sigma;

    nll = n * log(sigma) ...
        + sum(log(x)) ...
        + 0.5 * n * log(2*pi) ...
        + 0.5 * sum(z.^2);
end

