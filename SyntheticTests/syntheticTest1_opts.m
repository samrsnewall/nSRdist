%% Synthetic Test for recovering nSR distribution
%Synthetic test to evaluate whether we can recover some "true" nSR
%distribution from a number of cores with radiocarbon dates. We want to
%know if there is any bias is introduced or if there are consistency
%problems

%% Set up useful settings
S.evenSpacing = 'time'; %Set up whether the spacing is even in time ('time') or depth ('depth')
S.weighting = 'depth'; %Set up whether the weighting uses time ('time') or depth ('depth')
S.pdfMeanIsOne = 0; %Set up whether the mean of the lognormal pdf should be exactly 1;

%% Synthetic Set Up
%Set up the "true" nSR distribution
t_lognSR_var = 0.308524; %variance (in log space)
if S.pdfMeanIsOne
    t_lognSR_mu = -0.5*t_lognSR_var; %mean - This formula sets the mean of the lognormal to 1 in linear space
else
    t_lognSR_mu = 0.159458; %mean (in log space) from previous analyses
    t_lognSR_mu = 0.2;
end
t_lognSR_sigma = sqrt(t_lognSR_var);
t_nSR_pd = makedist("Lognormal", 'mu', t_lognSR_mu, 'sigma', t_lognSR_sigma);
t_lnSR_pd = makedist("Normal", "mu", t_lognSR_mu, 'sigma', t_lognSR_sigma);

%How many times to repeat synthetic test
J = 1000;

%Overall settings for synthetic cores
M = 100;
N = 10;
aveSR = 12;
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
    t_nSRsamps{m} = random(t_nSR_pd, N, 1);

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
r_SRs = cell(M,1);
r_nSRs = cell(M,1);
for m = 1:M
    r_aveSR(m) = (t_coreDepths{m}(end)-t_coreDepths{m}(1))/(t_coreAges{m}(end) - t_coreAges{m}(1));
    r_SRs{m} = diff(t_coreDepths{m})./diff(t_coreAges{m});
    r_nSRs{m} = r_SRs{m}./r_aveSR(m);
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
nSRdist = fitdist(weighted_nSRs, 'Lognormal');
lnSRdist = fitdist(log(weighted_nSRs), "Normal");

r_mu(j) = nSRdist.mu;
r_sigma(j) = nSRdist.sigma;

if j == 1
%% Compare true with fitted
%Compare nSRs against true distribution

subplot(4,3,[4,5,6])
yyaxis left
histogram(log(allnSRs), 'DisplayName', "Raw nSR")

yyaxis right
xtoplot = -4:0.01:4;
y_trueDist = pdf(t_lnSR_pd, xtoplot);
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

subplot(4,3,[10,11,12])
y_recDist = pdf(lnSRdist, xtoplot);
plot(xtoplot, y_recDist, 'Color', 'r', "DisplayName", "Reconstructed Distribution")
hold on
plot(xtoplot, y_trueDist, 'Color', 'blue',"DisplayName", "True Distribution")
legend()
xlabel("log nSR")
end
end

%% Compare all reconstructed parameters with true parameters

figure;
plot(r_mu, r_sigma, 'kx', 'DisplayName', 'Reconstructed Parameters')
hold on
plot(t_lognSR_mu, t_lognSR_sigma, 'Marker','o','MarkerEdgeColor', 'k',  'MarkerFaceColor', 'b', 'DisplayName', 'True Parameters', 'LineStyle','none')
xlabel("Mean")
ylabel("Sigma")
legend()

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
