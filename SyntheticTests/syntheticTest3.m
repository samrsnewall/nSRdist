%Synthetic test to evaluate whether we can recover some "true" nSR
%distribution from a number of cores with radiocarbon dates. We want to
%know if there is any bias is introduced or if there are consistency
%problems...

%% Synthetic Set Up
%Set up the "true" nSR distribution
t_lognSR_mu = 0.16; %mean (in log space)
t_lognSR_var = 0.31; %variance (in log space)
t_lognSR_sigma = sqrt(t_lognSR_var);
t_nSR_pd = makedist("Lognormal", 'mu', t_lognSR_mu, 'sigma', t_lognSR_sigma);
t_lnSR_pd = makedist("Normal", "mu", t_lognSR_mu, 'sigma', t_lognSR_sigma);

%How many times to repeat synthetic test
J = 1000;

%Create a set of M cores with N nSR estimates each
M = 100;
N = 10;
aveSR = 12;
t_depthSpacing = 15;

%Initialize vector to hold reconstructed values
r_mu = NaN(J,1);
r_sigma = NaN(J,1);

for j = 1:J


%Initialize some cells
t_nSRsamps = cell(M,1);
t_coreDepths = cell(M,1);
t_aveSR = NaN(M,1);
t_SRsamps = cell(M,1);
t_coreAgeIncrements = cell(M,1);
t_coreAges = cell(M,1);


for m = 1:M
    %Set up the depth spacing of each core
    t_coreDepths{m} = (0:t_depthSpacing:t_depthSpacing*(N))';

    %Sample some nSR values
    t_nSRsamps{m} = random(t_nSR_pd, N, 1);

    %Get the average SR for the core
    t_aveSR(m) = aveSR;
    t_SRsamps{m} = t_nSRsamps{m}.*t_aveSR(m);
    
    %Using nSR values and average SR, set up age values at each interval
    t_coreAgeIncrements{m} = diff(t_coreDepths{m})./t_SRsamps{m};
    t_coreAges{m} = [0; cumsum(t_coreAgeIncrements{m})];
end

if j == 1
%test construction of cores, plot to test
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
end

%% Perform Evaluation using our method
%Now, evaluate using methodology in our paper
allnSRs = [];
r_aveSR = NaN(M,1);
r_SRs = cell(M,1);
r_nSRs = cell(M,1);
for m = 1:M
    r_aveSR(m) = t_aveSR(m);
    r_SRs{m} = diff(t_coreDepths{m})./diff(t_coreAges{m});
    r_nSRs{m} = r_SRs{m}./r_aveSR(m);
    allnSRs = [allnSRs; r_nSRs{m}];
end

%Fit distribution
nSRdist = fitdist(allnSRs, 'Lognormal');
lnSRdist = fitdist(log(allnSRs), "Normal");

r_mu(j) = nSRdist.mu;
r_sigma(j) = nSRdist.sigma;

if j == 1
%% Compare true with fitted
%Compare nSRs against true distribution
figure;
subplot(2,1,1)
yyaxis left
histogram(log(allnSRs), 'DisplayName', "Reconstructed nSR")

yyaxis right
xtoplot = -4:0.01:4;
y_trueDist = pdf(t_lnSR_pd, xtoplot);
plot(xtoplot, y_trueDist, 'Color', 'blue', 'DisplayName', "True Distribution")
xlabel("log nSR")
legend()

subplot(2,1,2)
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






