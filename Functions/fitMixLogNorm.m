function[mixLogNormPDF, mixNormPDF, gmfit] = fitMixLogNorm(data_linear, numComponents)
%%% Create a mixed log normal to fit some dataset

%% Set the range of x values of interest 
% (use x values currently used in BIGMACS)
lognorm_BIGMACS = readtable("lognormal_BIGMACS.txt");
x = lognorm_BIGMACS.Var1';
lx = log(x);

%% Run gmfit
%Set the number of iterations (default=500 often doesn't converge)
options = statset('MaxIter', 2000);

%Use GMFit to fit a Gaussian Mixture distribution to the log of the data
if ~iscolumn(data_linear)
    data_linear = data_linear';
end
data_log = log(data_linear);
gmfit = fitgmdist(data_log, numComponents, "Options", options);

%% Construct Mixed Log Normal and Mixed Normal
%Find the means of fit
mu_fit = gmfit.mu;

%Initialise std deviation vector and prob density vector
stddev_fit = NaN(1,numComponents);
y = zeros(length(x), numComponents);
ly = zeros(length(lx), numComponents);

%Pull out sttdev of each component gaussian, and then find the lognormal
%pdf using the mean and stddev of the component gaussians as the input
%parameters.
for i = 1:numComponents
    stddev_fit(i) = sqrt(gmfit.Sigma(:,:,i)); %taking the square converts the variance to a stddev
    y(:,i) = lognpdf(x, mu_fit(i), stddev_fit(i)); % Find the y values of a lognormal using the mean and sttdeviations taken from the mixed gaussian
    ly(:,i) = normpdf(lx, mu_fit(i), stddev_fit(i)); % Find the y values of a Normal using...
end

%Find the mixing component of each individual Gaussian
mix_fit = gmfit.ComponentProportion;

%Find the mixed log normal, summing the individual gaussians after
%multiplying by their mixing component
z = sum(y.*mix_fit, 2);
lz = sum(ly.*mix_fit, 2);
mixLogNormPDF = [x', z];
mixNormPDF = [lx', lz];

%% Find the mean and variance of the distributions
%meanAnalytical = mix_fit(1)*exp(mu_fit(1) + (stddev_fit(1).^2)/2) + mix_fit(2) * exp(mu_fit(2) + (stddev_fit(2).^2)/2)

%create samples from each gaussian dist
numsamps = 10000;
samps1 = normrnd(mu_fit(1), stddev_fit(1), round(mix_fit(1)*numsamps), 1);
samps2 = normrnd(mu_fit(2), stddev_fit(2), round(mix_fit(2)*numsamps), 1);

%Calculate the mean of the Mixed Gaussian Distribution in LogSRspace
distMeanLog = (sum([samps1; samps2]))/(length([samps1;samps2]));

%Calculate the mean in Linear SR space
expsamps1 = exp(samps1);
expsamps2 = exp(samps2);
distMeanLinear = (sum([expsamps1; expsamps2]))/length([expsamps1; expsamps2]);

%Calculate Variance of mixed distributions in log and linear SR space
distVarLog = var([samps1;samps2]);
distVarLinear = var([expsamps1; expsamps2]);

%% Calculate mean and variance using sampling directly from Log Norm distribution to convince myself I get same answer
% %create samples from each lognormal dist
% lnsamps1 = lognrnd(mu_fit(1), stddev_fit(1), round(mix_fit(1)*numsamps), 1);
% lnsamps2 = lognrnd(mu_fit(2), stddev_fit(2), round(mix_fit(2)*numsamps), 1);
% meanOfLogNormalSamples = (sum([lnsamps1;lnsamps2])/length([lnsamps1;lnsamps2]));
% varLogNormalSamples = (var([lnsamps1; lnsamps2]));

%% Find the means and variances of the underlying data
dataMeanLinear = mean(data_linear);
dataVarLinear = var(data_linear);
dataMeanLog = mean(data_log);
dataVarLog = var(data_log);

 %% Plot to see what all the samples I've created look like
% figure
% subplot(2,2,1)
% histogram([samps1; samps2])
% xlabel("Log SR")
% title("Mean = " + num2str(meanInLogSpace) + "; Var = " + num2str(varInLogSpace))
% subplot(2,2,2)
% qqplot([samps1; samps2])
% xlabel("")
% subplot(2,2,[3 4])
% yyaxis left
% histogram([expsamps1; expsamps2])
% yyaxis right
% plot(x, z)
% xlabel("SR")
% xlim([0 10])
% title("Mean = " + num2str(meanExpOfGaussSamples) + "; Var = " + num2str(varExpOfGaussSamples))

%% Plot to see how data compare to estimated distributions
% figure
% subplot(2,2,1)
% yyaxis left
% histogram(data_log)
% yyaxis right
% plot(mixNormPDF(:,1), mixNormPDF(:,2))
% xlabel("Log SR")
% 
% title({"Data Mean = " + num2str(dataMeanLog) + "; Data Var = " + num2str(dataVarLog)...
%     , "Dist Mean = " + num2str(distMeanLog) + "; Dist Var = " + num2str(distVarLog)})
% subplot(2,2,2)
% qqplot(data_log)
% xlabel("Log SR Data")
% subplot(2,2,[3 4])
% yyaxis left
% histogram(data_linear)
% yyaxis right
% plot(x, z)
% xlabel("SR")
% xlim([0 10])
% title({"Data Mean = " + num2str(dataMeanLinear) + "; Data Var = " + num2str(dataVarLinear)...
%     , "Dist Mean = " + num2str(distMeanLinear) + "; Dist Var = " + num2str(distVarLinear)})
end