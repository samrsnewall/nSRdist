function[pdf2save, meanExpOfGaussSamples, varExpOfGaussSamples] = fitMixLogNorm(data_linear, numComponents)
%%% Create a mixed log normal to fit some dataset

%Set the range of x values of interest (use x values currently used in
%BIGMACS)
lognorm_BIGMACS = readtable("lognormal_BIGMACS.txt");
lx = lognorm_BIGMACS.Var1';

%Set the number of iterations (default=500 often doesn't converge)
options = statset('MaxIter', 2000);

%Use GMFit to fit a Gaussian Mixture distribution to the log of the data
if ~iscolumn(data_linear)
    data_linear = data_linear';
end
data_log = log(data_linear);
gmfit = fitgmdist(data_log, numComponents, "Options", options);

%Find the means of fit
mu_fit = gmfit.mu;

%Initialise std deviation vector and prob density vector
stddev_fit = NaN(1,numComponents);
ly = zeros(length(lx), numComponents);

%Pull out sttdev of each component gaussian, and then find the lognormal
%pdf using the mean and stddev of the component gaussians as the input
%parameters.
for i = 1:numComponents
    stddev_fit(i) = sqrt(gmfit.Sigma(:,:,i)); %taking the square converts the variance to a stddev
    ly(:,i) = lognpdf(lx, mu_fit(i), stddev_fit(i)); % Find the y values of a lognormal using the mean and sttdeviations taken from the mixed gaussian
end

%Find the mixing component of each individual Gaussian
mix_fit = gmfit.ComponentProportion;

%Find the mixed log normal, summing the individual gaussians after
%multiplying by their mixing component
lz = sum(ly.*mix_fit, 2);
pdf2save = [lx', lz];

%% Find the mean and variance of the log normal distribution
%meanAnalytical = mix_fit(1)*exp(mu_fit(1) + (stddev_fit(1).^2)/2) + mix_fit(2) * exp(mu_fit(2) + (stddev_fit(2).^2)/2)

%create samples from each gaussian dist
numsamps = 10000;
samps1 = normrnd(mu_fit(1), stddev_fit(1), round(mix_fit(1)*numsamps), 1);
samps2 = normrnd(mu_fit(2), stddev_fit(2), round(mix_fit(2)*numsamps), 1);
%Calculate the mean of the Gaussian Distribution in LogSRspace
meanInLogSRspace = (sum([samps1; samps2]))/(length([samps1;samps2]));

%Calculate the mean in Linear SR space
expsamps1 = exp(samps1);
expsamps2 = exp(samps2);
meanExpOfGaussSamples = (sum([expsamps1; expsamps2]))/length([expsamps1; expsamps2]);

%Calculate Variance of mixed distributions in log and linear SR space
varInLogSRspace = var([samps1;samps2]);
varExpOfGaussSamples = var([expsamps1; expsamps2]);

%% Calculate mean and variance using sampling directly from Log Norm distribution to convince myself I get same answer
% %create samples from each lognormal dist
% lnsamps1 = lognrnd(mu_fit(1), stddev_fit(1), round(mix_fit(1)*numsamps), 1);
% lnsamps2 = lognrnd(mu_fit(2), stddev_fit(2), round(mix_fit(2)*numsamps), 1);
% meanOfLogNormalSamples = (sum([lnsamps1;lnsamps2])/length([lnsamps1;lnsamps2]));
% varLogNormalSamples = (var([lnsamps1; lnsamps2]));

%% Plot to see what all the samples I've created look like
% figure
% subplot(3,1,1)
% histogram([samps1; samps2])
% xlabel("Log SR")
% title("Mean = " + num2str(meanInLogSRspace) + "; Var = " + num2str(varInLogSRspace))
% subplot(3,1,2)
% histogram([expsamps1; expsamps2])
% xlabel("SR")
% xlim([0 10])
% title("Mean = " + num2str(meanExpOfGaussSamples) + "; Var = " + num2str(varExpOfGaussSamples))
% subplot(3,1,3)
% histogram([lnsamps1; lnsamps2])
% xlabel("SR")
% xlim([0 10])
% title("Mean = " + num2str(meanOfLogNormalSamples) + "; Var = " + num2str(varLogNormalSamples))
end