function[pdf2save] = fitMixLogNorm(data_linear)
%%% Create a mixed log normal from some dataset
% Choose number of components
numcomps = 2;
%Set the range of x values of interest
%lx = -3:0.02:10;
lognorm_BIGMACS = readtable("../lognormal_BIGMACS.txt");
lx = lognorm_BIGMACS(:,1)';
%Set the number of iterations (default often doesn't converge)
options = statset('MaxIter', 1000);
%Use GMFit to fit a Gaussian Mixture distribution to the data (log of data
%used)
data_log = log(data_linear);
gmfit = fitgmdist(data_log, numcomps, "Options", options);
%Find the means of fit
mu_fit = gmfit.mu;
%Initialise std deviation vector and prob density vector
stddev_fit = NaN(1,numcomps);
ly = zeros(length(lx), numcomps);
%Pull out sttdev of each component gaussian, and then find the lognormal
%pdf using the mean and stddev of the component gaussians as the input
%parameters.
for i = 1:numcomps
    stddev_fit(i) = sqrt(gmfit.Sigma(:,:,i)); %taking the square converts the variance to a stddev
    ly(:,i) = lognpdf(lx, mu_fit(i), stddev_fit(i)); % Find the y values of a lognormal using the mean and sttdeviations taken from the mixed gaussian
end
%Find the mixing component of each individual Gaussian
mix_fit = gmfit.ComponentProportion;
%Find the mixed log normal, summing the individual gaussians after
%multiplying by their mixing component
lz = sum(ly.*mix_fit, 2);
pdf2save = [lx', lz];
end