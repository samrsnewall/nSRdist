function[h, p, chiStat] = chi2gof_vsMLN(gmfit, dataLog, numSRcalcs, fitS)

%Create cdf of distribution
%Get the parameters of the mixed gaussian in log space
mu1 = gmfit.mu(1); sigma1 = sqrt(gmfit.Sigma(1));
mu2 = gmfit.mu(2); sigma2 = sqrt(gmfit.Sigma(2));
w1  = gmfit.ComponentProportion(1);
w2  = gmfit.ComponentProportion(2);

%Create the cdf function handle
mlncdf = @(t) w1 * normcdf(t, mu1, sigma1) + w2 * normcdf(t, mu2, sigma2);

%Set up certain important parameters
%desiredSum = length(dataLog);
desiredSum = numSRcalcs;
binN = fitS.chi2binN;
numParams = 6;

%Perform chi2gof
[h,p,chiStat] = chi2gof_vsfunction(dataLog, mlncdf, numParams, desiredSum, binN, fitS);
end