function[mixLogNormPDF, mixNormPDF, mlnfit] = fitMixLogNorm(data_linear, x, numComponents, replicates, numObs)
% fitMixLogNorm  Fit a Mixture Log-Normal distribution to a dataset and
%                evaluate its PDF.
%
% Fits a numComponents-component Mixture Log-Normal (MLN) by taking the
% logarithm of data_linear and fitting a Gaussian Mixture Model (GMM) in
% log space using MATLAB's fitgmdist. The component means and standard
% deviations from the GMM are then used to construct the corresponding
% log-normal PDFs in the original linear space, which are combined with
% their mixing proportions to give the final MLN PDF.
%
% Passing numComponents = 1 fits a single Log-Normal (LN) distribution.
%
% If weighted replicates are used as input, the weighting must be applied
% before calling this function (i.e. data_linear should already be the
% replicated dataset).
%
% BIC and a replicate-corrected BIC (BICtaeheefix) are computed.
% BICtaeheefix divides the NLL by the replication factor
% (length(data_linear)/numObs) before computing BIC, so that model
% selection is based on the true number of observations rather than the
% inflated replicate count. The NLL includes a Jacobian correction term
% (+ sum(log(data_linear))) to express likelihood in the original linear
% space rather than log space.
%
% INPUTS
%   data_linear   - (numeric vector) Data on the linear (untransformed) scale.
%                   May be a weighted-replicate dataset; see makeWeightedReplicates.
%   x             - (numeric vector) Evaluation points (linear scale) at which
%                   to compute the fitted MLN PDF
%   numComponents - (integer) Number of log-normal components in the mixture;
%                   use 1 for a single Log-Normal, 2 for a two-component MLN
%   replicates    - (integer) Number of random initialisations passed to
%                   fitgmdist ('Replicates' option); higher values reduce
%                   the risk of converging to a local optimum
%   numObs        - (scalar) True number of observations before any replication,
%                   used to compute BICtaeheefix
%
% OUTPUTS
%   mixLogNormPDF - (N × 2 numeric matrix) Fitted MLN PDF in linear space:
%                   column 1 is x (evaluation points), column 2 is the
%                   probability density at each point
%   mixNormPDF    - (N × 2 numeric matrix) Fitted GMM PDF in log space:
%                   column 1 is log(x), column 2 is the probability density.
%                   Useful for diagnostic checking of the fit in log space.
%   mlnfit        - (struct) Fitting results:
%                     .NumParams              Number of free parameters
%                                             (2 for LN; 2 + 3*(K-1) for K-component MLN)
%                     .mu                     Component means in log space
%                     .Sigma                  Component variances in log space
%                     .ComponentProportion    Mixing weights (sum to 1)
%                     .NegativeLogLikelihood  NLL in linear space (with
%                                             Jacobian correction applied)
%                     .BIC                    Bayesian Information Criterion
%                     .BICtaeheefix           BIC corrected for replicate inflation
%
% See also: fitGamma, fitInvGamma, makeWeightedReplicates, ARfitdists, IRfitdists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run gmfit
options = statset('Display', 'off', 'MaxIter', 1000, 'TolFun', 1e-6);       %Set the number of iterations (default=500 often doesn't converge)
if ~iscolumn(data_linear)                                                   %Ensure the data is in column vector form
    data_linear = data_linear';
end
data_log = log(data_linear);                                                %Take logarithm of data
gmfit    = fitgmdist(data_log, numComponents, "Options", options, ...       %Fit mix normal to logarithm of data
     'Replicates', replicates, 'Start','randSample');


%Initialise std deviation vector and prob density vector
stddev_fit  = NaN(1,numComponents);
y           = zeros(length(x), numComponents);
lx          = log(x);
ly          = zeros(length(lx), numComponents);

%Pull out sttdev of each component gaussian, and then find the lognormal
%pdf using the mean and stddev of the component gaussians as the input
%parameters.
mu_fit = gmfit.mu;                                                          % Find the means of fit
for i = 1:numComponents
    stddev_fit(i)   = sqrt(gmfit.Sigma(:,:,i));                             % Taking the square converts the variance to a stddev
    y(:,i)          = lognpdf(x, mu_fit(i), stddev_fit(i));                 % Find the y values of a lognormal using the mean and sttdeviations taken from the mixed gaussian
    ly(:,i)         = normpdf(lx, mu_fit(i), stddev_fit(i));                % Find the y values of a Normal using mean and sttdeviations...
end


%Find the mixed log normal, summing the individual gaussians after
%multiplying by their mixing component
mix_fit         = gmfit.ComponentProportion;                                %Find the mixing component of each individual Gaussian
z               = sum(y.*mix_fit, 2);
if find(z == max(z)) == 1
    disp('the max value is the first index... This is suspicious, double check if this is meant to be the case')
end
lz              = sum(ly.*mix_fit, 2);
mixLogNormPDF   = [x', z];
mixNormPDF      = [lx', lz];

% Calculate the negative log likelihood in linear space
numParams = 2+(numComponents-1)*3;
nll_mixlogn = gmfit.NegativeLogLikelihood + sum(log(data_linear));          %addition of sum(log(data_linear)) is the Jacobian correction to convert back to linear space
BIC_mixlogn = 2*nll_mixlogn + numParams*log(length(data_linear));

%Calculate NLL and correct for number of repetitions, as suggested by
%Taehee (divide NLL by number of repetitions)

info_divisor = length(data_linear)./numObs;
nll_taeheeFix = nll_mixlogn./info_divisor;
BIC_taeheeFix = 2*nll_taeheeFix + numParams*log(numObs);


%Set up output structure
mlnfit.NumParams = numParams;
mlnfit.mu = gmfit.mu;
mlnfit.Sigma = gmfit.Sigma;
mlnfit.ComponentProportion = gmfit.ComponentProportion;
mlnfit.NegativeLogLikelihood = nll_mixlogn;
mlnfit.BIC = BIC_mixlogn;
mlnfit.BICtaeheefix = BIC_taeheeFix;
end