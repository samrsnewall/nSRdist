function[GamPDF, fitStruct] = fitGamma(data_linear, x, numObs)
% fitGamma  Fit a Gamma distribution to a dataset and evaluate its PDF.
%
% Fits a two-parameter Gamma distribution to data_linear using MATLAB's
% gamfit (maximum likelihood estimation), then evaluates the fitted PDF on
% the supplied x grid. If weighted replicates are used as input, the
% weighting must be applied before calling this function (i.e. data_linear
% should already be the replicated dataset).
%
% AIC, BIC, and a replicate-corrected BIC (BICtaeheefix) are computed.
% BICtaeheefix divides the NLL by the replication factor
% (length(data_linear)/numObs) before computing BIC, so that model
% selection is based on the true number of observations rather than the
% inflated replicate count.
%
% INPUTS
%   data_linear - (numeric vector) Data on the linear (untransformed) scale.
%                 May be a weighted-replicate dataset; see makeWeightedReplicates.
%   x           - (numeric vector) Evaluation points (linear scale) at which
%                 to compute the fitted Gamma PDF
%   numObs      - (scalar) True number of observations before any replication,
%                 used to compute BICtaeheefix
%
% OUTPUTS
%   GamPDF    - (numeric vector) Fitted Gamma PDF evaluated at each point
%               in x; same length as x
%   fitStruct - (struct) Fitting results:
%                 .alpha                  Shape parameter
%                 .beta                   Scale parameter
%                 .NumParams              Number of free parameters (2)
%                 .NegativeLogLikelihood  NLL of the fit
%                 .AIC                    Akaike Information Criterion
%                 .BIC                    Bayesian Information Criterion
%                 .BICtaeheefix           BIC corrected for replicate inflation
%
% See also: fitInvGamma, fitMixLogNorm, makeWeightedReplicates, ARfitdists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run gmfit
options = statset('Display', 'off', 'MaxIter', 200, 'TolFun', 1e-6);       %Set the number of iterations (default=100 often doesn't converge)
if ~iscolumn(data_linear)                                                   %Ensure the data is in column vector form
    data_linear = data_linear';
end
phat = gamfit(data_linear, "Options", options);                               %Fit gamma distribution to data

alpha = phat(1);
beta = phat(2);

%Create gamma pdf
GamPDF = gampdf(x, alpha, beta)';

%Find likelihood of fit distribution
nll = gamlike([alpha, beta], data_linear);

%Calculate NLL and correct for number of repetitions, as suggested by
%Taehee (divide NLL by number of repetitions)

info_divisor = length(data_linear)./numObs;
nll_taeheeFix = nll./info_divisor;
BIC_taeheeFix = 2*nll_taeheeFix + 2*log(numObs);

%Set up output structure with important values
fitStruct.alpha = alpha;
fitStruct.beta = beta;
fitStruct.NumParams = 2;
fitStruct.NegativeLogLikelihood = nll;
fitStruct.AIC = 2*fitStruct.NumParams + 2*nll;                             % AIC = 2k - 2ln(L) therefore = 2k +2*nll
fitStruct.BIC = log(length(data_linear))*fitStruct.NumParams + 2*nll;
fitStruct.BICtaeheefix = BIC_taeheeFix;



