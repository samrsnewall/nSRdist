function[invGamPDF, GamPDF, fitStruct] = fitInvGamma(data_linear, x, numObs)
% fitInvGamma  Fit an Inverse Gamma distribution to a dataset and evaluate
%              its PDF.
%
% Fits an Inverse Gamma distribution by applying the substitution
% y = 1/x and fitting a Gamma distribution to the reciprocal data using
% MATLAB's gamfit (MLE). The Gamma PDF on 1/x space is then transformed
% back to x space via px_to_pfx, which applies the Jacobian correction
% for the change of variables. The NLL includes a Jacobian correction term
% (+ 2*sum(log(data_linear))) to express likelihood in the original x space.
%
% If weighted replicates are used as input, the weighting must be applied
% before calling this function (i.e. data_linear should already be the
% replicated dataset).
%
% AIC, BIC, and a replicate-corrected BIC (BICtaeheefix) are computed.
% BICtaeheefix divides the NLL by the replication factor
% (length(data_linear)/numObs) before computing BIC, so that model
% selection is based on the true number of observations rather than the
% inflated replicate count.
%
% INPUTS
%   data_linear - (numeric vector) Data on the linear (untransformed) scale,
%                 i.e. without the inverse applied. May be a weighted-replicate
%                 dataset; see makeWeightedReplicates.
%   x           - (numeric vector) Evaluation points (linear scale) at which
%                 to compute the fitted Inverse Gamma PDF
%   numObs      - (scalar) True number of observations before any replication,
%                 used to compute BICtaeheefix
%
% OUTPUTS
%   invGamPDF - (numeric vector) Fitted Inverse Gamma PDF evaluated at each
%               point in x (linear scale); same length as x
%   GamPDF    - (numeric vector) Fitted Gamma PDF evaluated on the
%               reciprocal axis 1/x; same length as x. Intermediate output,
%               useful for diagnostic checking.
%   fitStruct - (struct) Fitting results (parameters refer to the fitted
%               Gamma on 1/x):
%                 .alpha                  Shape parameter
%                 .beta                   Scale parameter
%                 .NumParams              Number of free parameters (2)
%                 .NegativeLogLikelihood  NLL in original x space (with
%                                         Jacobian correction applied)
%                 .AIC                    Akaike Information Criterion
%                 .BIC                    Bayesian Information Criterion
%                 .BICtaeheefix           BIC corrected for replicate inflation
%
% See also: fitGamma, fitMixLogNorm, px_to_pfx, makeWeightedReplicates, ARfitdists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run gmfit
options = statset('Display', 'off', 'MaxIter', 200, 'TolFun', 1e-6);       %Set the number of iterations (default=100 often doesn't converge)
if ~iscolumn(data_linear)                                                  %Ensure the data is in column vector form
    data_linear = data_linear';
end
data_inv = 1./data_linear;                                                 %Take inverse of data
phat = gamfit(data_inv, "Options", options);                               %Fit gamma distribution to inverse of data

alpha = phat(1);
beta = phat(2);

% Create gamma on invx values
invx = 1./x;
GamPDF = gampdf(invx, alpha, beta)';

%Find likelihood of fit distribution
nll = gamlike([alpha, beta], data_inv) + 2*sum(log(data_linear));               %Has jacobian correction

%Find pdf vector in x space
f = @(x) 1./x;
[newx, invGamPDF_toInterpolate] = px_to_pfx(invx, GamPDF, f);

%Interpolate back to desired x vector spacing
invGamPDF = interp1(round(newx, 10), invGamPDF_toInterpolate, round(x, 10))';

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

% Check on likelihood
GamPDF1 = @(y) gampdf(y, alpha, beta);
invGamPDF1 = @(x) gampdf(1./x, alpha, beta) .* (1./x.^2);
x = data_linear;
nll1 = -sum(log(invGamPDF1(x)));
