function[invGamPDF, GamPDF, fitStruct] = fitInvGamma(data_linear, x, numObs)
%%% Create a Inverse Gamma to fit some dataset (if weighted, this must
%%% already be applied to data).

%%% INPUT VALUES
% data_linear   = data on linear scale (i.e. without inverse applied)
% x             = x values on linear scale for which you want to know the pdf
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
