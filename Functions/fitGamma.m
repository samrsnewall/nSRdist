function[GamPDF, fitStruct] = fitGamma(data_linear, x, numObs)
%%% Create a Gamma to fit some dataset (if weighted, this must
%%% already be applied to data). 

%%% INPUT VALUES
% data_linear   = data on linear scale
% x             = x values on linear scale for which you want to know the pdf
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



