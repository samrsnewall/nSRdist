function[GamPDF, phat] = fitGamma(data_linear, x)
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


