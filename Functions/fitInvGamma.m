function[invGamPDF, GamPDF, phat] = fitInvGamma(data_linear, x)
%%% Create a Inverse Gamma to fit some dataset (if weighted, this must
%%% already be applied to data). This function applies a regularization
%%% value. A regularization value of 0 is the same as not regularising.

%%% INPUT VALUES
% data_linear   = data on linear scale (i.e. without inverse applied)
% x             = x values on linear scale for which you want to know the pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run gmfit
options = statset('Display', 'off', 'MaxIter', 200, 'TolFun', 1e-6);       %Set the number of iterations (default=100 often doesn't converge)
if ~iscolumn(data_linear)                                                   %Ensure the data is in column vector form
    data_linear = data_linear';
end
data_inv = 1./data_linear;                                                %Take inverse of data
phat = gamfit(data_inv, "Options", options);                               %Fit gamma distribution to inverse of data

alpha = phat(1);
beta = phat(2);

% Create gamma on invx values
invx = 1./x;
GamPDF = gampdf(invx, alpha, beta)';

%Find pdf vector in x space
f = @(x) 1./x;
[newx, invGamPDF_toInterpolate] = px_to_pfx(invx, GamPDF, f);

%Interpolate back to desired x vector spacing
invGamPDF = interp1(round(newx, 10), invGamPDF_toInterpolate, round(x, 10))';

