function[mixLogNormPDF, mixNormPDF, mlnfit] = fitMixLogNorm(data_linear, x, numComponents, replicates, numObs)
%%% Create a mixed log normal to fit some dataset (if weighted, this must
%%% already be applied to data). This function applies a regularization
%%% value. A regularization value of 0 is the same as not regularising.

%%% INPUT VALUES
% data_linear   = data on linear scale (i.e. without logarithm applied)
% x             = x values on linear scale for which you want to know the pdf
% numComponents = number of components in mixture Log Normal desired
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
nll_mixlogn = gmfit.NegativeLogLikelihood + sum(log(data_linear));          %addition of sum(log(data_linear)) is the Jacobian correction to convert back to linear space
BIC_mixlogn = 2*nll_mixlogn + numComponents*log(length(data_linear));

%Calculate NLL and correct for number of repetitions, as suggested by
%Taehee (divide NLL by number of repetitions)

    info_divisor = length(data_linear)./numObs;
    nll_taeheeFix = nll_mixlogn./info_divisor;
    BIC_taeheeFix = 2*nll_taeheeFix + numComponents*log(numObs);


%Set up output structure
mlnfit.NumVariables = gmfit.NumVariables;
mlnfit.mu = gmfit.mu;
mlnfit.Sigma = gmfit.Sigma;
mlnfit.ComponentProportion = gmfit.ComponentProportion;
mlnfit.NegativeLogLikelihood = nll_mixlogn;
mlnfit.BIC = BIC_mixlogn;
mlnfit.BICtaeheefix = BIC_taeheeFix;
end