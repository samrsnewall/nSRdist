function[mu, var] = muVarPDFVec(pdf)
% muVarPDFVec  Compute the mean and variance of a PDF defined on a struct.
%
% Uses the trapezoid rule to numerically integrate x*p(x) and x²*p(x)
% over the supplied x-axis, giving the mean and variance of the
% distribution. The input PDF does not need to be normalised, but results
% are only meaningful if it approximately integrates to 1.
%
% INPUT
%   pdf  - (struct) with fields:
%            .x   (numeric vector) Evaluation points
%            .px  (numeric vector) Probability density values at each x
%
% OUTPUTS
%   mu   - (scalar) Mean of the distribution
%   var  - (scalar) Variance of the distribution
%
% See also: ARfitdists, IRfitdists, px_to_pfx

mu = trapz(pdf.x, pdf.x.*pdf.px);
var = trapz(pdf.x, pdf.x.*pdf.x.*pdf.px) - mu^2;
end