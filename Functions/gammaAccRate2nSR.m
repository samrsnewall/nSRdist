function[nSR, nSRprob] = gammaAccRate2nSR(accShape, accScale)
%This function takes in some value of alpha (shape factor) for a gamma
%distribution, creates the resulting distribution and converts its units
%from accumulation rate (yr/cm - as used in Bacon) to normalised sed rate
%(as used in BIGMACS). It outputs the nSR values and their associated
%probability density values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Use shape factor to get scale factor
% accMeanYr = 200; %(yr/cm) Note, accMeanYr has no impact on nSR distribution, so value is somewhat arbitrary
% accMean = accMeanYr/1000; %(kyr/cm)
% accScale = accMean/accShape;

%%Use shape and scale factor to get accMean
accMean = accScale*accShape;

%Find y values for Gamma distribution on normalised SR
x = 0:(accMean/500):accMean*60; %make x scale sufficiently large and high resolution based on accMean (This doesn't quite do its job - works for this given accMean)
yGamma = gampdf(x, accShape, accScale);

%Convert to SR from invSR
[SR, SRprob] = invSRtoSR(x,yGamma);

%Convert to normalised SR
nSR = SR/(1/accMean);
areaUnderCurve = trapz(nSR, SRprob);
nSRprob = SRprob/areaUnderCurve;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end