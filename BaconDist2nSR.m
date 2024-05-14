%Use Bacon defaults to get shape and scale factors
accShape = 1.5;
accMeanYr = 200; %(yr/cm)
accMean = accMeanYr/1000; %(kyr/cm)
accScale = accMean/accShape;
defaultGamma = makedist("Gamma", "a", accShape, "b", accScale);

%Find y values for Gamma distribution on normalised SR
x = 0:(accMean/50):accMean*60;
yGamma = gampdf(x, accShape, accScale);

%Convert to SR from invSR
[SR, SRprob] = invSRtoSR(x,yGamma);

%Convert to normalised SR
nSR = SR/(1/accMean);

% Plot these distributions
figure()
subplot(3,1,1)
plot(x,yGamma)
ylabel("probability density")
xlabel("invSR (kyr/cm)")
subplot(3,1,2)
plot(SR, SRprob)
ylabel("Probability Density")
xlabel("SR (cm/kyr)")
subplot(3,1,3)
plot(nSR, SRprob)
ylabel("Probability Density")
xlabel("nSR")