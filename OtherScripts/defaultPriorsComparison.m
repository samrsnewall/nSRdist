%Set up the Bacon default gamma distribution on inverse nSR
alpha = 1.5;
bacondefault = makedist("Gamma", "a", alpha, "b", 1./alpha);
xvals = 0:0.01:100;
bacondefaultline = pdf(bacondefault, xvals);

%Convert the default Bacon gamma to nSR
[SRx, BD_SR] = invSRtoSR(xvals, bacondefaultline);

%Load in the BIGMACS distribution
BMlog = load("lognormal_BIGMACS.txt");
BM.SRx = BMlog(:,1);
BM.pSR = BMlog(:,2);

%Convert the BIGMACS distribution to inverse nSR
[BM.invSRx, BM.pinvSR] = invSRtoSR(BM.SRx, BM.pSR);

%Plot both distributions on nSR
figure;
%plot(xvals, bacondefaultline)
plot(SRx, BD_SR)
hold on;
plot(BM.SRx, BM.pSR)
xlim([0 6])


% Plot both distributions on inverse nSR
figure;
plot(xvals, bacondefaultline)
hold on
plot(BM.invSRx, BM.pinvSR)
xlim([0 6])