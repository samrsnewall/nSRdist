function[h,p,chistats] = chi2gof_vsfunction(data, cdfFH, desiredSum, binN)
%% Chi2 GOF test of data against given distribution

%Define some bins
binEdges    = linspace(min(data), max(data), binN+1);
binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

%Find how many observations are in each bin
obsCounts     = histcounts(data, binEdges);
numObsCounts  = sum(obsCounts);
divisor       = numObsCounts./desiredSum;
obsCountsDown = round(obsCounts./divisor);

%Get the expected counts, given the hypothesised distribution
expProbs = cdfFH(binEdges(2:end)) - cdfFH(binEdges(1:end-1));
expCounts = expProbs.*sum(obsCountsDown);

%Perform chi2gof
[h, p, chistats] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts, 'Emin', 1);
if h == 0;
    disp("chi2gof: Accept H0 - p = " +  num2str(p))
else
    disp("chi2gof: Reject H0 - p = " + num2str(p))
end

%
%Get pdf numerically from cdf to use in plots
%THE PDF IS NOT COMING OUT CORRECTLY! FIX THIS
step = 1e-5;
pdfFH = @(t) (cdfFH(t+step)-cdfFH(t-step))/(2*step);
x = linspace(min(data), max(data), 100);

%Plot
figure;
hold on
yyaxis("left")
histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName", "Observed Counts")
histogram('BinCounts', expCounts, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName","Expected Counts")
ylabel("Counts")
yyaxis("right")
plot(x, pdfFH(x), 'k', 'LineWidth', 2, "DisplayName", "Fitted PDF")
ylabel("PDF")
xlabel("log nSR")
legend()
end