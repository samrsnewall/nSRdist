function[h,p,chiStat] = chi2gof_vsfunction(data, cdfFH, desiredSum, binN, fitS)
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
[h, p, chiStat] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts, 'Emin', 1);

%Get pdf from cdf by differentiating
step = 1e-5;
pdfFH = @(t) (cdfFH(t+step)-cdfFH(t-step))/(2*step);
x = linspace(min(data), max(data), 100);

%Plot
if fitS.dispChi2 == true;
    if h == 0
        disp("chi2gof: Accept H0; p = " +  num2str(p, 3) + "; chistat = " + num2str(chiStat.chi2stat, 3))
    else
        disp("chi2gof: Reject H0 - p = " + num2str(p, 3) + "; chistat = " + num2str(chiStat.chi2stat, 3))
    end

    figure;
    hold on

    yyaxis("left")
    hObs = histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName", "Observed Counts");
    hExp = histogram('BinCounts', expCounts, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName","Expected Counts");
    if divisor == 1
        ylabel("Counts")
    else
        ylabel("Downweighted Counts")

    yyaxis("right")
    lTested = plot(x, pdfFH(x), 'k', 'LineWidth', 2, "DisplayName", "Tested PDF");
    ylabel("PDF")

    lpval   = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "p = " + num2str(p, 3));
    lchi2stat = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "chi2stat = "  + num2str(chiStat.chi2stat, 3));
    xlabel("log nSR")
    legend()
end
end