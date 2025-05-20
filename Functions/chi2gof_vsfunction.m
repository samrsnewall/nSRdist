function[h,p,chiStat] = chi2gof_vsfunction(data, cdfFH, numParams, desiredSum, binN, fitS)
%% Chi2 GOF test of data against given distribution

%Define some bins
binEdges    = linspace(min(data), max(data), binN+1);
binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

%Find how many observations are in each bin
obsCounts     = histcounts(data, binEdges);
numObsCounts  = sum(obsCounts);
divisor       = numObsCounts/desiredSum;
obsCountsDown = round(obsCounts./divisor);

%Get the expected counts, given the hypothesised distribution
%Set up some preliminary bin edges
interiorEdges = binEdges(2:end-1);
expProbs = diff([0, cdfFH(interiorEdges), 1]);
expCounts = expProbs.*desiredSum;

if fitS.enforceBinSizeLimits
    minCountNumber = fitS.chi2MinCountNum;
    maxCountNumber = 50;
    [interiorEdges] = iterateBinSize(expCounts, cdfFH, interiorEdges, desiredSum, minCountNumber, maxCountNumber);
    binEdges = sort([min(data), interiorEdges, max(data)]);
    binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

    %Find how many observations are in each bin
    obsCounts     = histcounts(data, binEdges);
    numObsCounts  = sum(obsCounts);
    divisor       = numObsCounts/desiredSum;
    obsCountsDown = round(obsCounts./divisor);

    expProbs = diff([0, cdfFH(interiorEdges), 1]);
    expCounts = expProbs.*desiredSum;
end

%Perform chi2gof
[h, p, chiStat] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts, 'Emin', 5, 'NParams', numParams);

%If chistat has changed binedges, use those
if numel(chiStat.edges) ~= numel(binEdges)
    binEdges = chiStat.edges;
    binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

    %Find how many observations are in each bin
    obsCounts     = histcounts(data, binEdges);
    numObsCounts  = sum(obsCounts);
    divisor       = numObsCounts/desiredSum;
    obsCountsDown = round(obsCounts./divisor);

    %Get the expected counts, given the hypothesised distribution
    %Set up some preliminary bin edges
    interiorEdges = binEdges(2:end-1);
    expProbs = diff([0, cdfFH(interiorEdges), 1]);
    expCounts = expProbs.*desiredSum;

end

%Get pdf from cdf by differentiating
step = 1e-5;
pdfFH = @(t) (cdfFH(t+step) - cdfFH(t-step))/(2*step);

%Define x between the 0.01 and 99.99 percentile of distribution
lowpct = fzero(@(t) cdfFH(t) - 0.0001, -5);
highpct = fzero(@(t) cdfFH(t) - 0.9999, 5);
x = linspace(lowpct, highpct, 100);

%Plot figure, with function, observed histogram, expected histogram, pvalue
%and chistat value
if fitS.dispChi2 == true

    %Display Result of chi2gof to commandline
    if h == 0
        disp("chi2gof: Accept H0; p = " +  num2str(p, 3) + "; chistat = " + num2str(chiStat.chi2stat, 3))
    else
        disp("chi2gof: Reject H0 - p = " + num2str(p, 3) + "; chistat = " + num2str(chiStat.chi2stat, 3))
    end

    %Plot figure
    figure;
    hold on

    %Plot histograms
    yyaxis("left")
    hObs = histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName", "Observed Counts");
    hExp = histogram('BinCounts', expCounts, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName","Expected Counts");
    if divisor == 1
        ylabel("Counts")
    else
        ylabel("Downweighted Counts")
    end

    if fitS.enforceBinSizeLimits == false
        %Plot function
        yyaxis("right")
        lTested = plot(x, pdfFH(x), 'k', 'LineWidth', 2, "DisplayName", "Tested PDF"); %#ok<*NASGU>
        ylabel("PDF")
    end

    %Add pval and chi2stat to legend
    lpval   = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "p = " + num2str(p, 3));
    lchi2stat = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "chi2stat = "  + num2str(chiStat.chi2stat, 3));
    %xlabel("log nSR") Note, it's not always log nSR...
    legend()
end
end