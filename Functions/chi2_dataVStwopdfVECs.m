function [h1, p1, chiStat1, h2, p2, chiStat2] = chi2_dataVStwopdfVECs(data, desiredSum, binN, pdfVEC1, pdfVEC2, fitS)
    %Define some bins
binEdges    = linspace(min(data), max(data), binN+1);
binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

%Find how many observations are in each bin
obsCounts     = histcounts(data, binEdges);
numObsCounts  = sum(obsCounts);
divisor       = numObsCounts/desiredSum;
obsCountsDown = round(obsCounts./divisor);

%Get the expected counts given pdfVEC1
%Set up some preliminary bin edges
interiorEdges = binEdges(2:end-1);
cumsumpdf1 = cumsum(pdfVEC1.fx);
pdfVEC1.cdf_x = cumsumpdf1./max(cumsumpdf1);
cdf1_at_edges = interp1(pdfVEC1.x, pdfVEC1.cdf_x, interiorEdges);
expProbs1 = diff([0, cdf1_at_edges, 1]);
expCounts1 = expProbs1.*desiredSum;

%Get the expected counts, given pdfVEC2
cumsumpdf2 = cumsum(pdfVEC2.fx);
pdfVEC2.cdf_x = cumsumpdf2./max(cumsumpdf2);
cdf2_at_edges = interp1(pdfVEC2.x, pdfVEC2.cdf_x, interiorEdges);
expProbs2 = diff([0, cdf2_at_edges, 1]);
expCounts2 = expProbs2.*desiredSum;

fitS.enforceBinSizeLimits = true;
if fitS.enforceBinSizeLimits
    minCountNumber = 5;
    maxCountNumber = desiredSum/3;
    [interiorEdges] = iterateBinSize_2pdfs(expCounts1, pdfVEC1, expCounts2, pdfVEC2, interiorEdges, desiredSum, minCountNumber, maxCountNumber);
    binEdges = [min(data), interiorEdges, max(data)];
    binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

    %Find how many observations are in each bin
    obsCounts     = histcounts(data, binEdges);
    numObsCounts  = sum(obsCounts);
    divisor       = numObsCounts/desiredSum;
    obsCountsDown = round(obsCounts./divisor);

    cdf1_at_edges = interp1(pdfVEC1.x, pdfVEC1.cdf_x, interiorEdges);
    expProbs1 = diff([0, cdf1_at_edges, 1]);
    expCounts1 = expProbs1.*desiredSum;
    cdf2_at_edges = interp1(pdfVEC2.x, pdfVEC2.cdf_x, interiorEdges);
    expProbs2 = diff([0, cdf2_at_edges, 1]);
    expCounts2 = expProbs2.*desiredSum;
end

%Perform chi2gof
[h1, p1, chiStat1] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts1, 'Emin', 5, 'NParams', pdfVEC1.numParams);
[h2, p2, chiStat2] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts2, 'Emin', 5, 'NParams', pdfVEC2.numParams);

%If chistat has changed binedges, use those
if numel(chiStat1.edges) ~= numel(binEdges)
    disp('a')
    binEdges = chiStat1.edges;
    binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

    %Find how many observations are in each bin
    obsCounts     = histcounts(data, binEdges);
    numObsCounts  = sum(obsCounts);
    divisor       = numObsCounts/desiredSum;
    obsCountsDown = round(obsCounts./divisor);

    %Get the expected counts, given the hypothesised distribution
    %Set up some preliminary bin edges
    interiorEdges = binEdges(2:end-1);
    cdf1_at_edges = interp1(pdfVEC1.x, pdfVEC1.cdf_x, interiorEdges);
    expProbs1 = diff([0, cdf1_at_edges, 1]);
    expCounts1 = expProbs1.*desiredSum;
    cdf2_at_edges = interp1(pdfVEC2.x, pdfVEC2.cdf_x, interiorEdges);
    expProbs2 = diff([0, cdf2_at_edges, 1]);
    expCounts2 = expProbs2.*desiredSum;
end

%Plot figure, with function, observed histogram, expected histogram, pvalue
%and chistat value
if fitS.dispChi2 == true

    %Display Result of chi2gof to commandline
    if h1 == 0
        disp("chi2gof: Accept H0; p = " +  num2str(p1, 3) + "; chistat = " + num2str(chiStat1.chi2stat, 3))
    else
        disp("chi2gof: Reject H0 - p = " + num2str(p1, 3) + "; chistat = " + num2str(chiStat1.chi2stat, 3))
    end

    %Plot figure
    figure;
    subplot(2,1,1)
    hold on
    %Plot histograms
    yyaxis("left")
    hObs = histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName", "Observed Counts");
    hExp = histogram('BinCounts', expCounts1, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName","Expected Counts");
    if divisor == 1
        ylabel("Counts")
    else
        ylabel("Downweighted Counts")
    end

    if fitS.enforceBinSizeLimits == false
        %Plot function
        yyaxis("right")
        lTested = plot(pdfVEC1.x, pdfVEC1.fx, 'k', 'LineWidth', 2, "DisplayName", "Tested PDF"); %#ok<*NASGU>
        ylabel("PDF")
    end

    %Add pval and chi2stat to legend
    lpval   = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "p = " + num2str(p1, 3));
    lchi2stat = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "chi2stat = "  + num2str(chiStat1.chi2stat, 3));
    %xlabel("log nSR") Note, it's not always log nSR...
    legend()
    xlim([-4 4])
    subplot(2,1,2)
    hold on
    %Plot histograms
    yyaxis("left")
    hObs = histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName", "Observed Counts");
    hExp = histogram('BinCounts', expCounts2, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName","Expected Counts");
    if divisor == 1
        ylabel("Counts")
    else
        ylabel("Downweighted Counts")
    end

    if fitS.enforceBinSizeLimits == false
        %Plot function
        yyaxis("right")
        lTested = plot(pdfVEC2.x, pdfVEC2.fx, 'k', 'LineWidth', 2, "DisplayName", "Tested PDF"); %#ok<*NASGU>
        ylabel("PDF")
    end

    %Add pval and chi2stat to legend
    lpval2   = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "p = " + num2str(p2, 3));
    lchi2stat2 = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "chi2stat = "  + num2str(chiStat2.chi2stat, 3));
    %xlabel("log nSR") Note, it's not always log nSR...
    legend()
    xlim([-4 4])

end

end