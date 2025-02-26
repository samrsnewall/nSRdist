function [h1, p1, chiStat1, h2, p2, chiStat2] = chi2_dataVStwopdfVECs(data, desiredSum, binN, pdfVEC1, pdfVEC2, fitS)
%Define some initial bins
binEdges    = linspace(min(data), max(data), binN+1);
interiorEdges = binEdges(2:end-1);
binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

%Get cdf of pdfs
cumsumpdf1 = cumsum(pdfVEC1.px);
pdfVEC1.cdf_x = cumsumpdf1./max(cumsumpdf1);
cumsumpdf2 = cumsum(pdfVEC2.px);
pdfVEC2.cdf_x = cumsumpdf2./max(cumsumpdf2);

%Calculate observed and expected bin counts with given binEdges
[obsCountsDown, expCounts1, expCounts2] = calcBinCounts(data, binEdges, desiredSum, pdfVEC1, pdfVEC2);

fitS.enforceBinSizeLimits = true;
if fitS.enforceBinSizeLimits
    minCountNumber = 5;
    maxCountNumber = desiredSum/3;
    % tic
    % [interiorEdges] = iterateBinSize_2pdfs(expCounts1, pdfVEC1, expCounts2, pdfVEC2, interiorEdges, desiredSum, minCountNumber, maxCountNumber);
    % toc
    % binEdges = [min(data), interiorEdges, max(data)];
    % binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;
    % 
    % %Calculate observed and expected bin counts
    % [obsCountsDown, expCounts1, expCounts2] = calcBinCounts(data, binEdges, desiredSum, pdfVEC1, pdfVEC2);

    %APPROACH #2 - Start with N equally spaced bins
    startN = binN+1;
    binEdges3    = linspace(min(data), max(data), startN);
    [~, expCounts13, expCounts23] = calcBinCounts(data, binEdges3, desiredSum, pdfVEC1, pdfVEC2);

    %Check whether the bin counts are too high, and increase number of bins
    newN = startN;
    while sum(expCounts13 > maxCountNumber | expCounts23 > maxCountNumber) > 0
        newN = newN+1;
        binEdges3    = linspace(min(data), max(data), newN);
        [~, expCounts13, expCounts23] = calcBinCounts(data, binEdges3, desiredSum, pdfVEC1, pdfVEC2);
    end

    %Start to decrease number of bins, until there are too few
    while sum(expCounts13 < maxCountNumber | expCounts23 < maxCountNumber) == 0
        %Try with less bins
        newN = newN-1;
        binEdges4    = linspace(min(data), max(data), newN);
        [~, expCounts13, expCounts23] = calcBinCounts(data, binEdges4, desiredSum, pdfVEC1, pdfVEC2);
    end
    
    %Add one more bin number (highest number of bins that keeps all bins
    %below maxCountsNumber;
    if sum(expCounts13 > maxCountNumber | expCounts23 > maxCountNumber) > 0
        newN = newN+1;
        binEdges3    = linspace(min(data), max(data), startN);
        [~, expCounts13, expCounts23] = calcBinCounts(data, binEdges3, desiredSum, pdfVEC1, pdfVEC2);
    end

    %Need to check which bins are too small now
    %If first bin is too small
    if expCounts13(1) < minCountNumber || expCounts23(1) < minCountsNumber
        %Find what a first bin size that would have at least 5 counts
        %in both pdfs
        intEdgeMin = min(binEdges(2:end-1));
        if intEdgeMin >0
            intEdgeMinEst = intEdgeMin/10;
        else
            intEdgeMinEst = intEdgeMin*10;
        end
        minVals2Try = linspace(intEdgeMinEst, max(binEdges(2:end-1)), 100);
        expCountsMin1 = desiredSum.*(interp1(pdfVEC1.x, pdfVEC1.cdf_x, minVals2Try));
        expCountsMin2 = desiredSum.*(interp1(pdfVEC2.x, pdfVEC2.cdf_x, minVals2Try));
        minEdgeInd = find(...
            expCountsMin1 > minCountNumber &...
            expCountsMin1 < maxCountNumber &...
            expCountsMin2 > minCountNumber &...
            expCountsMin2 < maxCountNumber, 1 );
        if sum(minEdgeInd) <1
            figure()
            yyaxis left
            hold on
            plot(pdfVEC1.x, pdfVEC1.fx, 'r')
            plot(pdfVEC2.x, pdfVEC2.fx, 'b')
            error("No bin edges found that work for both distributions")
        end
        minPossBinEdge = minVals2Try(minEdgeInd);   %Find a possible bin edge that would provide at least 5 counts to minimum bin
        minBinEdge = min(binEdges3(binEdges3>minPossBinEdge)); %Choose even spaced bin edge that is next biggest than minimum possible bin edge
        binEdges3 = [binEdges3(1), binEdges3(binEdges3 >=minBinEdge)];
        [~,expCounts13, expCounts23] = calcBinCounts(data, binEdges3, desiredSum, pdfVEC1, pdfVEC2);
    end
    %If last bin is too small
    if expCounts13(end) < minCountNumber || expCounts23(end) < minCountNumber
        %Find what a first bin size that would have at least 5 counts
        %in both pdfs
        intEdgeMax = max(binEdges(2:end-1));
        if intEdgeMax >0
            intEdgeMaxEst = intEdgeMax*10;
        else
            intEdgeMaxEst = intEdgeMax/10;
        end
        maxVals2Try = linspace(min(binEdges(2:end-1)), intEdgeMaxEst, 100);
        expCountsMax1 = desiredSum.*(interp1(pdfVEC1.x, pdfVEC1.cdf_x, maxVals2Try));
        expCountsMax2 = desiredSum.*(interp1(pdfVEC2.x, pdfVEC2.cdf_x, maxVals2Try));
        maxEdgeInd = find(...
            expCountsMax1 < desiredSum - minCountNumber &...
            expCountsMax1 > desiredSum - maxCountNumber &...
            expCountsMax2 < desiredSum - minCountNumber &...
            expCountsMax2 > desiredSum - maxCountNumber, 1, 'last' );
        if sum(maxEdgeInd) <1
            figure()
            yyaxis left
            hold on
            plot(pdfVEC1.x, pdfVEC1.px, 'r')
            plot(pdfVEC2.x, pdfVEC2.px, 'b')
            error("No bin edges found that work for both distributions")
        end
        maxPossBinEdge = maxVals2Try(maxEdgeInd);   %Find a possible bin edge that would provide at least 5 counts to minimum bin
        maxBinEdge = max(binEdges3(binEdges3<maxPossBinEdge)); %Choose even spaced bin edge that is next biggest than minimum possible bin edge
        binEdges3 = [binEdges3(binEdges3 <= maxBinEdge), binEdges3(end)];
        [obsCountsDown,expCounts1, expCounts2] = calcBinCounts(data, binEdges3, desiredSum, pdfVEC1, pdfVEC2);
        binEdges = binEdges3;
        interiorEdges = binEdges(2:end-1);
        binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;
    end
end

%Perform chi2gof
[h1, p1, chiStat1] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts1, 'Emin', 5, 'NParams', pdfVEC1.numParams);
[h2, p2, chiStat2] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts2, 'Emin', 5, 'NParams', pdfVEC2.numParams);

%If chistat has changed binedges, use those
if numel(chiStat1.edges) ~= numel(binEdges)
    disp('a')
    binEdges = chiStat1.edges;
    binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;
    
    %Find observed and expected counts
    [obsCountsDown, expCounts1, expCounts2] = calcBinCounts(data, binEdges, desiredSum, pdfVEC1, pdfVEC2);
end

%Plot figure, with function, observed histogram, expected histogram, pvalue
%and chistat value
if fitS.dispChi2 == true

    % %Display Result of chi2gof to commandline
    % if h1 == 0
    %     disp("chi2gof: Accept H0; p = " +  num2str(p1, 3) + "; chistat = " + num2str(chiStat1.chi2stat, 3))
    % else
    %     disp("chi2gof: Reject H0 - p = " + num2str(p1, 3) + "; chistat = " + num2str(chiStat1.chi2stat, 3))
    % end

    %Plot figure
    figure;
    subplot(2,1,1)
    hold on
    %Plot histograms
    yyaxis("left")
    hObs = histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName", "Observed Counts");
    hExp = histogram('BinCounts', expCounts1, 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName","Expected Counts");
    if numel(data) == desiredSum
        ylabel("Counts")
    else
        ylabel("Downweighted Counts")
    end

    if fitS.enforceBinSizeLimits == false
        %Plot function
        yyaxis("right")
        lTested = plot(pdfVEC1.x, pdfVEC1.px, 'k', 'LineWidth', 2, "DisplayName", "Tested PDF"); %#ok<*NASGU>
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
    hObs = histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName", "Observed Counts");
    hExp = histogram('BinCounts', expCounts2, 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName","Expected Counts");
    if numel(data) == desiredSum
        ylabel("Counts")
    else
        ylabel("Downweighted Counts")
    end

    if fitS.enforceBinSizeLimits == false
        %Plot function
        yyaxis("right")
        lTested = plot(pdfVEC2.x, pdfVEC2.px, 'k', 'LineWidth', 2, "DisplayName", "Tested PDF"); %#ok<*NASGU>
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