function [h, p, chiStat] = chi2_dataVStwopdfVECs(data, desiredSum, binN, pdfVs, fitS)
%Define some initial bins
binEdges    = linspace(min(data), max(data), binN+1);
interiorEdges = binEdges(2:end-1);
binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

%%% get cdfs of each pdf
for i = 1:size(pdfVs,1)
    cumsumpdfi = cumsum(pdfVs{i}.px);
    pdfVs{i}.cdf_x = cumsumpdfi./max(cumsumpdfi);
end

%Calculate observed and expected bin counts with given binEdges
[obsCountsDown, expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);

fitS.enforceBinSizeLimits = true;
if fitS.enforceBinSizeLimits
    minCountNumber = fitS.chi2MinCountNum;
    %minCountNumber = 5;
    maxCountNumber = desiredSum/3;

    %APPROACH #2 - Start with N equally spaced bins
    startN = binN+1;
    binEdges3    = linspace(min(data), max(data), startN);
    [~, expCountsN] = calcBinCounts(data, binEdges3, desiredSum, pdfVs);

    %Check whether the bin counts are too high, and increase number of bins
    newN = startN;

    while any(expCountsN(:,:) > maxCountNumber, "all")
        newN = newN+1;
        binEdges3    = linspace(min(data), max(data), newN);
        [~, expCountsN] = calcBinCounts(data, binEdges3, desiredSum, pdfVs);
    end

    %Start to decrease number of bins, until there are too few
    while sum(expCountsN(:,:) < maxCountNumber, "all") == 0
        %Try with less bins
        newN = newN-1;
        binEdges4    = linspace(min(data), max(data), newN);
        [~, expCountsN] = calcBinCounts(data, binEdges4, desiredSum, pdfVs);
    end
    
    %Add one more bin number (highest number of bins that keeps all bins
    %below maxCountNumber;
    if sum(expCountsN > maxCountNumber, "all") > 0
        newN = newN+1;
        binEdges3    = linspace(min(data), max(data), startN);
        [~, expCountsN] = calcBinCounts(data, binEdges3, desiredSum, pdfVs);
    end

    %Need to check which bins are too small now
    %If first bin is too small
    if any(expCountsN(:,1) < minCountNumber)
        %Find what a first bin size that would have at least 5 counts
        %in both pdfs
        intEdgeMin = min(binEdges(2:end-1));
        if intEdgeMin >0
            intEdgeMinEst = intEdgeMin/10;
        else
            intEdgeMinEst = intEdgeMin*10;
        end
        minVals2Try = linspace(intEdgeMinEst, max(binEdges(2:end-1)), 100);
        for i = 1:size(pdfVs,1)
            pdfi = pdfVs{i};
            expCountsMini = desiredSum.*(interp1(pdfi.x, pdfi.cdf_x, minVals2Try));
            expCountsMins(i, :) = expCountsMini;
        end
        expCountsMinsCompat = expCountsMins > minCountNumber & expCountsMins < maxCountNumber;
        minEdgeInd = find(all(expCountsMinsCompat,1), 1);
        if sum(minEdgeInd) <1
            figure()
            yyaxis left
            hold on
            for i = 1:size(pdfVs,1)
                pdfi = pdfVs{i};
                plot(pdfi.x, pdfi.px)
            end
            error("No bin edges found that work for both distributions")
        end
        minPossBinEdge = minVals2Try(minEdgeInd);   %Find a possible bin edge that would provide at least 5 counts to minimum bin
        minBinEdge = min(binEdges3(binEdges3>minPossBinEdge)); %Choose even spaced bin edge that is next biggest than minimum possible bin edge
        binEdges3 = [binEdges3(1), binEdges3(binEdges3 >=minBinEdge)];
        [~,expCountsN] = calcBinCounts(data, binEdges3, desiredSum, pdfVs);
    end
    %If last bin is too small
    if any(expCountsN(:,end) < minCountNumber)
        %Find what a first bin size that would have at least 5 counts
        %in both pdfs
        intEdgeMax = max(binEdges(2:end-1));
        if intEdgeMax >0
            intEdgeMaxEst = intEdgeMax*10;
        else
            intEdgeMaxEst = intEdgeMax/10;
        end
        maxVals2Try = linspace(min(binEdges(2:end-1)), intEdgeMaxEst, 100);
        for i = 1:size(pdfVs,1)
            pdfi = pdfVs{i};
            expCountsMaxi = desiredSum.*(interp1(pdfi.x, pdfi.cdf_x, maxVals2Try));
            expCountsMaxs(i, :) = expCountsMaxi;
        end
        expCountsMaxsCompat = expCountsMaxs > desiredSum - maxCountNumber & expCountsMaxs < desiredSum - minCountNumber;
        maxEdgeInd = find(all(expCountsMaxsCompat,1), 1, 'last');
        if sum(maxEdgeInd) <1
            figure()
            yyaxis left
            hold on
            for i = 1:size(pdfVs,1)
                pdfi = pdfVs{i};
                plot(pdfi.x, pdfi.px)
            end
            error("No bin edges found that work for both distributions")
        end
        maxPossBinEdge = maxVals2Try(maxEdgeInd);   %Find a possible bin edge that would provide at least 5 counts to maximum bin
        maxBinEdge = max(binEdges3(binEdges3<maxPossBinEdge)); %Choose even spaced bin edge that is next smallest than maximum possible bin edge
        binEdges3 = [binEdges3(binEdges3 <= maxBinEdge), binEdges3(end)];
        [obsCountsDown,expCounts] = calcBinCounts(data, binEdges3, desiredSum, pdfVs);
        binEdges = binEdges3;
        interiorEdges = binEdges(2:end-1);
        binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;
    end
end

%Perform chi2gof
for i = 1:size(pdfVs, 1)
[h(i), p(i), chiStat{i}] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts(i,:), 'Emin', 5, 'NParams', pdfVs{i}.numParams);
end

%If chistat has changed binedges, use those
if numel(chiStat{1}.edges) ~= numel(binEdges)
    binEdges = chiStat{1}.edges;
    binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;
    
    %Find observed and expected counts
    [obsCountsDown, expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);
end

%Plot figure, with function, observed histogram, expected histogram, pvalue
%and chistat value
if fitS.dispChi2 == true

    %loop for figure
    figure;
    for i = 1:size(pdfVs, 1)
        subplot(size(pdfVs, 1), 1, i)
        hold on
        %Plot histograms
        %yyaxis("left")
        hObs = histogram('BinCounts',obsCountsDown, 'BinEdges', binEdges, 'FaceColor', 'k', "DisplayName", "Observed Counts");
        hExp = histogram('BinCounts', expCounts(i,:), 'BinEdges', binEdges, 'FaceColor', 'r', "DisplayName","Expected Counts");
        if numel(data) == desiredSum
            ylabel("Counts")
        else
            ylabel("Downweighted Counts")
        end

        if fitS.enforceBinSizeLimits == false
            %Plot function
            yyaxis("right")
            lTested = plot(pdfVs{i}.x, pdfVs{i}.px, 'k', 'LineWidth', 2, "DisplayName", "Tested PDF"); %#ok<*NASGU>
            ylabel("PDF")
        end

        %Add pval and chi2stat to legend
        lpval   = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "p = " + num2str(p(i), 3));
        lchi2stat = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "chi2stat = "  + num2str(chiStat{i}.chi2stat, 3));
        %xlabel("log nSR") Note, it's not always log nSR...
        legend()
        xlim([-3 3])
        title(pdfVs{i}.pdfName)
    end
end
end