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
    [binEdges, expCounts, obsCountsDown, binCenters] = iterateBinSizeNVec(data, pdfVs, binN, desiredSum, fitS);
end

%Perform chi2gof
for i = 1:size(pdfVs, 1)
    [h(i), p(i), chiStat{i}] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts(i,:), 'Emin', 5, 'NParams', pdfVs{i}.numParams);
end

%If chistat has changed binedges, use those
if numel(chiStat{1}.edges) ~= numel(binEdges)
    disp('chi2gof function modified bins')
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
        %Set x limits
        binEdgesDiff = (binEdges(3) - binEdges(2));
        xlim([binEdges(2)-binEdgesDiff binEdges(end-1)+binEdgesDiff])

        %Plot observed counts
        hObs = histogram('BinCounts',obsCountsDown(2:end-1), 'BinEdges', binEdges(2:end-1), 'FaceColor', 'k', "DisplayName", "Observed Counts");
        plot([-100 binEdges(2)], [obsCountsDown(1) obsCountsDown(1)], 'Color', [0.3 0.3 0.3], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")
        plot([100 binEdges(end-1)], [obsCountsDown(end) obsCountsDown(end)], 'Color', [0.3 0.3 0.3], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")

        %Plot expected counts
        hExp = histogram('BinCounts', expCounts(i,2:end-1), 'BinEdges', binEdges(2:end-1), 'FaceColor', 'r', "DisplayName","Expected Counts");
        plot([-100 binEdges(2)], [expCounts(1) expCounts(1)], 'Color', [0.9 0 0], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")
        plot([100 binEdges(end-1)], [expCounts(end) expCounts(end)], 'Color', [0.9 0 0], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")
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
        lchi2stat = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "chi2stat = "  + num2str(chiStat{i}.chi2stat, 3));
        lpval   = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "p = " + num2str(p(i), 3) + "; h = " + num2str(h(i), 1));
        %lhnought= plot(nan, nan, 'LineStyle', 'none', "Marker", "none", 'DisplayName', "h = " + num2str(h(i), 1));
        
        legend()
        % xlim([-3 3])
        title(pdfVs{i}.pdfName)
    end
    xlabel("log(nSR)") %Note, it's not always log nSR...
end
end