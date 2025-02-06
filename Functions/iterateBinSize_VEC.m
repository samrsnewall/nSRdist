function[interiorEdges] = iterateBinSize_VEC(expCounts, pdfVEC, interiorEdges, desiredSum, minCountNumber, maxCountNumber)
%If the bins have less than 5 counts or more than 12, recalculate bin
%counts until all bins fit this criteria
smallbins = expCounts < minCountNumber;
bigbins = expCounts > maxCountNumber;
binSizeTester = smallbins| bigbins;

while(sum(binSizeTester) ~=0)
    %If the first bin is too big, find a good size
    if binSizeTester(1) == 1
        %Find a good lowest bin boundary
        intEdgeMin = min(interiorEdges);
        if intEdgeMin >0
            intEdgeMinEst = intEdgeMin/50;
        else
            intEdgeMinEst = intEdgeMin*50;
        end
        minVals2Try = linspace(intEdgeMinEst, max(interiorEdges), 500);
        expCountsMin = desiredSum.*(interp1(pdfVEC.x, pdfVEC.cdf_x, minVals2Try));
        minEdgeInd = find(expCountsMin > minCountNumber & expCountsMin <maxCountNumber, 1 );
        minimumBinEdge = minVals2Try(minEdgeInd);
        interiorEdges = sort([minimumBinEdge, interiorEdges]); %#ok<AGROW>
        interiorEdges = interiorEdges(interiorEdges >= minimumBinEdge);
        expProbs = diff([0, (interp1(pdfVEC.x, pdfVEC.cdf_x, interiorEdges)), 1]);
        expCounts = expProbs.*desiredSum;
        %Test again
        smallbins = expCounts < minCountNumber;
        bigbins = expCounts > maxCountNumber;
        binSizeTester = smallbins| bigbins;
    end
    %If the last bin is too big, find a good size
    if binSizeTester(end) == 1
        %Find a good lowest bin boundary
        intEdgeMax = max(interiorEdges);
        if intEdgeMax >0
            intEdgeMaxEst = intEdgeMax*50;
        else
            intEdgeMaxEst = intEdgeMax/50;
        end
        maxVals2Try = linspace(min(interiorEdges),intEdgeMaxEst, 500);
        expCountsmax = desiredSum.*(interp1(pdfVEC.x, pdfVEC.cdf_x, maxVals2Try));
        maxEdgeInd = find(expCountsmax < desiredSum-minCountNumber & expCountsmax > desiredSum - maxCountNumber, 1 );
        maximumBinEdge = maxVals2Try(maxEdgeInd);
        interiorEdges = sort([interiorEdges, maximumBinEdge]); %#ok<AGROW>
        expProbs = diff([0, (interp1(pdfVEC.x, pdfVEC.cdf_x, interiorEdges)), 1]);
        expCounts = expProbs.*desiredSum;
        %Test again
        smallbins = expCounts < minCountNumber;
        bigbins = expCounts > maxCountNumber;
        binSizeTester = smallbins| bigbins;
    end
    if sum(smallbins) ~=0
        %Remove all bin edges that create bins that are too small
        if sum(smallbins) == numel(smallbins)
            %If all bin edges would be removed, then remove every other bin
            %edge
            interiorEdges = interiorEdges(logical(accumarray([1:2:numel(smallbins)]',1)')); %#ok<NBRAK1>
        end
        interiorEdges = interiorEdges(~smallbins(2:end));
        expProbs = diff([0, (interp1(pdfVEC.x, pdfVEC.cdf_x, interiorEdges)), 1]);
        expCounts = expProbs.*desiredSum;
        %Test again
        smallbins = expCounts < minCountNumber;
        bigbins = expCounts > maxCountNumber;
        binSizeTester = smallbins | bigbins;
    end
    if sum(bigbins) ~=0
        %split big bins into two
        bigBinMinEdge = logical([bigbins(2:end-1), 0]);
        bigBinMaxEdge = logical([0,bigbins(2:end-1)]);
        bigBinsSplit = (interiorEdges(bigBinMinEdge)+interiorEdges(bigBinMaxEdge))./2;
        %Add the remaining edges and halfway edges
        interiorEdges = sort([interiorEdges,bigBinsSplit]);
        %Calculate expected counts now
        expProbs = diff([0, (interp1(pdfVEC.x, pdfVEC.cdf_x, interiorEdges)), 1]);
        expCounts = expProbs.*desiredSum;
        %Test again
        smallbins = expCounts < minCountNumber;
        bigbins = expCounts > maxCountNumber;
        binSizeTester = smallbins| bigbins;
    end
end


