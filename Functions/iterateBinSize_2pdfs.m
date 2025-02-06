function[interiorEdges] = iterateBinSize_2pdfs(expCounts1, pdfVEC1, expCounts2, pdfVEC2, interiorEdges, desiredSum, minCountNumber, maxCountNumber)
%If the bins have less than 5 counts or more than 12, recalculate bin
%counts until all bins fit this criteria
smallbins1 = expCounts1 < minCountNumber;
bigbins1 = expCounts1 > maxCountNumber;
binSizeTester1 = smallbins1| bigbins1;

smallbins2 = expCounts2 < minCountNumber;
bigbins2 = expCounts2 > maxCountNumber;
binSizeTester2 = smallbins2| bigbins2;

smallbins = smallbins1 | smallbins2;
bigbins = bigbins1 | bigbins2;
binSizeTester = binSizeTester1 | binSizeTester2;

counter = 0;

while(sum(binSizeTester) ~=0)
    counter = counter+1;

    %If the first bin is too big or small, find a good size
    if binSizeTester(1) == 1
        %Find a good lowest bin boundary
        intEdgeMin = min(interiorEdges);
        if intEdgeMin >0
            intEdgeMinEst = intEdgeMin/10;
        else
            intEdgeMinEst = intEdgeMin*10;
        end
        minVals2Try = linspace(intEdgeMinEst, max(interiorEdges), 100);
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
        minimumBinEdge = minVals2Try(minEdgeInd);
        interiorEdges = sort([minimumBinEdge, interiorEdges]); %#ok<AGROW>
        interiorEdges = interiorEdges(interiorEdges >= minimumBinEdge);
        expProbs1 = diff([0, (interp1(pdfVEC1.x, pdfVEC1.cdf_x, interiorEdges)), 1]);
        expCounts1 = expProbs1.*desiredSum;
        expProbs2 = diff([0, (interp1(pdfVEC2.x, pdfVEC2.cdf_x, interiorEdges)), 1]);
        expCounts2 = expProbs2.*desiredSum;
        %Test again
        smallbins1 = expCounts1 < minCountNumber;
        bigbins1 = expCounts1 > maxCountNumber;
        binSizeTester1 = smallbins1| bigbins1;

        smallbins2 = expCounts2 < minCountNumber;
        bigbins2 = expCounts2 > maxCountNumber;
        binSizeTester2 = smallbins2| bigbins2;

        smallbins = smallbins1 | smallbins2;
        bigbins = bigbins1 | bigbins2;
        binSizeTester = binSizeTester1 | binSizeTester2;
    end
    %If the last bin is too big, find a good size
    if binSizeTester(end) == 1
        %Find a good lowest bin boundary
        intEdgeMax = max(interiorEdges);
        if intEdgeMax >0
            intEdgeMaxEst = intEdgeMax*10;
        else
            intEdgeMaxEst = intEdgeMax/10;
        end
        maxVals2Try = linspace(min(interiorEdges),intEdgeMaxEst, 100);
        expCountsmax1 = desiredSum.*(interp1(pdfVEC1.x, pdfVEC1.cdf_x, maxVals2Try));
        expCountsmax2 = desiredSum.*(interp1(pdfVEC2.x, pdfVEC2.cdf_x, maxVals2Try));
        try
            maxEdgeInd = find(...
                expCountsmax1 < desiredSum - minCountNumber &...
                expCountsmax1 > desiredSum - maxCountNumber &...
                expCountsmax2 < desiredSum - minCountNumber &...
                expCountsmax2 > desiredSum - maxCountNumber, 1 );
        catch
            error("No bin edges found that work for both distributions")
        end
        maximumBinEdge = maxVals2Try(maxEdgeInd);
        interiorEdges = sort([interiorEdges, maximumBinEdge]); %#ok<AGROW>
        expProbs1 = diff([0, (interp1(pdfVEC1.x, pdfVEC1.cdf_x, interiorEdges)), 1]);
        expCounts1 = expProbs1.*desiredSum;
        expProbs2 = diff([0, (interp1(pdfVEC2.x, pdfVEC2.cdf_x, interiorEdges)), 1]);
        expCounts2 = expProbs2.*desiredSum;
        %Test again
        smallbins1 = expCounts1 < minCountNumber;
        bigbins1 = expCounts1 > maxCountNumber;
        binSizeTester1 = smallbins1| bigbins1;

        smallbins2 = expCounts2 < minCountNumber;
        bigbins2 = expCounts2 > maxCountNumber;
        binSizeTester2 = smallbins2| bigbins2;

        smallbins = smallbins1 | smallbins2;
        bigbins = bigbins1 | bigbins2;
        binSizeTester = binSizeTester1 | binSizeTester2;
    end
    if sum(smallbins) ~=0
        %Remove all bin edges that create bins that are too small
        if sum(smallbins) == numel(smallbins)
            %If all bin edges would be removed, then remove every other bin
            %edge
            interiorEdges = interiorEdges(logical(accumarray([1:2:numel(smallbins)]',1)')); %#ok<NBRAK1>
        end
        interiorEdges = interiorEdges(~smallbins(2:end));
        expProbs1 = diff([0, (interp1(pdfVEC1.x, pdfVEC1.cdf_x, interiorEdges)), 1]);
        expCounts1 = expProbs1.*desiredSum;
        expProbs2 = diff([0, (interp1(pdfVEC2.x, pdfVEC2.cdf_x, interiorEdges)), 1]);
        expCounts2 = expProbs2.*desiredSum;
        %Test again
        smallbins1 = expCounts1 < minCountNumber;
        bigbins1 = expCounts1 > maxCountNumber;
        binSizeTester1 = smallbins1| bigbins1;

        smallbins2 = expCounts2 < minCountNumber;
        bigbins2 = expCounts2 > maxCountNumber;
        binSizeTester2 = smallbins2| bigbins2;

        smallbins = smallbins1 | smallbins2;
        bigbins = bigbins1 | bigbins2;
        binSizeTester = binSizeTester1 | binSizeTester2;
    end
    if sum(bigbins) ~=0
        %split big bins into two
        bigBinMinEdge = logical([bigbins(2:end-1), 0]);
        bigBinMaxEdge = logical([0,bigbins(2:end-1)]);
        bigBinsSplit = (interiorEdges(bigBinMinEdge)+interiorEdges(bigBinMaxEdge))./2;
        %Add the remaining edges and halfway edges
        interiorEdges = sort([interiorEdges,bigBinsSplit]);
        %Calculate expected counts now
        expProbs1 = diff([0, (interp1(pdfVEC1.x, pdfVEC1.cdf_x, interiorEdges)), 1]);
        expCounts1 = expProbs1.*desiredSum;
        expProbs2 = diff([0, (interp1(pdfVEC2.x, pdfVEC2.cdf_x, interiorEdges)), 1]);
        expCounts2 = expProbs2.*desiredSum;
        %Test again
        smallbins1 = expCounts1 < minCountNumber;
        bigbins1 = expCounts1 > maxCountNumber;
        binSizeTester1 = smallbins1| bigbins1;

        smallbins2 = expCounts2 < minCountNumber;
        bigbins2 = expCounts2 > maxCountNumber;
        binSizeTester2 = smallbins2| bigbins2;

        smallbins = smallbins1 | smallbins2;
        bigbins = bigbins1 | bigbins2;
        binSizeTester = binSizeTester1 | binSizeTester2;
    end

    if counter >= 1001
        error("Bin Iteration taking too many iterations")
    end
end
