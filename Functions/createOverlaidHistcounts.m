function[countsBottomBar, countsTopBar] = createOverlaidHistcounts(data1, data2, binedges)
%This function takes two vectors of data observations, creates counts for
%them within a given set of bin edges, and then sums them to create a sum
%of the counts. This summed count can be plotted as a stacked bar chart,
%where the underlayed data is visible at the top of the bar chart.

%data1 = the data to be plotted as the first, lowest part of the histogram
%data2 = the data to be plotted as the second, highest part of the
%histogram

counts1 = histc(data1, binedges);
counts2 = histc(data2, binedges);
countssum = counts1 + counts2;

countsTopBar = countssum;
countsBottomBar = counts1;
end
