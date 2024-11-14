%% Dist2Coast
%%% This function queries how far away a given lat/long pairing is from the
%%% nearest coastline using data from https://upwell.pfeg.noaa.gov/erddap/griddap/dist2coast_1deg.html

%Must add longs and lats to the workspace for this to work!

queryF = readmatrix("dist2coast1pt.csv","NumHeaderLines",2);

d2coast3 = NaN(length(lats), 1);

for i = 1:length(lats)
    latI = find(queryF(:,1) == round(lats(i),2));
    longI = find(queryF(:,2) == round(longs(i), 2));
    latlongI = intersect(latI, longI);
    d2coast3(i) = queryF(latlongI, 3);
end


d2coast4 = NaN(length(lats), 1);

for i = 1:length(lats)
    latI = find(queryF(:,1) == round(lats(i)+0.01,2));
    longI = find(queryF(:,2) == round(longs(i)+0.01, 2));
    latlongI = intersect(latI, longI);
    d2coast4(i) = queryF(latlongI, 3);
end
