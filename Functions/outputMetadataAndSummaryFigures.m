function [] = outputMetadataAndSummaryFigures(subsetLog, inputDataT)
% outputMetadataAndSummaryFigures  Produce summary map and metadata
%                                   histograms for a subset of cores.
%
% Generates two figures and prints summary statistics to the command
% window for the cores identified by subsetLog in inputDataT:
%
%   Figure 1 — Robinson-projection world map with core locations marked.
%   Figure 2 — Four-panel histogram of mean SR, water depth, and mean
%               sampling intervals (by age and by depth).
%
% Also prints to the command window: total number of 14C date pairs used
% and total length of sediment analysed.
%
% INPUTS
%   subsetLog  - (logical vector) Rows of inputDataT to include
%   inputDataT - (table) Core metadata table, as assembled by calcData.
%                Required columns: lats, longs, meanSR, depths,
%                MSI_byage, MSI_bydepth, nSRcounts, num14cpairs,
%                sedimentlength
%
% See also: calcData

%% Select cores of interest
dataTbl  = inputDataT(subsetLog, :);
numCores = height(dataTbl);

%% Figure 1: World map of core locations
figure;
worldmap('World')
setm(gca, 'mapprojection', 'robinson')
geoshow('landareas.shp', 'FaceColor', '[0.7 0.7 0.7]', 'EdgeColor', '[0.7 0.7 0.7]')
load coastlines coastlat coastlon
plotm(coastlat, coastlon, 'Color', 'k')
hold on
plotm(dataTbl.lats, dataTbl.longs, 'rs')

%% Figure 2: Metadata histograms
figure;

subplot(4,1,1)
histogram(dataTbl.meanSR, 0:2:90, 'FaceColor', 'k')
xlabel('Mean SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
title("Cores' Mean SR")

subplot(4,1,2)
histogram(dataTbl.depths ./ 1000, 0:0.25:6, 'FaceColor', 'k')
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:0.5:6)
title("Core Depths")

subplot(4,1,3)
histogram(dataTbl.MSI_byage, 20, 'FaceColor', 'k')
xlabel('Sampling Frequency by age (kyr/date)')
ylabel('Counts')
title("Cores' Mean Sampling Interval By Age")

subplot(4,1,4)
histogram(dataTbl.MSI_bydepth, 20, 'FaceColor', 'k')
xlabel('Sampling Frequency by depth (cm/date)')
ylabel('Counts')
title("Cores' Mean Sampling Interval By Depth")

%% Print summary statistics
% Compute total sediment length from nSRcounts Row 2 (dep_diffs), excluding NaN header columns.
lengthsed_core = nan(numCores, 1);
for i = 1:numCores
    if ~isempty(dataTbl.nSRcounts{i})
        nSRi = dataTbl.nSRcounts{i};
        NaN_logi = ~isnan(nSRi(1,:));
        lengthsed_core(i) = sum(nSRi(2, NaN_logi), 'omitmissing');
    end
end

disp("The total number of 14C pairs used is")
disp(sum(dataTbl.num14cpairs, 'omitmissing'))
disp("The total length of sediment used is")
disp(sum(lengthsed_core, 'omitmissing'))
disp(sum(dataTbl.sedimentlength, 'omitmissing'))

end
