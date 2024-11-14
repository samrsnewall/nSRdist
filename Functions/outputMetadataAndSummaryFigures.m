function[] = outputMetadataAndSummaryFigures(subsetLog, inputDataT)
%% Take data of interest
%----- Find corenames, lats, longs, depths of cores included
dataTbl = inputDataT(subsetLog, :);

%% Make Map
%----- Make map with locations denoted as red squares
figure;
worldmap('World')
setm(gca, 'mapprojection', 'robinson')
geoshow('landareas.shp', 'FaceColor','[0.7 0.7 0.7]', 'EdgeColor', '[0.7 0.7 0.7]')
load coastlines coastlat coastlon
plotm(coastlat, coastlon, 'Color', 'k')
hold on
plotm(dataTbl.lats, dataTbl.longs, 'rs')
%% Plot Metadata Histograms

%Plot histograms of Mean SR, Depths, Resolution by Age, Resolution by Depth
figure;
subplot(4,1,1)
histogram(dataTbl.meanSR, 0:2:90, 'FaceColor','k')
xlabel('Mean SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
title("Cores' Mean SR")

subplot(4,1,2)
histogram(dataTbl.depths./1000, 0:0.25:6, 'FaceColor', 'k')
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:0.5:6)
title("Core Depths")

subplot(4,1,3)
histogram(dataTbl.MSI_byage, 20, 'FaceColor','k')
xlabel('Sampling Frequency by age (kyr/date)')
ylabel('Counts')
title("Cores' Mean Sampling Interval By Age")

subplot(4,1,4)
histogram(dataTbl.MSI_bydepth, 20, 'FaceColor','k')
%xlim([0 150])
xlabel('Sampling Frequency by depth (cm/date)')
ylabel('Counts')
%xticks(0:25:150)
title("Cores' Mean Sampling Interval By Depth")

%% Display summary information

%Input values from BIGMACS transition matrix
BIGMACSTM = [0.7215, 0.0940, 0.1846; 0.4328, 0.2687, 0.2985; 0.2670, 0.1041, 0.6290];

%Find out how much sediment the cores used constitute
numCores = length(dataTbl.lats);
lengthsed_core = nan(numCores,1);
for i = 1:numCores
    if ~isempty(dataTbl.nSRcounts{i})
        lengthsed_core(i) = sum(dataTbl.nSRcounts{i}(2,:), 'omitmissing');
    end
end

%Display total number fo 14C pairs used and total length of sediment used
disp("The total number of 14C pairs used is")
disp(sum(dataTbl.num14cpairs, 'omitmissing'))
disp("The total length of sediment used is")
disp(sum(lengthsed_core, 'omitmissing'))
disp(sum(dataTbl.sedimentlength, 'omitmissing'))

