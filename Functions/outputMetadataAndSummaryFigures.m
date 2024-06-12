function[] = outputMetadataAndSummaryFigures(subsetIND, cores, lats, longs, depths, meanSR, MSI_byage, MSI_bydepth, nSRcounts, sedimentlength, num14cpairs)
%% Take data of interest
%----- Find corenames, lats, longs, depths of cores included
core_inc = cores(subsetIND);
lat_inc = lats(subsetIND);
lon_inc = longs(subsetIND);
dep_inc = depths(subsetIND);
meanSR_inc = meanSR(subsetIND);
MSI_byage_inc = MSI_byage(subsetIND);
MSI_bydepth_inc = MSI_bydepth(subsetIND);
sedimentlength_inc = sedimentlength(subsetIND);
num14cpairs_inc = num14cpairs(subsetIND);

%% Make Map
%----- Make map with locations denoted as red stars
figure(23)
worldmap('World')
setm(gca, 'mapprojection', 'robinson')
geoshow('landareas.shp', 'FaceColor','[0.7 0.7 0.7]', 'EdgeColor', '[0.7 0.7 0.7]')
load coastlines
plotm(coastlat, coastlon, 'Color', 'k')
hold on
plotm(lat_inc, lon_inc, 'r*')
%% Plot Metadata Histograms

%Plot histograms of Mean SR, Depths, Resolution by Age, Resolution by Depth
figure(29)
subplot(4,1,1)
histogram(meanSR_inc, 0:2:90, 'FaceColor','k')
xlabel('Mean SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
title("Cores' Mean SR")

subplot(4,1,2)
histogram(dep_inc./1000, 0:0.25:6, 'FaceColor', 'k')
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:0.5:6)
title("Core Depths")

subplot(4,1,3)
histogram(MSI_byage_inc, 20, 'FaceColor','k')
xlabel('Sampling Frequency by age (kyr/date)')
ylabel('Counts')
title("Cores' Mean Sampling Interval By Age")

subplot(4,1,4)
histogram(MSI_bydepth_inc, 20, 'FaceColor','k')
%xlim([0 150])
xlabel('Sampling Frequency by depth (cm/date)')
ylabel('Counts')
%xticks(0:25:150)
title("Cores' Mean Sampling Interval By Depth")

%% Display summary information

%Input values from BIGMACS transition matrix
BIGMACSTM = [0.7215, 0.0940, 0.1846; 0.4328, 0.2687, 0.2985; 0.2670, 0.1041, 0.6290];

%Find out how much sediment the cores used constitute
lengthsed_core = nan(length(nSRcounts),1);
for i = 1:length(nSRcounts)
    if ~isempty(nSRcounts{i})
        lengthsed_core(i) = sum(nSRcounts{i}(2,:), 'omitmissing');
    end
end

%Display total number fo 14C pairs used and total length of sediment used
disp("The total number of 14C pairs used is")
disp(sum(num14cpairs_inc, 'omitmissing'))
disp("The total length of sediment used is")
disp(sum(sedimentlength_inc, 'omitmissing'))
disp(sum(lengthsed_core, 'omitmissing'))

