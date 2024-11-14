function[] = outputMetadataAndSummaryFiguresComparison(subset1Log, subset2Log, inputDataTable1, inputDataTable2, subset1Color, subset2Color, subset1String, subset2String)
%% Take data of interest
%----- Find corenames, lats, longs, depths of cores included
dataTbl1 = inputDataTable1(subset1Log, :);

%----- Find corenames, lats, longs, depths of cores included
dataTbl2 = inputDataTable2(subset2Log,:);
%% Make Map

%----- Make map with locations denoted as red stars
figure;
worldmap('World')
setm(gca, 'mapprojection', 'robinson')
geoshow('landareas.shp', 'FaceColor','[0.7 0.7 0.7]', 'EdgeColor', '[0.7 0.7 0.7]')
load coastlines coastlon coastlat
plotm(coastlat, coastlon, 'Color', 'k')
hold on
p1 = plotm(dataTbl1.lats, dataTbl1.longs, 'Color', subset1Color, 'Marker', 'square', "LineStyle","none", "LineWidth", 1, "DisplayName", subset1String);
p2 = plotm(dataTbl2.lats, dataTbl2.longs, 'Color', subset2Color, 'Marker', 'square',"LineStyle","none","LineWidth", 1, "DisplayName", subset2String);
legend([p1 p2])

%% Plot Metadata Histograms

%Create lowSR and highSR bin counts for meanSR
mSR_topbinedge = ceil(max([dataTbl1.meanSR; dataTbl2.meanSR])/10)*10;
mSR_binwidth = 2;
mSR_binedges = 0:mSR_binwidth:mSR_topbinedge;
[mSRbottom, mSRtop] = createOverlaidHistcounts(dataTbl1.meanSR, dataTbl2.meanSR, mSR_binedges);

%Plot histograms of Mean SR, Depths, Resolution by Age, Resolution by Depth
figure;
subplot(4,1,1)
hold on
bar(mSR_binedges+(mSR_binwidth/2), mSRtop, "FaceColor", subset2Color, "BarWidth",1)
bar(mSR_binedges+(mSR_binwidth/2), mSRbottom, "FaceColor", subset1Color, "BarWidth",1)
xlabel('Mean SR (cm/kyr)')
ylabel('Counts')
xlim([0 mSR_topbinedge])
%title("Cores' Mean SR")

%Create lowSR and highSR bin counts for depth
depth_binwidth = 0.25;
depth_binedges = 0:depth_binwidth:6;
[depth_bottom, depth_top] = createOverlaidHistcounts(dataTbl1.depths./1000, dataTbl2.depths./1000, depth_binedges);

subplot(4,1,2)
hold on
bar(depth_binedges+(depth_binwidth/2), depth_top, "FaceColor", subset2Color, "BarWidth",1)
bar(depth_binedges+(depth_binwidth/2), depth_bottom, "FaceColor", subset1Color, "BarWidth",1)
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:0.5:6)
%title("Core Depths")

% %Create lowSR and highSR bin counts for d2coast
% d2coast_binwidth = 50;
% d2coast_binedges = 0:d2coast_binwidth:ceil(max(distance2coast)/d2coast_binwidth)*d2coast_binwidth;
% [d2coast_bottom, d2coast_top] = createOverlaidHistcounts(d2coast1, d2coast2, d2coast_binedges);
% 
% subplot(4,1,3)
% hold on
% bar(d2coast_binedges+(d2coast_binwidth/2), d2coast_top, "FaceColor", subset2Color, "BarWidth",1)
% bar(d2coast_binedges+(d2coast_binwidth/2), d2coast_bottom, "FaceColor", subset1Color, "BarWidth",1)
% ylabel("Counts")
% xlabel("Distance to coastline (km)")
% xlim([0 max(d2coast_binedges)+d2coast_binwidth])

%Create lowSR and highSR bin counts for MSI_byage0
MSI_byage_binwidth = 500;
MSI_byage_binedges = 0:500:ceil(max([dataTbl1.MSI_byage; dataTbl2.MSI_byage])/500)*500;
[MSI_byage_bottom, MSI_byage_top] = createOverlaidHistcounts(dataTbl1.MSI_byage, dataTbl2.MSI_byage, MSI_byage_binedges);

subplot(4,1,4)
hold on
bar(MSI_byage_binedges+(MSI_byage_binwidth/2), MSI_byage_top, "FaceColor", subset2Color, "BarWidth",1)
bar(MSI_byage_binedges+(MSI_byage_binwidth/2), MSI_byage_bottom, "FaceColor", subset1Color, "BarWidth",1)
xlabel('Mean Sampling Interval (Kyrs)')
ylabel('Counts')
xlim([0 max(MSI_byage_binedges)])
%title("Cores' Mean Sampling Interval By Age")

