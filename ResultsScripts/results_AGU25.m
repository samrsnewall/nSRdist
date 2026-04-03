%% Add necessary paths
addpath("../Functions/")

%% Set up histogram settings
%Create binedges
maxbinedge   = 20;
bw           = 0.1; %bin width
binEdges     =  0:bw:maxbinedge; %nSR
invBinEdges  =  0:bw:maxbinedge; %invnSR
logBinEdges  = -5:bw:5;          %lognSR
logBinCenters = (logBinEdges(1:end-1)+logBinEdges(2:end)).*0.5;

%% Fit Mix Log Norm to BIGMACS data for comparison reasons
%Load BIGMACS files
BMlognorm = readtable("../BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
BM.nSR.x = BMlognorm.Var1';
BM.nSR.px = BMlognorm.Var2;
BM.hist = readmatrix("../BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BM.hist = BM.hist(:,4);
BM.TM = readmatrix("../BIGMACSdata/transition_parameter.txt");
[MLN_BIGMACS, ~, BM.gmfit] = fitMixLogNorm(BM.hist, BM.nSR.x, 2, 3);
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);
BMpdf.numParams = 6;
BMpdf.x = BM.lnSR.x;
BMpdf.px = BM.lnSR.px;

%% Load data and fits
%dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Jun2_fitJun24_depthweight.mat");
dA = load("../Results/dataT_All1_RLGtrue_DS0p05_Dec9_fitDec9_depthweight_400R.mat");

%% Take data of interest
%----- Find corenames, lats, longs, depths of cores included
dataTbl = dA.d.dataT(dA.d.S1.chooseLog, :);
%----- Make map with locations denoted as red squares
figure;
subplot(2,2,[1 2])
worldmap('World')
setm(gca, 'mapprojection', 'robinson')
geoshow('landareas.shp', 'FaceColor','[0.7 0.7 0.7]', 'EdgeColor', '[0.7 0.7 0.7]')
load coastlines coastlat coastlon
plotm(coastlat, coastlon, 'Color', 'k')
hold on
plotm(dataTbl.lats, dataTbl.longs, 'rs')
%Plot histograms of Mean SR, Depths, Resolution by Age, Resolution by Depth

subplot(2,2,3)
histogram(dataTbl.depths./1000, 0:0.25:6, 'FaceColor', 'k')
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:0.5:6)
title("Core Depths")

subplot(2,2,4)
histogram(dataTbl.meanSR, 0:2:90, 'FaceColor','k')
xlabel('Mean SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
title("Cores' Mean SR")

%% Plot age modes 
%Find out which cores passed filtering
usedLog = ~isnan(dA.d.dataT.meanSR);

plotAgeModes(dA.d.S1.chooseLog,dA.d.S1.chooseLog, dA.d.dataT.ageModes, dA.d.dataT.cores)

%% %Create table of useful information
dA.d.S1.BMode.weightedC = dA.d.S1.BMode.weightedC;
dA.d.S1.BMode.MLN.chiStats = dA.d.S1.BMode.MLN.chiStats;
dA.d.S1.BMode.invGam.chiStats = dA.d.S1.BMode.invGam.chiStats;
dStrus = {dA.d.S1.BMode, dA.d.S1.BMedian, dA.d.S1.BSampIR, dA.d.S1.New0IR, dA.d.S1.New500IR, dA.d.S1.New1000IR, dA.d.S1.New1500IR};
dStrusStrings = ["BMode", "BMedian", "BSamp", "RSR0", "RSR500", "RSR1000", "RSR1500"];
numRunsEva = length(dA.d.S1.BSampIR.MLN.chiStats.h);

MeanAgePairsT  = NaN(length(dStrus),1);
MeanSedLength  = NaN(length(dStrus),1);
MeanSedTime    = NaN(length(dStrus),1);
nSR_median     = NaN(length(dStrus),1);
nSR_min     = NaN(length(dStrus),1);
nSR_max     = NaN(length(dStrus),1);
nSR_95lo     = NaN(length(dStrus),1);
nSR_95hi     = NaN(length(dStrus),1);
MLNacc      = NaN(length(dStrus),1);
LNacc      = NaN(length(dStrus),1);
InvGamAcc   = NaN(length(dStrus),1);
GamAcc   = NaN(length(dStrus),1);
BMacc    = NaN(length(dStrus),1);

for i = 1:length(dStrus)
    dStru = dStrus{i};
    if isa(dStru.weightedC, 'cell')
        allWC = cell2mat(dStru.weightedC);
            MeanSedLength(i,1) = mean(dStru.sedLength);
    MeanSedTime(i,1) = mean(dStru.sedTimeSpan);
    else
        allWC = dStru.weightedC;
        MeanSedLength(i,1) = mean(dStru.sedLength);
        MeanSedTime(i,1)   = mean(dStru.sedTimeSpan);
    end
    numRunsEva = length(dStru.MLN.chiStats.h);
    
    MeanAgePairsT(i,1) = mean(dStru.numCpairs);
    allWCsorted = sort(allWC);
    numWC = length(allWC);
    
    nSR_median(i,1) = median(allWC);
    nSR_min(i,1) = min(allWC);
    nSR_max(i,1) = max(allWC);
    nSR_95lo(i,1) = allWCsorted(ceil(0.025*numWC));
    nSR_95hi(i,1) = allWCsorted(ceil(0.975*numWC));
    MLNacc(i,1) = sum(dStru.MLN.chiStats.h == 0)./numRunsEva;
    InvGamAcc(i,1) = sum(dStru.invGam.chiStats.h == 0)./numRunsEva;
    GamAcc(i,1) = sum(dStru.Gam.chiStats.h == 0)./numRunsEva;
    LNacc(i,1) = sum(dStru.LN.chiStats.h == 0)./numRunsEva;
    BMacc(i,1) = sum(dStru.chiStatTvsBM.h == 0)./numRunsEva;  
end
dA.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median,...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi, LNacc, MLNacc, GamAcc, InvGamAcc,...
    BMacc, 'RowNames', dStrusStrings);

%% BIGMACS Histogram and fit
figure;
histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 1500*BM.lnSR.px, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 2000])
ylabel("cm")
xlabel("log(nSR)")
title("BIGMACS - Lee et al., 2023")

%% Plot of BIGMACS and BMedian histograms with BIGMACS distribution superimposed
nsubs = 2;
figure
subplot(nsubs,1,2-1)
histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 1500*BM.lnSR.px, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BIGMACS: var = " + num2str(var(log(BM.hist)), 3))
subplot(nsubs,1,3-1)
histogram(log(dA.d.S1.BMedian.weightedC{1}),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian: var = " + num2str(var(log(dA.d.S1.BMedian.weightedC{1})), 3))

%% Plot of above with normalized histogram

nsubs = 2;
fA = figure;
subplot(nsubs,1,2-1)
BM.logHC = histcounts(log(BM.hist), logBinEdges);
normF = 1./sum((BM.logHC.*uniquetol(diff(logBinEdges), 1e-5)));
hold on
histogram('BinCounts', BM.logHC, 'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BIGMACS")
subplot(nsubs,1,3-1)
BMedian_logHC = histcounts(log(dA.d.S1.BMedian.weightedC{1}), logBinEdges);
normF = 1./sum((BMedian_logHC.*uniquetol(diff(logBinEdges), 1e-5)));
histogram('BinCounts', BMedian_logHC,'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1);
hold on
plot(BM.lnSR.x, BM.lnSR.px*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian")
xlabel("log(nSR)")

%saveas(fA, "../Figures/BMvsBMedian_normalizedPDF.jpg")

%% 
nsubs = 2;
figure
subplot(nsubs,1,2-1)
histogram(log(dA.d.S1.BMode.weightedC{1}),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMode")
subplot(nsubs,1,3-1)
histogram(log(dA.d.S1.BMedian.weightedC{1}),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian; var =" + num2str(var(log(dA.d.S1.BMedian.weightedC{1}))))

%% All Bchron and RSR500 samplings, histogram with 68th percentile bars
nsubs = 3;
numruns = length(dA.d.S1.BSampIR.OneRunDatas);
 figure;
% subplot(nsubs,1,1)
% histogram(log(BM.hist),'BinEdges', logBinEdges, 'FaceColor', 'b', 'FaceAlpha', 0.1)
% hold on
% plot(BM.lnSR.x, 1500*BM.lnSR.px, 'LineWidth', 2)
% box on
% xlim([-2.5 2.5])
% ylim([0 3500])
% ylabel("cm")
% title("BIGMACS - Lee et al., 2023")
subplot(nsubs,1,2-1)
histogram(log(dA.d.S1.BMedian.weightedC{1}),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
box on
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BMedian [500-4000yr]")
subplot(nsubs,1,3-1)
BSamp_hists = sort(dA.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', BSamp_hists(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(nsubs,1,4-1)
MCSamp_hists = sort(dA.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', MCSamp_hists(numruns*0.5, :), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
hold on
plot(BM.lnSR.x, 2700*BM.lnSR.px, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
ylim([0 3500])
ylabel("cm")
xlabel("log(nSR)")
title("MCSamp")

%% Compare the variance of each approach
commonxlim = [0 5];

figure;
subplot(4,1,1)
xline(var(BM.hist))
xlim(commonxlim)
subplot(4,1,2)
xline(var(dA.d.S1.BMedian.weightedC{1}))
xlim(commonxlim)
subplot(4,1,3)
histogram(cellfun(@var, dA.d.S1.BSampIR.weightedC), 20)
xlim(commonxlim)
subplot(4,1,4)
histogram(cellfun(@var, dA.d.S1.New500IR.weightedC), 20)
xlim(commonxlim)

%% Estimate the variance of the median histograms

counts1 = MCSamp_hists(numruns*0.5, :)';      % ensure column vector
edges   = logBinEdges(:);

% Precompute midpoints and widths (shared by all histograms)
mid    = 0.5 * (edges(1:end-1) + edges(2:end));
widths = diff(edges);

N1  = sum(counts1);
mu1 = sum(counts1 .* mid) / N1;

var_MCSamp = sum(counts1 .* ( (mid - mu1).^2 + widths.^2 / 12 )) / N1;
std1 = sqrt(var1);

counts2 = BSamp_hists(numruns*0.5, :)';       % ensure column
N2  = sum(counts2);
mu2 = sum(counts2 .* mid) / N2;

var_BSamp = sum(counts2 .* ( (mid - mu2).^2 + widths.^2 / 12 )) / N2;
std2 = sqrt(var2);


%% Plot the many different fits

figure;
hold on
plot(dA.d.S1.New500IR.MLN.lnSR.x, dA.d.S1.New500IR.MLN.lnSR.px, 'k')
plot(dA.d.S1.New500IR.invGam.lnSR.x, dA.d.S1.New500IR.invGam.lnSR.px, 'r')


%% Load data of just A7 core
d1C = load("../Results/dataT_RLGtrue_R200M20_Apr23_A7_fitDec8_depthweight_1000R.mat");

figure;
subplot(3,1,1)
histogram(log(d1C.d.S1.BMedian.weightedC{1}),'BinEdges', logBinEdges, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'DisplayName', 'Replicated nSR counts');
box on
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
title("BMedian [500-4000yr]")
subplot(3,1,2)
A = sort(d1C.d.S1.BSampIR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
title("BSamp [500-4000yr]")
subplot(3,1,3)
A = sort(d1C.d.S1.New500IR.lnSRHistCounts, 1);
hold on
box on
histogram('BinCounts', mean(A), 'BinEdges', logBinEdges, 'FaceAlpha', 0.1)
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.16, :)),(A(numruns*0.84, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '68% bounds')
%errorbar(logBinCenters, A(numruns*0.5, :),(A(numruns*0.5, :)-A(numruns*0.025, :)),(A(numruns*0.975, :)-A(numruns*0.5, :)), 'LineStyle', 'none', 'DisplayName', '95% bounds')
xlim([-2.5 2.5])
%ylim([0 100])
ylabel("cm")
xlabel("log(nSR)")
title("MCSamp")


%% Plot BMedian for another core
%Get nSR vector
median_nSRs = dA.d.dataT.bchronMedian;

%Get median ages and depths
median_ages = cumsum(median_nSRs{31}(3,:))/1000;
median_depths = cumsum(median_nSRs{31}(2,:));

%Construct repeated vectors for plotting with NaNs
median_depths2 = repelem(median_depths, 2);
median_ages2 = repelem(median_ages, 2);

figure;
subplot(3,1,1)
hold on
plot(median_depths, median_ages, 'ko', 'LineStyle', 'none','DisplayName', "Median Age", 'LineWidth', 1)
plot(median_depths2, median_ages2, '-r', 'Marker', '.')
xlabel("Depth (cm)")
ylabel("Age (kyr)")
xlims = [220 340];
ylims = [3 7.5];
xlim(xlims)
ylim(ylims)

%Find indices of nSRs that have dt < 500;
filtOutLogi = median_nSRs{31}(3,2:end) < 500;
filtOutInd = find(filtOutLogi);
filtOutInd2 = filtOutInd*2 - 1;

%Put NaN in repeated vectors so they nSRs of dt <500 do not plot
median_ages2NaN = median_ages2;
%median_ages2NaN(filtOutInd2) = NaN;
median_ages2NaN((median_depths2 == 252.5 | median_depths2 == 260.5 | median_depths2 == 334.5)) = NaN;

subplot(3,1,2)
hold on
plot(median_depths, median_ages, 'ko', 'LineStyle', 'none','DisplayName', "Median Age", 'LineWidth', 1)
plot(median_depths2, median_ages2NaN, '-r', 'Marker', '.')
xlabel("Depth (cm)")
ylabel("Age (kyr)")
xlim([xlims])
ylim([ylims])

%Create new vector that combines SR estimate instead of removing
% filtOutInd3 = [filtOutInd2, filtOutInd2+1];
% nLogi = length(median_ages2);
% filtOutLogi3 = false(1,nLogi); filtOutLogi3(filtOutInd3) = true;
% median_depths3 = median_depths2(~filtOutLogi3); 
% median_ages3 = median_ages2(~filtOutLogi3); 
median_depths3 = median_depths2(~(median_depths2 == 252.5 | median_depths2 == 260.5 | median_depths2 == 334.5));
median_ages3 = median_ages2(~(median_depths2 == 252.5 | median_depths2 == 260.5 | median_depths2 == 334.5)); 

subplot(3,1,3)
hold on
plot(median_depths, median_ages, 'ko', 'LineStyle', 'none','DisplayName', "Median Age", 'LineWidth', 1)
plot(median_depths3, median_ages3, '-r', 'Marker', '.')
xlabel("Depth (cm)")
ylabel("Age (kyr)")
xlim([xlims])
ylim([ylims])

%Can I do this with a step plot showing the actual SRs?
SRs = median_nSRs{31}(2,2:end)./(median_nSRs{31}(3,2:end)./1000);
ages = cumsum(median_nSRs{31}(3,:))/1000;
aveSR = (median_depths(end)-median_depths(1))/(((median_ages(end)-median_ages(1))));
nSRs = SRs./aveSR;

%Get removed SRs vector
SRsNaN = SRs; SRsNaN(filtOutLogi) = NaN;
nSRsNaN = SRsNaN./aveSR;

%Get combined SRs vector
depth_and_age_gaps = median_nSRs{31}(2:3,2:end);
if sum(depth_and_age_gaps(2,:) < 500) > 0
    depth_and_age_gaps2 = NaN(size(depth_and_age_gaps));
    j = 1;
    skip_i = 0;
    for i = 1:size(depth_and_age_gaps, 2)
        if skip_i == 0;
            if depth_and_age_gaps(2,i) > 500
                depth_and_age_gaps2(:,j) = depth_and_age_gaps(:,i);
                j = j+1;
            else
                depth_and_age_gaps2(:,j) = sum(depth_and_age_gaps(:,[i i+1]), 2);
                
                skip_i = skip_i+1;
                while depth_and_age_gaps2(2,j) < 500
                    depth_and_age_gaps2(:,j) = sum(depth_and_age_gaps(:,[i:i+skip_i]), 2);
                    skip_i = skip_i+1;
                end
                j = j+1;
            end
        else
            skip_i = skip_i-1;
        end
        if i == 10
            a = 1;
        end
    end
end

depth_and_age_gaps = depth_and_age_gaps2;

combSRs = depth_and_age_gaps(1,:)./(depth_and_age_gaps(2,:)./1000);
combnSRs = combSRs./aveSR;
combages = cumsum([median_nSRs{31}(3,1), depth_and_age_gaps(2,:)])./1000;


figure
subplot(3,1,1)
stairs(ages, [nSRs nSRs(end)], '-k', 'LineWidth', 1)
ylabel('nSR')
xlabel('Age (kyr)')
xlim([0 20])

subplot(3,1,2)
stairs(ages, [nSRsNaN nSRsNaN(end)], '-k', 'LineWidth', 1)
ylabel('nSR')
xlabel('Age (kyr)')
xlim([0 20])
subplot(3,1,3)
stairs(combages, [combnSRs combnSRs(end)], '-k', 'LineWidth', 1)
ylabel('nSR')
xlabel('Age (kyr)')
xlim([0 20])
