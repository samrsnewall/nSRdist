%% A7 Explanatory Figures
%This script provides some example sets of figures about the sampling
%procedures, with core A7 and SHAK06-5K as examples. 

%% Plot A7 age depth model probabilistic
dA7 = load("../Results/dataT_RLGtrue_R200M20_Apr23_A7_fitApr23_depthweight.mat");
A7Bchron.input = readtable("../BchronFolders/Bchron_RLGtrue_R200M20_May14_DS0p05/Outputs/A7/inputData.txt");
A7Bchron.theta = readmatrix("../BchronFolders/Bchron_RLGtrue_R200M20_May14_DS0p05/Outputs/A7/theta.csv");

%% Plot calibrated radiocarbon ages
figure("Position", [100 100 400 300]);
hold on
plot(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none','DisplayName', "Mode Age", 'LineWidth', 1)
errorbar(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, A7Bchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none', 'DisplayName', "95% Uncertainty", 'LineWidth', 1)
xlabel("Depth (cm)")
ylabel("Age (cal kyr)")
xlim([0 500])
%plot([A7Bchron.input.Depth(1), A7Bchron.input.Depth(end)], [d1C.d.dataT.ageModes{1}{1}(1)./1000, d1C.d.dataT.ageModes{1}{1}(end)./1000], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
%plot(A7Bchron.input.Depth, cumsum(d1C.d.dataT.bchronMode{1}(3,:))./1000, 'LineStyle','none', 'LineWidth', 1, "DisplayName", "BchronMode", "Marker","square", "Color",'k', "MarkerFaceColor","r")
legend("Location", "SouthEast")
% plot(d1C.d.dataT)

%% Plot BMedian
%Get nSR vector
median_nSRs = dA7.d.dataT.bchronMedian;

%Get median ages and depths
median_ages = cumsum(median_nSRs{1}(3,:));
median_depths = cumsum(median_nSRs{1}(2,:));

%Construct repeated vectors for plotting with NaNs
median_depths2 = repelem(median_depths, 2);
median_ages2 = repelem(median_ages, 2);

figure;
subplot(3,1,1)
hold on
plot(median_depths, median_ages, 'ko', 'LineStyle', 'none','DisplayName', "Median Age", 'LineWidth', 1)
plot(median_depths2, median_ages2, 'r', 'Marker', 'none')
xlabel("Depth (cm)")
ylabel("Age (kyr)")
xlim([175 300])
ylim([11000 14500])

%Find indices of nSRs that have dt < 500;
filtOutInd = find(median_nSRs{1}(3,:) < 500);
filtOutInd2 = filtOutInd*2 - 1;

%Put NaN in repeated vectors so they nSRs of dt <500 do not plot
median_ages2NaN = median_ages2;
median_ages2NaN(filtOutInd2) = NaN;

subplot(3,1,2)
hold on
plot(median_depths, median_ages, 'ko', 'LineStyle', 'none','DisplayName', "Median Age", 'LineWidth', 1)
plot(median_depths2, median_ages2NaN, '-r', 'Marker', 'none')
xlabel("Depth (cm)")
ylabel("Age (kyr)")
xlim([175 300])
ylim([11000 14500])

%Create new vector that combines SR estimate instead of removing
filtOutInd3 = [filtOutInd2, filtOutInd2+1];
nLogi = length(median_ages2);
filtOutLogi3 = false(1,nLogi); filtOutLogi3(filtOutInd3) = true;
median_depths3 = median_depths2(~filtOutLogi3); 
median_ages3 = median_ages2(~filtOutLogi3); 

subplot(3,1,3)
hold on
plot(median_depths, median_ages, 'ko', 'LineStyle', 'none','DisplayName', "Median Age", 'LineWidth', 1)
plot(median_depths3, median_ages3, '-r', 'Marker', 'none')
xlabel("Depth (cm)")
ylabel("Age (kyr)")
xlim([175 300])
ylim([11000 14500])

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



%% Plot 4 potential age histories with BchronProb (subplots)
figure;
for i = 1:4
subplot(2,2,i)
hold on
plot(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none','DisplayName', "Mode Age", 'LineWidth', 1)
errorbar(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, A7Bchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none', 'DisplayName', "95% Uncertainty", 'LineWidth', 1)
xlabel("Depth (cm)")
ylabel("Age (cal kyr)")
MCMCagehistory = cumsum(dA7.d.dataT.bchronProb{1}(4,i.*length(A7Bchron.input.Depth)+1:(i+1).*length(A7Bchron.input.Depth)))./1000;
plot([A7Bchron.input.Depth(1), A7Bchron.input.Depth(end)], [MCMCagehistory(1), MCMCagehistory(end)], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
plot(A7Bchron.input.Depth, MCMCagehistory, 'LineStyle','-', 'LineWidth', 1, "DisplayName", "BchronMCMC" + num2str(i), "Marker","square", "Color",'b', "MarkerFaceColor","k")
legend("Location", "SouthEast")
end
% plot(d1C.d.dataT)

%% Plot 20 individual BchronProb age histories (individual plots)

for i = 1:20
figure("Position", [100 100 400 300]);
hold on
plot(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none','DisplayName', "Mode Age", 'LineWidth', 1)
errorbar(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, A7Bchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none', 'DisplayName', "95% Uncertainty", 'LineWidth', 1)
xlabel("Depth (cm)")
ylabel("Age (cal kyr)")
MCMCagehistory = cumsum(dA7.d.dataT.bchronProb{1}(4,i.*length(A7Bchron.input.Depth)+1:(i+1).*length(A7Bchron.input.Depth)))./1000;
%plot([A7Bchron.input.Depth(1), A7Bchron.input.Depth(end)], [MCMCagehistory(1), MCMCagehistory(end)], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
plot(A7Bchron.input.Depth, MCMCagehistory, 'LineStyle','none', 'LineWidth', 1, "DisplayName", "BchronMCMC" + num2str(i), "Marker","square", "Color",'k', "MarkerFaceColor","r")
legend("Location", "SouthEast")
imageName = "/Users/samnewall/Desktop/A7BSamp/fig" + num2str(i) + ".jpg";
%saveas(gcf, imageName)
end
% plot(d1C.d.dataT)

%% Plot 20 individual RSR0 age histories (individual plots)
for i = 1:20
figure("Position", [100 100 400 300]);
hold on
plot(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none','DisplayName', "Mode Age", 'LineWidth', 1)
errorbar(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, A7Bchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none', 'DisplayName', "95% Uncertainty", 'LineWidth', 1)
xlabel("Depth (cm)")
ylabel("Age (cal kyr)")
MCMCagehistory = cumsum(dA7.d.dataT.nSRcounts{1}(4,i.*length(A7Bchron.input.Depth)+1:(i+1).*length(A7Bchron.input.Depth)))./1000;
%plot([A7Bchron.input.Depth(1), A7Bchron.input.Depth(end)], [MCMCagehistory(1), MCMCagehistory(end)], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
plot(A7Bchron.input.Depth, MCMCagehistory, 'LineStyle','none', 'LineWidth', 1, "DisplayName", "RSR0: " + num2str(i), "Marker","square", "Color",'k', "MarkerFaceColor","r")
legend("Location", "SouthEast")
imageName = "/Users/samnewall/Desktop/A7RSR0/fig" + num2str(i) + ".jpg";
%saveas(gcf, imageName)
end
% plot(d1C.d.dataT)

%% Plot 20 individual RSR1000 age histories (individual plots)
for i = 1:20
figure("Position", [100 100 400 300]);
hold on
plot(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none','DisplayName', "Mode Age", 'LineWidth', 1)
errorbar(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, A7Bchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none', 'DisplayName', "95% Uncertainty", 'LineWidth', 1)
xlabel("Depth (cm)")
ylabel("Age (cal kyr)")
MCMChistorybreaks = find(isnan(dA7.d.dataT.nSRcounts1000{1}(1,:)));
MCMChistoryindexes = [MCMChistorybreaks(i), MCMChistorybreaks(i+1)-1];
MCMCdepthhistory = cumsum(dA7.d.dataT.nSRcounts1000{1}(3,MCMChistoryindexes(1):MCMChistoryindexes(2)));
MCMCagehistory = cumsum(dA7.d.dataT.nSRcounts1000{1}(4,MCMChistoryindexes(1):MCMChistoryindexes(2)))./1000;
%plot([A7Bchron.input.Depth(1), A7Bchron.input.Depth(end)], [MCMCagehistory(1), MCMCagehistory(end)], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
plot(MCMCdepthhistory, MCMCagehistory, 'LineStyle','none', 'LineWidth', 1, "DisplayName", "RSR1000: " + num2str(i), "Marker","square", "Color",'k', "MarkerFaceColor","r")
legend("Location", "SouthEast")
imageName = "/Users/samnewall/Desktop/A7RSR1000/fig" + num2str(i) + ".jpg";
%saveas(gcf, imageName)
end
% plot(d1C.d.dataT)

%% Same for Shak06-5K
dSHAK = load("../Results/dataT_All1_RLGtrue_DS0p05_SHAK065K_Jul16_fitJul16_depthweight.mat");
dSHAKmindt0 = load("../Results/dataT_All1_RLGtrue_DS0p05_SHAK065K_Jul16_fitJul17_depthweight_mindt0.mat");
SHAKBchron.input = readtable("../BchronFolders/All1_RLGtrue_DS0p05_Jun2/Outputs/SHAK06-5K/inputData.txt");

%% Plot calibrated radiocarbon ages
figure("Position", [100 100 400 300]);
hold on
plot(SHAKBchron.input.Depth, dSHAK.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none','DisplayName', "Mode Age", 'LineWidth', 1)
errorbar(SHAKBchron.input.Depth, dSHAK.d.dataT.ageModes{1}{1}./1000, SHAKBchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none', 'DisplayName', "95% Uncertainty", 'LineWidth', 1)
xlabel("Depth (cm)")
ylabel("Age (cal kyr)")
xlim([0 350])
%plot([SHAKBchron.input.Depth(1), SHAKBchron.input.Depth(end)], [dSHAK.d.dataT.ageModes{1}{1}(1)./1000, dSHAK.d.dataT.ageModes{1}{1}(end)./1000], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
plot(SHAKBchron.input.Depth, cumsum(dSHAK.d.dataT.bchronMode{1}(3,:))./1000, 'LineStyle','none', 'LineWidth', 1, "DisplayName", "BchronMode", "Marker","square", "Color",'k', "MarkerFaceColor","r")
legend("Location", "SouthEast")

%% Plot A7 and SHAK06-5K bchronMode age history
% plot(d1C.d.dataT)
figure("Position", [100 100 400 300]);
hold on
plot(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none', 'LineWidth', 1, 'HandleVisibility','off')
errorbar(A7Bchron.input.Depth, dA7.d.dataT.ageModes{1}{1}./1000, A7Bchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none', 'LineWidth', 1, 'HandleVisibility','off')
xlabel("Depth (cm)")
ylabel("Age (cal kyr)")
xlim([0 500])
%plot([A7Bchron.input.Depth(1), A7Bchron.input.Depth(end)], [d1C.d.dataT.ageModes{1}{1}(1)./1000, d1C.d.dataT.ageModes{1}{1}(end)./1000], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
plot(A7Bchron.input.Depth, cumsum(dA7.d.dataT.bchronMode{1}(3,:))./1000, 'LineStyle','none', 'LineWidth', 1, "DisplayName", "A7", "Marker","square", "Color",'k', "MarkerFaceColor","r")
plot(SHAKBchron.input.Depth, dSHAK.d.dataT.ageModes{1}{1}./1000, 'ko', 'LineStyle', 'none', 'LineWidth', 1, 'HandleVisibility','off')
errorbar(SHAKBchron.input.Depth, dSHAK.d.dataT.ageModes{1}{1}./1000, SHAKBchron.input.Error./100, 'vertical', 'k', 'LineStyle' ,'none',  'LineWidth', 1, 'HandleVisibility','off')
%plot([SHAKBchron.input.Depth(1), SHAKBchron.input.Depth(end)], [dSHAK.d.dataT.ageModes{1}{1}(1)./1000, dSHAK.d.dataT.ageModes{1}{1}(end)./1000], 'r--', 'LineWidth', 1, 'DisplayName',  "MeanSR")
plot(SHAKBchron.input.Depth, cumsum(dSHAK.d.dataT.bchronMode{1}(3,:))./1000, 'LineStyle','none', 'LineWidth', 1, "DisplayName", "SHAK06-5K", "Marker","square", "Color",'k', "MarkerFaceColor","b")
legend("Location", "SouthEast")

%% Plot histogram of nSR counts from both cores BchronMode runs
SHAK0andA7 = [dA7.d.S1.BMode.weightedC, dSHAK.d.S1.BMode.weightedC];
figure; 
hold on
histogram(SHAK0andA7, "BinEdges",binEdges, 'FaceColor', 'b')
histogram(dA7.d.S1.BMode.weightedC, "BinEdges",binEdges, 'FaceColor', 'r', 'FaceAlpha', 1)
xlim([0 6])
xlabel("nSR")
ylabel("cm")
title("nSR weighted counts histogram")

%% Plot histograms of nSR from BMode run, showing impact of removing data where dt < 500 (as done in Lin2014)
figure; 
subplot(2,1,1)
histogram(dA7.d.S1.BMode.weightedC, "BinEdges",binEdges, 'FaceColor', 'r')
xlim([0 6])
ylim([0 100])
xlabel("nSR")
ylabel("cm")
subplot(2,1,2)
hold on
histogram(dSHAKmindt0.d.S1.BMode.weightedC, "BinEdges",binEdges, 'FaceColor', 'k', 'FaceAlpha', 0.1)
histogram(dSHAK.d.S1.BMode.weightedC, "BinEdges",binEdges, 'FaceColor', 'b')
xlim([0 6])
ylim([0 100])
xlabel("nSR")
ylabel("cm")

figure; 
subplot(2,1,1)
histogram(dA7.d.S1.BMode.weightedC, "BinEdges",binEdges, 'FaceColor', 'r')
xlim([0 6])
ylim([0 100])
xlabel("nSR")
ylabel("cm")
subplot(2,1,2)
histogram(dSHAK.d.S1.BMode.weightedC, "BinEdges",binEdges, 'FaceColor', 'b')
xlim([0 6])
ylim([0 100])
xlabel("nSR")
ylabel("cm")