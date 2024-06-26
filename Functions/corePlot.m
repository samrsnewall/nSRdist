function[] = corePlot(corename, LabIDs, incDepths, excLabIDs, excDepths)
%% Read in Radiocarbon Data
%Read in some radiocarbon data from a net cdf file
WA_path = "/Applications/PaleoDataView/WA_Foraminiferal_Isotopes_2022";
fnm = fullfile(WA_path, "Age/", corename + ".age");

%Read in radiocarbon data from the core
depth_m = ncread(fnm, "Depth"); %(meters)
depth = depth_m.*100; %convert to cm
age = ncread(fnm, "Age dated"); %(14C kyrs BP)
error = ncread (fnm, "Age +Error"); %(14C kyrs BP)
label = ncread(fnm, "Label"); %(Lab ID)
label = string(label);

%% Filter for plotting
[age, depth, error, ~, ~,~, EM, ManE, NotCCR, ageGapTooHigh] = filtering(age, depth, error, label, LabIDs, incDepths, excLabIDs, excDepths);


%% Plot all radiocarbon data
figure
hold on
errorbar(depth, age, error, "vertical", 'o', "color", 'k', 'DisplayName', "Accepted")
errorbar(EM.depth, EM.age, EM.error, "vertical", '.', 'color', 'r', 'DisplayName', 'Material Not Chosen')
errorbar(ManE.depth, ManE.age, ManE.error, "vertical", 'o', 'color', 'r', 'DisplayName', "Manually Excluded")
errorbar(NotCCR.depth, NotCCR.age, NotCCR.error, "vertical", 'x', 'color', 'r', 'DisplayName', "Outside CC Range")
plot(ageGapTooHigh.depth, ageGapTooHigh.age, 'r-')
%ylim([0 age(end)*1.1])
set(gca, 'YTickLabel',get(gca,'YTick'))
xlabel("Depth (cm)")
ylabel(["Radiocarbon Age","(14C kyr BP)"])
legend("location", "southeast")
title(corename)

end