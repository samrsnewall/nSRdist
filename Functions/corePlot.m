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
[ChosenMSPF, EM, ManE, NotCCR, ~,~] = filteringForPlotting(age, depth, error, label, LabIDs, incDepths, excLabIDs, excDepths);

%% Plot all radiocarbon data
figure
errorbar(ChosenMSPF.depth, ChosenMSPF.age, ChosenMSPF.error, "vertical", 'o', "color", 'k', 'DisplayName', "Accepted")
errorbar(EM.depth, EM.age, EM.error, "vertical", '.', 'color', 'r', 'DisplayName', 'Not Chosen MSPF')
errorbar(ManE.depth, ManE.age, ManE.error, "vertical", 'o', 'color', 'r', 'DisplayName', "Manually Excluded")
errorbar(NotCCR.depth, NotCCR.age, NotCCR.error, "vertical", 'x', 'color', 'r', 'DisplayName', "Outside CC Range")
set(gca, 'YTickLabel',get(gca,'YTick'))
xlabel("Depth (cm)")
ylabel(["Radiocarbon Age","(14C kyr BP)"])
legend()
title(corename)

end