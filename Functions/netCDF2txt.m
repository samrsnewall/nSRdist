function[] = netCDF2txt(corename, LabIDs, incDepths, excLabIDs, excDepths, folderName, separateBySR)
%% Read in Radiocarbon Data
%Set up path to Mulitza et al 2022 World Atlas Dataset
WA_path = "/Applications/PaleoDataView/WA_Foraminiferal_Isotopes_2022";
fnm = fullfile(WA_path, "Age/", corename + ".age");
%Read in radiocarbon data from the core
depth_m = ncread(fnm, "Depth"); %(meters)
depth_cm = depth_m.*100; %convert to cm
age = ncread(fnm, "Age dated"); %(14C kyrs BP)
error = ncread (fnm, "Age +Error"); %(14C kyrs BP)
label = ncread(fnm, "Label"); %(Lab ID)
label = string(label);

%% Filtering
%Filter for MSPF dates, remove manually determined outliers, only keep
%dates between 1 and 42 14C ky BP, and only keep cores with 4 or more
%accepted dates
[age, depth_cm, error, ~, emptybreak1, emptybreak2] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths);

if emptybreak1 ==1 || emptybreak2 == 1
    return
end

%% Filter for SR>8 cores
if separateBySR ==1
    meanSR = (depth_cm(end) - depth_cm(1))/(age(end) - age(1));
    if meanSR < 8
        return
    end
end

%% Output in csv file in format Taehee wants
%convert to correct units (depth in m, age in ky)
depth = depth_cm/100;
%set up dR, dSTD and cc columns
dR = zeros(length(age), 1);
dSTD = ones(length(age), 1)*0.2;
cc = ones(length(age), 1)*2;
%Set up table
outputTable = table(depth, age, error, dR, dSTD, cc);
%Output to txt file with tab delimiters

outputFilename = fullfile(folderName, corename + ".txt");
writetable(outputTable, outputFilename, "Delimiter", '\t')

end