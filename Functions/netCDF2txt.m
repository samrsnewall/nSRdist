function[] = netCDF2txt(corename, LabIDs, incDepths, excLabIDs, excDepths, dataLoc, separateBySR, S)
%% Read in Radiocarbon Data
%%Read in radiocarbon data from the core
if dataLoc == "WA"
    [age, depth_cm, error, label] = getDataWA(corename, S);
elseif dataLoc == "Lin2014"
    [age, depth_cm, error, label] = getDatatxt(corename, S);
end

%% Filtering
%Filter for MSPF dates, remove manually determined outliers, only keep
%dates between 1 and 42 14C ky BP, and only keep cores with 4 or more
%accepted dates
[age, depth_cm, error, ~, emptybreak1, emptybreak2] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename,  S);

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

outputFilename = fullfile(S.folderName, corename + ".txt");
writetable(outputTable, outputFilename, "Delimiter", '\t')

end