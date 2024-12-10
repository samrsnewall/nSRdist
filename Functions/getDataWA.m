function[age, depth_cm, error, label] = getDataWA(corename, S)
fnm = fullfile(S.WApath, "Age/", corename + ".age");
%Read in radiocarbon data from the core
depth_m     = ncread(fnm, "Depth"); %(meters)
depth_cm    = depth_m.*100; %convert to cm
depth_cm    = round(depth_cm, 2); %Round to the nearest 100th of a cm to avoid floating point number problems

age         = ncread(fnm, "Age dated"); %(14C kyrs BP)

error       = ncread (fnm, "Age +Error"); %(14C kyrs BP)

label       = ncread(fnm, "Label"); %(Lab ID)
label       = string(label);
end