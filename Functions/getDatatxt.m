function[age, depth_cm, error, label] = getDatatxt(corename)

%Set path to folder with all txt files
txt_path = "/Users/samnewall/Documents/MATLAB/nSRdist_code/Lin2014Cores";

%Get names of all files
openfile = dir(txt_path);
filenames = strings(length(openfile),1);
for i = 1:length(openfile)
    filenames(i) = string(openfile(i).name);
end

%Find file that has core name
fileIndex   = contains(filenames, corename);
filename    = filenames(fileIndex);

%Open file and read in data, ensuring the id column is read in as a string
opts = detectImportOptions(fullfile(txt_path, filename));
opts = setvartype(opts, "id", "string");
dataT = readtable(fullfile(txt_path, filename), opts);

%Ensure data is ordered by depth, shallowest first
dataT = sortrows(dataT, "Depth");

%Pick out data of interest
age_yrs     = dataT.Age;
depth_cm    = dataT.Depth;
error_yrs   = dataT.Error;
label       = string(dataT.id);

%Perform any necessary unit conversion
age         = age_yrs./1000;
error       = error_yrs./1000;

%Fix for floating point errors
depth_cm = round(depth_cm, 2);
end
