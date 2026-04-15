function[age, depth_cm, error, label] = getDatatxt(corename, S)
% getDatatxt  Read radiocarbon age data from a Bchron-formatted .txt file.
%
% Locates and reads the radiocarbon input file for a named sediment core
% from the Lin2014 dataset. The file is expected to be a tab-delimited
% table with columns: id (string), Age (yr), Error (yr), Depth (cm). Ages
% and errors are converted from years to kyr on output. Rows are sorted by
% ascending depth, and depths are rounded to 2 decimal places to avoid
% floating-point inconsistencies.
%
% INPUTS
%   corename  - (string) Sediment core identifier used to locate the file.
%               The function searches for any file in the BchronInput folder
%               whose name contains corename.
%   S         - (struct) Settings struct. Field used:
%                 .sandboxPath  Root path of the repository; the data file
%                               is expected at:
%                               <sandboxPath>/Lin2014Cores/BchronInput/
%
% OUTPUTS
%   age       - (numeric vector) Conventional radiocarbon ages (kyr BP),
%               sorted by ascending depth
%   depth_cm  - (numeric vector) Sample depths (cm), sorted ascending,
%               rounded to 2 decimal places
%   error     - (numeric vector) 1-sigma age uncertainties (kyr)
%   label     - (string vector) Laboratory IDs for each date
%
% See also: getDataWA, filtering, oneCoreScenarios, nSRBchron

%Set path to folder with all txt files
txt_path = fullfile(S.sandboxPath, "Lin2014Cores", "BchronInput");

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
