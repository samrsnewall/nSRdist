function[age, depth_cm, error, label] = getDataWA(corename, S)
% getDataWA  Read radiocarbon age data from a World Atlas .age NetCDF file.
%
% Reads the radiocarbon dates for a named sediment core from the World
% Atlas (WA) dataset. Depths are converted from metres to centimetres and
% rounded to 2 decimal places to avoid floating-point inconsistencies.
% Ages and errors are already in units of 14C kyr BP in the source files
% and are returned as-is.
%
% INPUTS
%   corename  - (string) Sediment core identifier (e.g. "RC13-228"). The
%               function constructs the filename as:
%               <S.WApath>/Age/<corename>.age
%   S         - (struct) Settings struct. Field used:
%                 .WApath  Path to the World Atlas data directory
%
% OUTPUTS
%   age       - (numeric vector) Conventional radiocarbon ages (kyr BP)
%   depth_cm  - (numeric vector) Sample depths (cm), rounded to 2 decimal
%               places
%   error     - (numeric vector) 1-sigma age uncertainties (kyr BP)
%   label     - (string vector) Laboratory IDs for each date
%
% See also: getDatatxt, filtering, oneCoreScenarios, nSRBchron

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