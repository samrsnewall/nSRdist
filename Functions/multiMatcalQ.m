function [ageprob, calAge] = multiMatcalQ(age, error, date_is, S)
% multiMatcal  Calibrate multiple radiocarbon ages using the Marine20
%               calibration curve (via MatCal).
%
% Calls matcal once per date and assembles the resulting probability
% density functions into a single matrix. A reservoir age correction of
% zero is applied (i.e. the Marine20 curve's built-in reservoir correction
% is used as-is), and the user-specified DeltaR uncertainty is propagated
% through calibration.
%
% INPUTS
%   age      - (numeric vector, n×1) Conventional radiocarbon ages (kyr BP)
%   error    - (numeric vector, n×1) 1-sigma uncertainties on each age
%              (kyr)
%   date_is  - (numeric vector) Indices into age/error selecting which
%              dates to calibrate. Allows a subset to be calibrated without
%              re-passing the full arrays.
%   S        - (struct) Settings struct. Field used:
%                .DeltaRError  (numeric) 1-sigma uncertainty on the
%                              reservoir age correction (years). A
%                              DeltaR of 0 is assumed; only the error
%                              is propagated.
%                .matcalFast   (boolean) whether to use the editted matcal
%                version for speed
%
% OUTPUTS
%   ageprob - (matrix, 55001 × numel(date_is)) Each column holds the
%             calibrated probability density function (PDF) for one
%             date, on the shared calAge axis. Row i corresponds to
%             calAge(i) calendar years BP.
%   calAge  - (numeric vector, 55001×1) Calendar year axis (years BP)
%             corresponding to the rows of ageprob, as returned by
%             matcal for the last calibrated date.
%
% See also: agediffcalc, scenariopdfNorm, oneCoreScenarios

%% Set up reservoir age correction
resAge   = 0;                  % DeltaR is assumed zero; only the error is used
resError = S.DeltaRError;      % 1-sigma uncertainty on DeltaR (years)

%% Calibrate first age to get length of cal age vector
% Input ages must be converted from kyr to years before passing to matcal.
if ~S.matcalFast
    [~, ~, holder, ~] = matcal(age(date_is(1)) * 1000, error(date_is(1)) * 1000, ...
        'Marine20', 'CalBP', 'resage', resAge, 'reserr', resError, 'plot', 0);
else
    [~, ~, holder, ~] = matcalFast(age(date_is(1)) * 1000, error(date_is(1)) * 1000, ...
        'Marine20', 'CalBP', 'resage', resAge, 'reserr', resError, 'plot', 0);
end

% calAge is identical for every matcal call (same calibration curve and
% axis), so it is read from the first call only.
calAge = holder(:,1);

% Initialize vector to hold ages from all other dates
lengthAgeV = length(calAge);
numAges     = length(date_is);
ageprob     = zeros(lengthAgeV, numAges);
ageprob(:,1)= holder(:,2);

%Calibrate the rest of the ages
for i = 2:numAges
    if ~S.matcalFast
        [~, ~, holder, ~] = matcal(age(date_is(i)) * 1000, error(date_is(i)) * 1000, ...
            'Marine20', 'CalBP', 'resage', resAge, 'reserr', resError, 'plot', 0);
    else
        [~, ~, holder, ~] = matcalFast(age(date_is(1)) * 1000, error(date_is(1)) * 1000, ...
            'Marine20', 'CalBP', 'resage', resAge, 'reserr', resError, 'plot', 0);
    end
    ageprob(:,i) = holder(:,2);
end
end
