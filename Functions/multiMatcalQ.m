function [ageprob, calAge] = multiMatcalQ(age, error, date_is, S)
% multiMatcalQ  Calibrate multiple radiocarbon ages using the Marine20
%               calibration curve (via MatCal).
%
% Calls matcalq once per date and assembles the resulting probability
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
%
% OUTPUTS
%   ageprob - (matrix, 18334 × numel(date_is)) Each column holds the
%             calibrated probability density function (PDF) for one
%             date, on the shared calAge axis. Row i corresponds to
%             calAge(i) calendar years BP.
%             Note: 18334 is the fixed length of the Marine20 output
%             vector returned by matcalq and should not be changed
%             manually.
%   calAge  - (numeric vector, 18334×1) Calendar year axis (years BP)
%             corresponding to the rows of ageprob, as returned by
%             matcalq for the last calibrated date.
%
% See also: agediffcalc, scenariopdfNorm, oneCoreScenarios

%% Set up reservoir age correction
resAge   = 0;                  % DeltaR is assumed zero; only the error is used
resError = S.DeltaRError;      % 1-sigma uncertainty on DeltaR (years)

%% Calibrate each selected date
% matcalq returns a fixed-length vector (18334 rows) for the Marine20 curve.
% Input ages must be converted from kyr to years before passing to matcalq.
lengthAgeVQ = 18334;
numAges     = length(date_is);
ageprob     = zeros(lengthAgeVQ, numAges);

for i = 1:numAges
    [~, ~, holder, ~] = matcalq(age(date_is(i)) * 1000, error(date_is(i)) * 1000, ...
        'Marine20', 'CalBP', 'resage', resAge, 'reserr', resError, 'plot', 0);
    ageprob(:,i) = holder(:,2);
end

% calAge is identical for every matcalq call (same calibration curve and
% axis), so it is read from the last call only.
calAge = holder(:,1);

end
