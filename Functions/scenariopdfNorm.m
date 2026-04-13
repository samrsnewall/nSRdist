function [interp_invSR, iProbs_wdN, aveSR, reversalpairs, numpairs, ageMode, depdiff, agediff, MSI_byage, MSI_bydepth, IDpairs, agediffV] = scenariopdfNorm(depth, date_is, label, ageprob, calAge, IDpairs, agediffV, S, plotfigs)
% scenariopdfNorm  Compute the normalised inverse-SR PDF for a single
%                  age-depth scenario, and detect age reversals.
%
% This is the core PDF-method computation. For a given set of dates (one
% age-depth scenario), it:
%   1. Identifies which consecutive date pairs show probable age reversals.
%   2. If any reversals are found, returns early (outputs zeroed) so the
%      caller (scenariosDealWithReversals) can split the scenario.
%   3. Otherwise, calls agediffcalc to compute the age-difference PDF for
%      each pair, converts to inverse SR, normalises by the average SR of the
%      scenario, interpolates onto a common grid, then depth-weights and
%      sums the per-pair PDFs into a single scenario PDF.
%
% NORMALISATION
% Inverse SR values are normalised by the average SR (computed as total
% sediment length divided by total time span, in cm/kyr) so that the
% resulting PDF is in dimensionless units of normalised inverse SR (nSR = 1
% at the average SR). This allows PDFs from different cores and time periods
% to be combined.
%
% INPUTS
%   depth    - (numeric vector) Depth of every date in the core (cm).
%              Indexed by date_is to get depths for this scenario.
%   date_is  - (numeric vector) Indices of the dates in this scenario,
%              pointing into depth, label, and ageprob.
%   label    - (string vector) Lab IDs for all dates in the core
%   ageprob  - (matrix, n × numel(date_is)) Calibrated age PDFs for the
%              dates in this scenario (rows = calendar years, one column
%              per date). Produced by multiMatcalQ.
%   calAge   - (numeric vector, n×1) Calendar year axis (years BP) for the
%              rows of ageprob
%   IDpairs  - (string vector) Age-difference computation cache (see
%              agediffcalc for details)
%   agediffV - (cell array) Age-difference result cache
%   S        - (struct) Settings struct. Fields used:
%                .pdfMinVal        Probability threshold for trimming PDFs
%                .reversalCriteria Fraction of probability mass below 0 that
%                                  classifies a pair as a reversal (default 0.75)
%                .pdfMethod        (logical) Whether to build the full
%                                  interpolated PDF or just compute summary
%                                  statistics
%                .weighting        "depth" | "age" | "none" — how to weight
%                                  per-pair PDFs before summing
%   plotfigs - (logical) Produce diagnostic plots if true
%
% OUTPUTS
%   interp_invSR  - (numeric vector) Common nSR grid onto which all per-pair
%                   PDFs were interpolated. Empty if S.pdfMethod is false or
%                   reversals were found.
%   iProbs_wdN    - (numeric vector) Normalised, depth-weighted sum of all
%                   per-pair PDFs on the interp_invSR grid. Empty if
%                   S.pdfMethod is false or reversals were found.
%   aveSR        - (numeric scalar) Average SR for this scenario (cm/kyr),
%                   computed as total depth span / total time span using
%                   mode ages.
%   reversalpairs - (logical vector, 1×numpairs) Entry m is 1 if the pair
%                   (m, m+1) is classified as a reversal.
%   numpairs      - (integer) Number of consecutive date pairs
%   ageMode       - (numeric vector) Mode calendar age for each date (yr BP)
%   depdiff       - (numeric scalar) Total depth span of the scenario (cm)
%   agediff       - (numeric scalar) Total age span using mode ages (yr)
%   MSI_byage     - (numeric scalar) Mean sampling interval by age (yr/date)
%   MSI_bydepth   - (numeric scalar) Mean sampling interval by depth (cm/date)
%   IDpairs       - Updated age-difference cache
%   agediffV      - Updated age-difference cache
%
% See also: agediffcalc, multiMatcalQ, combinepdfs, scenariosDealWithReversals

%% Convert calendar age axis to kyr
calAgekyrs = calAge ./ 1000;

%% Mode age of each date and summary statistics
% Find the most-probable calendar age for each date.
ageMode = zeros(length(date_is), 1);
for i = 1:length(date_is)
    [~, modeInd] = max(ageprob(:,i));
    ageMode(i)   = calAge(modeInd);
end

dep_is  = depth(date_is);    % Depths of the dates in this scenario (cm)
agediff = ageMode(end) - ageMode(1);    % Total time span (yr)
depdiff = dep_is(end) - dep_is(1);      % Total depth span (cm)

% Mean sampling intervals
MSI_byage   = agediff  / length(date_is);   % yr/date
MSI_bydepth = depdiff  / length(date_is);   % cm/date

%% Compute age-difference PDFs for each consecutive date pair
numpairs   = length(date_is) - 1;
label_is   = label(date_is);
[agediff_vals, agediff_probsums, IDpairs, agediffV] = agediffcalc(ageprob, calAge, numpairs, label_is, IDpairs, agediffV, S);

deldep = diff(dep_is);    % Depth difference for each pair (cm)

%% Detect age reversals
% A pair is classified as a reversal if more than S.reversalCriteria of its
% age-difference probability mass falls at or below zero.
reversalpairs = zeros(1, numpairs);
for m = 1:numpairs
    if sum(agediff_probsums{m}(agediff_vals{m} <= 0)) >= S.reversalCriteria
        reversalpairs(m) = 1;
    end
end

% Return early if any reversals exist: the caller will split this scenario.
if sum(reversalpairs) > 0
    interp_invSR = NaN;
    iProbs_wdN   = NaN;
    aveSR = NaN;
    return
end

%% Restrict age differences to positive values only
% Any residual probability at non-positive age differences is discarded and
% the remaining distribution is renormalised.
agediff_probsumsPos = cell(1, numpairs);
agediff_valsPos     = cell(1, numpairs);
for m = 1:numpairs
    gt0log                  = agediff_vals{m} > 0;
    agediff_valsPos{m}      = agediff_vals{m}(gt0log);
    agediff_probsumsPos{m}  = agediff_probsums{m}(gt0log);
    agediff_probsumsPos{m}  = agediff_probsumsPos{m} ./ sum(agediff_probsumsPos{m}, 'all');
end

%% Convert age differences to inverse SR (yr/cm) and re-normalise PDFs
invSR_vals     = cell(1, numpairs);
invSR_probsums = cell(1, numpairs);
invSRAUC       = NaN(1, numpairs);
for m = 1:numpairs
    invSR_vals{m}     = agediff_valsPos{m} ./ deldep(m);   % yr/cm
    invSRAUC(m)       = trapz(invSR_vals{m}, agediff_probsumsPos{m});
    invSR_probsums{m} = agediff_probsumsPos{m} ./ invSRAUC(m);
end

%% Compute average SR for this scenario (cm/kyr)
% Uses the mode ages of the outermost dates and the total depth span.
[~, topage_ind]    = max(ageprob(:,1));
[~, bottomage_ind] = max(ageprob(:,end));
topage_mode     = calAgekyrs(topage_ind);
bottomage_mode  = calAgekyrs(bottomage_ind);
total_agediff   = bottomage_mode - topage_mode;   % kyr
total_depthdiff = dep_is(end) - dep_is(1);        % cm
aveSR          = total_depthdiff ./ total_agediff;   % cm/kyr

%% Normalise inverse SR by average SR
% Multiplying inverse SR (yr/cm) by average SR (cm/kyr) and a unit-conversion
% factor (1/1000 yr/kyr) gives dimensionless nSR = invSR / (1/aveSR).
invSR_normvals = cellfun(@(x) x * (aveSR .* (1/1000)), invSR_vals, 'un', 0);

%% Diagnostic plots (individual pairwise PDFs)
if plotfigs == 1
    pdfCreation = figure("Name", "pdfCreation");

    subplot(4,1,1)
    for n = 1:numpairs
        plot(invSR_vals{1,n} ./ 1000, invSR_probsums{1,n}, 'LineWidth', 1)
        hold on
    end
    ylabel("Probability")
    xlabel("Inverse of Sed Rate ky/cm")

    subplot(4,1,2)
    hold on
    for n = 1:numpairs
        plot(invSR_normvals{1,n}, invSR_probsums{1,n})
    end
    ylabel("Probability")
    xlabel("Normalised Inverse Sed Rate")
end

%% Interpolate, align, weight and sum per-pair PDFs into a scenario PDF
if S.pdfMethod

    % Interpolate each pair's PDF onto a common grid with fixed spacing.
    spacing = 0.00005;
    vector_invSRvals = vertcat(invSR_normvals{:});
    lowerbound   = floor(min(vector_invSRvals * (1/spacing))) * spacing;
    upperbound   = max(vector_invSRvals);
    interp_invSR = lowerbound:spacing:(upperbound + spacing * 0.9);

    interp_invSR_indy  = cell(1, numpairs);
    interp_invSRp_indy = cell(1, numpairs);
    for nd = 1:numpairs
        if ~isempty(invSR_normvals{nd})
            min_val = min(invSR_normvals{nd});
            max_val = max(invSR_normvals{nd});
            ind1    = find(interp_invSR > min_val, 1, 'first');
            ind2    = find(interp_invSR < max_val, 1, 'last');
            interp_invSR_indy{nd}      = interp_invSR(ind1:ind2);
            interp_invSRp_indy_hold    = interp1(invSR_normvals{nd}, invSR_probsums{nd}, interp_invSR(ind1:ind2));
            interp_invSRp_indy{nd}     = interp_invSRp_indy_hold ./ sum(interp_invSRp_indy_hold);
        end
    end

    % Place each pair's interpolated PDF into a shared matrix aligned to the
    % master grid, then apply depth weighting.
    mat_iProbs = zeros(length(interp_invSR), numpairs);
    for i = 1:numpairs
        if ~isempty(interp_invSR_indy{i})
            idx = find(abs(interp_invSR - interp_invSR_indy{i}(1)) <= 0.00001);
            mat_iProbs(idx:idx + length(interp_invSR_indy{i}(:)) - 1, i) = interp_invSRp_indy{i}(:);
        end
    end

    if S.weighting == "depth"
        mat_iProbs_wd = mat_iProbs .* deldep';
    else
        % "age" weighting is not yet implemented; "none" uses equal weights.
        mat_iProbs_wd = mat_iProbs;
    end

    if plotfigs == 1
        figure(pdfCreation)
        subplot(4,1,3)
        hold on
        plot(interp_invSR, mat_iProbs_wd)
        xlabel("Individual invSR Ratios weighted")
        ylabel("Probability")
    end

    % Sum across pairs and normalise.
    iProbs_wd  = sum(mat_iProbs_wd, 2);
    AUC        = trapz(interp_invSR, iProbs_wd);
    iProbs_wdN = iProbs_wd ./ AUC;

    if plotfigs == 1
        figure(pdfCreation)
        subplot(4,1,4)
        hold on
        plot(interp_invSR, iProbs_wdN)
        ylabel("Probability")
        xlabel("Inverse of Sed Rate Ratio (Weighted by depth)")
    end

else
    interp_invSR = [];
    iProbs_wdN   = [];
end

end
