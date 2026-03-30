function [agediff_vals, agediff_probsums, IDpairs, agediffV] = agediffcalc(ageprob, ageVector, numpairs, labels, IDpairs, agediffV, S)
% agediffcalc  Compute the age-difference PDF for each consecutive pair of
%              radiocarbon dates in a scenario.
%
% For each neighbouring date pair (m, m+1), this function convolves their
% individual calibrated PDFs to produce a PDF over possible age differences
% (age_{m+1} - age_m). The convolution is performed analytically: every
% pairwise product of individual ages is enumerated, age differences are
% computed, and probabilities with the same age difference are summed.
%
% To reduce computation time, values in each calibrated PDF below
% S.pdfMinVal are first set to NaN so that the convolution skips
% negligibly-probable ages.
%
% CACHING
% Computing the age-difference PDF for a given date pair is expensive. To
% avoid repeating this work when the same pair appears in multiple
% scenarios, previously computed results are cached in IDpairs and
% agediffV and re-used on subsequent calls. The cache is passed in and out
% as function arguments so it persists across calls within one core.
%
% INPUTS
%   ageprob    - (matrix, n × (numpairs+1)) Calibrated age PDFs for the
%                dates in this scenario. Column i is the PDF for date i.
%   ageVector  - (numeric vector, n×1) Calendar year axis shared by all
%                columns of ageprob (years BP)
%   numpairs   - (integer) Number of consecutive date pairs, i.e.
%                (number of dates in scenario) - 1
%   labels     - (string vector, (numpairs+1)×1) Lab IDs for the dates in
%                this scenario, used as cache keys
%   IDpairs    - (string vector) Cache: lab ID pair strings already
%                computed, e.g. "ODP1234,ODP1235"
%   agediffV   - (cell array) Cache: each entry is a two-column matrix
%                [agediff_vals, agediff_probsums] for the corresponding
%                entry in IDpairs
%   S          - (struct) Settings struct. Field used:
%                  .pdfMinVal  Values below this threshold in each
%                              calibrated PDF are treated as zero to
%                              speed up convolution.
%
% OUTPUTS
%   agediff_vals     - (cell array, 1×numpairs) Each cell contains a vector
%                      of unique age-difference values (years) for one pair
%   agediff_probsums - (cell array, 1×numpairs) Each cell contains the
%                      summed probability for each entry in agediff_vals
%   IDpairs          - Updated cache of computed pair IDs
%   agediffV         - Updated cache of computed age-difference distributions
%
% See also: scenariopdfNorm, multiMatcalQ

%% Initialise output cells
agediff_vals     = cell(1, numpairs);
agediff_probsums = cell(1, numpairs);

%% Threshold calibrated PDFs to skip negligibly-probable ages
ind1           = ageprob(:,:) <= S.pdfMinVal;
ageprob_Nans   = ageprob;
ageprob_Nans(ind1) = NaN;

%% Compute (or retrieve from cache) the age-difference PDF for each pair
for m = 1:numpairs

    labelPair    = labels(m) + "," + labels(m+1);
    pairCheckLog = ismember(IDpairs, labelPair);

    if sum(pairCheckLog) ~= 0
        %------------------------------------------------------------------
        % Cache hit: retrieve previously computed result
        %------------------------------------------------------------------
        collector            = agediffV{pairCheckLog};
        agediff_vals{m}      = collector(:,1);
        agediff_probsums{m}  = collector(:,2);

    else
        %------------------------------------------------------------------
        % Cache miss: compute age-difference PDF from scratch
        %------------------------------------------------------------------
        % Extract the non-negligible ages and probabilities for each date.
        notNans1 = find(~isnan(ageprob_Nans(:,m)));
        ages1    = ageVector(notNans1);
        probs1   = ageprob(notNans1, m);

        notNans2 = find(~isnan(ageprob_Nans(:,m+1)));
        ages2    = ageVector(notNans2);
        probs2   = ageprob(notNans2, m+1);

        % Enumerate all pairwise age differences and their joint probabilities.
        agediffs   = ages2' - ages1;          % outer difference matrix
        pair_probs = probs1 .* probs2';       % joint probability matrix

        % Collapse to a 1-D distribution over unique age-difference values.
        [agediff_vals{m}, ~, k2] = unique(agediffs);
        agediff_probsums{m}      = accumarray(k2, pair_probs(:));

        % Store result in cache for re-use across scenarios.
        IDpairs  = [IDpairs;  labelPair];                               %#ok<AGROW>
        agediffV = [agediffV; {[agediff_vals{m}, agediff_probsums{m}]}]; %#ok<AGROW>
    end
end

end
