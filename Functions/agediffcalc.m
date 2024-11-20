function[agediff_vals, agediff_probsums, IDpairs, agediffV] = agediffcalc(ageprob, ageVector, numpairs, labels, IDpairs, agediffV, S)
%This function calculates the age difference probability distribution
%functions for neighbouring pairs of radiocarbon dates. In order to reduce
%computational time, all calibrated radiocarbon pdfs have values reduced to
%0 (by entering NaN) where they are below some threshold value S.pdfMinVal,
%basically assuming those values are negligibly likely.

%Inputs
% ageVector - a 1D vector of length n, which holds the calendar year values
% that were interpolated onto by matcalq
% ageprob   - an n x m array where the mth column holds the calibration pdf of the
% mth radiocarbon date in this scenario
% S         - settings structure, holds value of pdfMinVal

%Do some cell initialisation to store each iterations values
agediff_vals = cell(1,numpairs);
agediff_probsums = cell(1,numpairs);

% input Nan where each pdf is less than S.pdfMinVal to reduce complexity
ind1 = ageprob(:,:)<=(S.pdfMinVal);
ageprob_Nans = ageprob;
ageprob_Nans(ind1) = NaN;

for m = 1:(numpairs)

    labelPair = labels(m) + "," + labels(m+1);
    pairCheckLog = ismember(IDpairs, labelPair);
    if sum(pairCheckLog) ~=0
        collector = agediffV{pairCheckLog};
        agediff_vals{m} = collector(:,1);
        agediff_probsums{m} = collector(:,2);
    else
        %Find vector of years for which the two radiocarbon dates are greater
    %than S.pdfMinVal (using ageprob_Nans)

    notNans1 = find(~isnan(ageprob_Nans(:,m))); %Find indices where there are no Nans in ageprob_Nans
    ages1 = ageVector(notNans1);
    probs1 = ageprob(notNans1, m);
    notNans2 = find(~isnan(ageprob_Nans(:,m+1)));
    ages2 = ageVector(notNans2);
    probs2 = ageprob(notNans2, m+1);

    %Calculate an agediff for each piecewise pairing of ages and the probability of that pairing (product
    %of probs of each age)
    agediffs = ages2'-ages1;
    pair_probs = probs1.*probs2';

    %Calculate the probability of each individual possibility of invSR_vals by
    %summing the probabilities in the pair_probs matrix with indices that
    %have that SR_val within the invSR matrix
    [agediff_vals{m},~,k2] = unique(agediffs);
    agediff_probsums{m} = accumarray(k2,pair_probs(:));

    IDpairs = [IDpairs; labelPair]; %#ok<AGROW>
    agediffV = [agediffV; {[agediff_vals{m}, agediff_probsums{m}]}]; %#ok<AGROW>

    % %Save the trouble of repeating construction of notNans1 next iteration
    % ages1 = ages2;
    % probs1 = probs2;
    end
end
end