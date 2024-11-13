function[agediff_vals, agediff_probsums, deldep] = agediffcalc(ageprob, ageVector, dep_is)
%Calculate age difference and probability for each pairing of ages - PDF METHOD
%Do some cell initialisation to store each iterations values
numpairs = length(dep_is)-1;

agediff_vals = cell(1,numpairs);
agediff_probsums = cell(1,numpairs);
deldep = zeros(numpairs,1);

% input Nan where each pdf is lesser than 1e-7
ind1 = ageprob(:,:)<=(1e-6);
ageprob_Nans = ageprob;
ageprob_Nans(ind1) = NaN;

for m = 1:(numpairs)

    %Find vector of years for which the two radiocarbon dates are greater
    %than 1e-10 (using ageprob_Nans)
    notNans1 = find(~isnan(ageprob_Nans(:,m))); %Find indices where there are no Nans in ageprob_Nans
    notNans2 = find(~isnan(ageprob_Nans(:,m+1)));
    ages1 = ageVector(notNans1); %Utilise fact that the indices from the MatCal output are the ages +1 (ages start at 0, increase by 1 year, indices start at 1).
    ages2 = ageVector(notNans2);
    probs1 = ageprob(notNans1, m);
    probs2 = ageprob(notNans2, m+1);

    %find the depth difference between the two ages in cm
    deldep(m) = (dep_is(m+1)-dep_is(m));

    %Calculate an agediff for each piecewise pairing of ages and the probability of that pairing (product
    %of probs of each age)
    agediffs = ages2'-ages1;
    pair_probs = probs1.*probs2';

    %Calculate the probability of each individual possibility of invSR_vals by
    %summing the probabilities in the pair_probs matrix with indices that
    %have that SR_val within the invSR matrix
    [agediff_vals{m},~,k2] = unique(agediffs);
    agediff_probsums{m} = accumarray(k2,pair_probs(:));
end
end