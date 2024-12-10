function[transnums, numCSE2x, coreTM, TMweighted] = TMcalculation(nSRcountsCell, coreSubsetLogical, S)
%This function takes in nSR counts and calculates a transition matrix from
%it.
%% Combine all counts into one array
nSRcountsArray = countsCell2Array(nSRcountsCell, coreSubsetLogical);
nSRcounts = nSRcountsArray(1,:);

%% --- Calculate TM without weights
%Categories are Steady, Expansion, Contraction (S, E, C)
char_categ = '';
for i = 1:(length(nSRcounts))
    if nSRcounts(i) >= 0.922 && nSRcounts(i) <1.085
        char_categ(i) = 'S';

    elseif nSRcounts(i)>=1.085 && nSRcounts(i) < inf
        char_categ(i) = 'E';

    elseif nSRcounts(i) >= 0 && nSRcounts(i) < 0.922
        char_categ(i) = 'C';

    else %If the normSR is not one of those values it'll be a NaN, which signifies the end of a run of ages.
        char_categ(i) = 'b';

    end
end

%----- Convert this to a transition matrix
%Find the number of times each transition occurs
transnums      = NaN(3,3);
transnums(1,1) = length(strfind(char_categ, "CC")); %Probability of Contraction to Contraction
transnums(1,2) = length(strfind(char_categ, "CS")); %Probability of Contraction to Steady
transnums(1,3) = length(strfind(char_categ, "CE")); %Probability of Contraction to Expansion
transnums(2,1) = length(strfind(char_categ, "SC")); %Probability of Steady to Contraction
transnums(2,2) = length(strfind(char_categ, "SS")); %Probability of Steady to Steady
transnums(2,3) = length(strfind(char_categ, "SE")); %Probability of Steady to Expansion
transnums(3,1) = length(strfind(char_categ, "EC")); %Probability of Expansion to Contraction
transnums(3,2) = length(strfind(char_categ, "ES")); %Probability of Expansion to Steady
transnums(3,3) = length(strfind(char_categ, "EE")); %Probability of Expansion to Expansion


%Find the number of times transitions occur from each state to any
%other state
numCSE2x    = zeros(3,1);
numCSE2x(1) = length(strfind(char_categ(1:end-1), "C")) - length(strfind(char_categ, "Cb"));
numCSE2x(2) = length(strfind(char_categ(1:end-1), "S"))- length(strfind(char_categ, "Sb"));
numCSE2x(3) = length(strfind(char_categ(1:end-1), "E"))- length(strfind(char_categ, "Eb"));

%calculate the transition matrix for this core
coreTM      = NaN(4,3);
coreTM(1,:) = numCSE2x./sum(numCSE2x);
coreTM(2,:) = transnums(1,:)./numCSE2x(1);
coreTM(3,:) = transnums(2,:)./numCSE2x(2);
coreTM(4,:) = transnums(3,:)./numCSE2x(3);

%% Calculate TM with weights

%Decide how to weight
if S.weighting == "depth"
    weights   = nSRcountsArray(2, :);
elseif S.weighting == "age"
    weights = nSRcountsArray(4,:)./1000; % dividing by 1000 to convert units from yr to kyr
else
    weights = ones(size(nSRcountsArray(2,:)));
end

%Set up empty character and weight vectors
charV = '';
weightV = ones(0,0);
%Decide which category the nSR count belongs to
for i = 1:length(nSRcounts)
    if nSRcounts(i) >= 0.922 && nSRcounts(i) <1.085
        charAdd = 'S';
    elseif nSRcounts(i)>=1.085 && nSRcounts(i) < inf
        charAdd = 'E';
    elseif nSRcounts(i) >= 0 && nSRcounts(i) < 0.922
        charAdd = 'C';
    else
        charAdd = 'b';
    end
    %If the nSR is in the same category as the previous one (i.e., the last
    %one added to the charV vector), just add the weights. If not, add a
    %new entry to the charV and weightV vectors.
    if i ~= 1 && charAdd == charV(end)
        weightV(end) = weightV(end) + weights(i);
    else
        charV = append(charV, charAdd);
        
        weightV = [weightV; weights(i)];
    end
end

% Find the number of times for each transition
transCounts = NaN(3,3);
transCounts(1,1) = sum(weightV(charV == 'C')) - sum(charV == 'C');
transCounts(1,2) = length(strfind(charV, "CS")); %Probability of Contraction to Steady
transCounts(1,3) = length(strfind(charV, "CE")); %Probability of Contraction to Expansion
transCounts(2,1) = length(strfind(charV, "SC")); %Probability of Steady to Contraction
transCounts(2,2) = sum(weightV(charV == 'S')) - sum(charV == 'S');
transCounts(2,3) = length(strfind(charV, "SE")); %Probability of Steady to Expansion
transCounts(3,1) = length(strfind(charV, "EC")); %Probability of Expansion to Contraction
transCounts(3,2) = length(strfind(charV, "ES")); %Probability of Expansion to Steady
transCounts(3,3) = sum(weightV(charV == 'E')) - sum(charV == 'E');

% Get the total number of transitions from any given state
numTrans = sum(transCounts, 2);

% Calculate the transition matrix
TM = NaN(4,3);
TM(1,:) = numTrans./sum(numTrans);
TM(2,:) = transCounts(1,:)./numTrans(1);
TM(3,:) = transCounts(2,:)./numTrans(2);
TM(4,:) = transCounts(3,:)./numTrans(3);

TMweighted = TM;

end