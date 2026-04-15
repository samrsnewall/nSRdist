function[transnums, numCSE2x, TM, TMw] = TMcalculationAR(nSRcountsAR,  fitS)
% TMcalculationAR  Compute a 3-state sedimentation transition matrix from
%                  an nSR counts array.
%
% Classifies each nSR value into one of three states — Contraction (C),
% Steady (S), or Expansion (E) — using the thresholds from Lin et al.
% (2014): C = [0, 0.922), S = [0.922, 1.085), E = [1.085, inf). NaN
% values (run-break markers) are assigned the state 'b' and are
% excluded from transition counts. The function then counts all consecutive-
% pair transitions between non-break states to build the unweighted 4×3
% transition matrix TM.
%
% Row layout of TM:
%   Row 1: state frequency vector [fC, fS, fE]   (fraction of transitions
%          originating from each state)
%   Row 2: transition probabilities out of C      [P(C→C), P(C→S), P(C→E)]
%   Row 3: transition probabilities out of S      [P(S→C), P(S→S), P(S→E)]
%   Row 4: transition probabilities out of E      [P(E→C), P(E→S), P(E→E)]
%
% INPUTS
%   nSRcountsAR - (3 × N numeric matrix) nSR matrix in the standard format
%                 (see README "Internal Data Formats"). Only row 1 (nSR
%                 values) is used for the unweighted TM; rows 2–3
%                 (depth/age differences) are available for the weighted
%                 calculation below.
%   fitS        - (struct) Fitting settings struct. Reserved for the
%                 weighted transition matrix calculation (see commented-out
%                 code below); not used by the active code path. Relevant
%                 field when the weighted section is re-enabled:
%                   .weighting  "depth" | "age" | "none"
%
% OUTPUTS
%   transnums  - (3 × 3 numeric matrix) Raw counts of each pairwise
%                transition (rows = from-state, cols = to-state; order: C, S, E)
%   numCSE2x   - (3 × 1 numeric vector) Total number of transitions
%                originating from each state [nC; nS; nE]
%   TM         - (4 × 3 numeric matrix) Unweighted transition matrix; see
%                row layout above
%   TMw        - (4 × 3 or empty) Weighted transition matrix. Currently
%                returned as [] because it is unclear from Lin et al. (2014)
%                whether weighting was applied before or after TM
%                construction. The weighted implementation is retained below
%                in commented form for re-verification if needed.
%
% See also: ARfitdists, fitData
%% Get NSR values
nSRcounts = nSRcountsAR(1,:);

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
TM      = NaN(4,3);
TM(1,:) = numCSE2x./sum(numCSE2x);
TM(2,:) = transnums(1,:)./numCSE2x(1);
TM(3,:) = transnums(2,:)./numCSE2x(2);
TM(4,:) = transnums(3,:)./numCSE2x(3);

%% Calculate TM with weights
TMw = [];
% 
% %Set up weights choice
% if fitS.weighting == "depth"
%     weights   = nSRcountsAR(2, :);
% elseif fitS.weighting == "age"
%     weights = nSRcountsAR(3,:)./1000; % dividing by 1000 to convert units from yr to kyr
% else
%     TMw = TM;
% end
% 
% if fitS.weighting ~= "none"
%     %Set up empty character and weight vectors
%     charV = '';
%     weightV = ones(0,0);
%     %Decide which category the nSR count belongs to
%     for i = 1:length(nSRcounts)
%         if nSRcounts(i) >= 0.922 && nSRcounts(i) <1.085
%             charAdd = 'S';
%         elseif nSRcounts(i)>=1.085 && nSRcounts(i) < inf
%             charAdd = 'E';
%         elseif nSRcounts(i) >= 0 && nSRcounts(i) < 0.922
%             charAdd = 'C';
%         else
%             charAdd = 'b';
%         end
%         %If the nSR is in the same category as the previous one (i.e., the last
%         %one added to the charV vector), just add the weights. If not, add a
%         %new entry to the charV and weightV vectors.
%         if i ~= 1 && charAdd == charV(end)
%             weightV(end) = weightV(end) + weights(i);
%         else
%             charV = append(charV, charAdd);
%             weightV = [weightV; weights(i)]; %#ok<AGROW>
%         end
%     end
% 
%     % Find the number of times for each transition
%     transCounts = NaN(3,3);
%     transCounts(1,1) = sum(weightV(charV == 'C')) - sum(charV == 'C');
%     transCounts(1,2) = length(strfind(charV, "CS")); %Probability of Contraction to Steady
%     transCounts(1,3) = length(strfind(charV, "CE")); %Probability of Contraction to Expansion
%     transCounts(2,1) = length(strfind(charV, "SC")); %Probability of Steady to Contraction
%     transCounts(2,2) = sum(weightV(charV == 'S')) - sum(charV == 'S');
%     transCounts(2,3) = length(strfind(charV, "SE")); %Probability of Steady to Expansion
%     transCounts(3,1) = length(strfind(charV, "EC")); %Probability of Expansion to Contraction
%     transCounts(3,2) = length(strfind(charV, "ES")); %Probability of Expansion to Steady
%     transCounts(3,3) = sum(weightV(charV == 'E')) - sum(charV == 'E');
% 
%     % Get the total number of transitions from any given state
%     numTrans = sum(transCounts, 2);
% 
%     % Calculate the transition matrix (weighted)
%     TMw = NaN(4,3);
%     TMw(1,:) = numTrans./sum(numTrans);
%     TMw(2,:) = transCounts(1,:)./numTrans(1);
%     TMw(3,:) = transCounts(2,:)./numTrans(2);
%     TMw(4,:) = transCounts(3,:)./numTrans(3);
% end

end