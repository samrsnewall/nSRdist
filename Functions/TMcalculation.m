function[transnums, numCSE2x, coreTM] = TMcalculation(normSRs)
%This function takes in nSR counts and calculates a transition matrix from
%it

%---- Categorise normalised sed rates
%Categories are Steady, Expansion, Contraction (S, E, C)
char_categ = '';
for i = 1:(length(normSRs))
    if normSRs(i) >= 0.922 && normSRs(i) <1.085
        char_categ(i) = 'S';

    elseif normSRs(i)>=1.085 && normSRs(i) < inf
        char_categ(i) = 'E';

    elseif normSRs(i) >= 0 && normSRs(i) < 0.922
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
end