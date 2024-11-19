function[ageprob, calAge] = multiMatcalQ(age, error, date_is, S)
%This function calibrates multiple radiocarbon ages and returns an array
%ageprob that holds the probability densities for each ith age in its ith
%column. The array has 55001 rows, each row corresponds to the year row-1 
%so the 13th row is the probability density of that age at calendar year 12.

%Note, the input 'date_is' contains which dates, from all dates contained in
%the input 'age' to be calibrated and returned. This allows for different
%scenarios to be considered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up reservoir age information
resAge = 0;
resError = S.DeltaRError;

%Set up length of age vector (known from before)
lengthAgeVQ = 18334;

%Initialise ageprob array
numAges = length(date_is);
ageprob = zeros(lengthAgeVQ, numAges);

%Go through loop to calibrate each age, choosing how many ages you
for i = 1:numAges
    [~,~,holder,~] = matcalq(age(date_is(i))*1000, error(date_is(i))*1000,  'Marine20', 'CalBP', 'resage', resAge, 'reserr', resError, 'plot', 0);
    ageprob(:,i) = holder(:,2);
end
calAge = holder(:,1);
end