function[ageprob, calAge] = multiMatcal(age, error, date_is, S)
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

%Set up age vector
% calAge = 0:55000;
calAge = (0+1:3:55000-1)'; %Create vector of every 3 years, to reduce computational time

%Initialise ageprob array
numAges = length(date_is);
ageprob = zeros(length(calAge), numAges);

%Go through loop to calibrate each age, choosing how many ages you
for i = 1:numAges
    [~,~,holder,~] = matcal(age(date_is(i))*1000, error(date_is(i))*1000,  'Marine20', 'CalBP', 'resage', resAge, 'reserr', resError, 'plot', 0);
    matCalProbs = zeros(length(calAge), 1);
    for j = 1:55001/3
        matCalProbs(j) = sum(holder(j*3-2 : j*3, 2));
    end
    ageprob(:,i) = matCalProbs;
end


end