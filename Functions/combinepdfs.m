function[master_invSRvals, master_invSRprobs] = combinepdfs(cells_invSRvals, cells_invSRprobs, lengths)
%This function takes an array of cells, containing invSR and corresponding
%probs for a number of different cores/scenarios, and combines them into a
%single pdf.

%Get rid of scenarios/cores that didn't pass a filter
ind = ~cellfun('isempty',cells_invSRvals); %Create index of scenarios that passed filters
cells_invSRvalsNE= cells_invSRvals(ind); %Apply index to vals
cells_invSRprobs= cells_invSRprobs(ind); %Apply index to probs
numscenarios = length(cells_invSRvalsNE);%How many scenarios/cores passed the filters

if numscenarios == 0
    error("None of the cores used passed all the criteria, hence there are no pdfs of invSR to combine")
end

%Remove any nan values from lengths
lengths = lengths(~isnan(lengths));

%Create master vector of invSRvalues, which all scenarios will be placed
%onto
min_invSRval = min(cellfun(@min,cells_invSRvalsNE)); %find min value imported
max_invSRval = max(cellfun(@max,cells_invSRvalsNE)); %find max value imported
interval = cells_invSRvalsNE{1}(2)-cells_invSRvalsNE{1}(1); %define an interval
master_invSRvals = (min_invSRval:interval:(max_invSRval+interval/2))'; %create master vector

%Fit the probabilities of each individual scenario to the matching values on
%the master vector
m_invSRprobs = zeros(length(master_invSRvals),numscenarios); %initiate matrix to create master vector of probabilities
for i = 1:length(cells_invSRvalsNE)
    %idx1 = find(abs(master_invSRvals - cells_invSRvals{i}(1)) <=0.00001);
    idx = knnsearch(master_invSRvals, cells_invSRvalsNE{i}(1)); %Find where the invSRvalues of each scenario fit onto the master vector
    m_invSRprobs(idx:idx+length(cells_invSRvalsNE{i}(:))-1,i) = cells_invSRprobs{i}(:).*lengths(i); %Transfer probability values of each scenario into the matrix, such that their index corresponds to the index of their invSR value on the master vector of invSRvalues.
end

%Sum all invSRprobs up and divide by number of scenarios to get an average
%across all cores.
master_invSRprobs = sum(m_invSRprobs,2)./length(cells_invSRvalsNE);
end