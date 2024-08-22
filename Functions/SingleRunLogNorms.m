function[SR_MixLogNorm1Run] = SingleRunLogNorms(nSRcounts, coreLog, numruns, x, weightQ, weightRepDP, weightRepInflator, regularizationVal)

%Initialise vector
SR_MixLogNorm1Run = NaN(length(x), numruns);

%Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove    = cellfun(@isempty,nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~cores2remove);

% Take 1 run's worth of data from each core and calculate mix log normals.
for i = 1:numruns
OneRunData = [];
    % For each core, go through and get the ith run's data
    for j = 1:length(nSRcountsClean)

      %Open the core's cell
      nSRcount_opencell = nSRcountsClean{j};

      %Find the NaN's that split runs
      split_index       = find(isnan(nSRcount_opencell(1,:)));
      split_indexStart  = split_index + 1;
      split_indexEnd    = split_index - 1;

      % use the NaN indeces to choose the data of that run for each core
      % and concatenate it into OneRunData
      nSRcount_1core1run = nSRcount_opencell(1:2,split_indexStart(i):split_indexEnd(i+1));
      OneRunData         = cat(2, OneRunData, nSRcount_1core1run);
    end

    %% MUST ADD THIS WEIGHTING COMPONENT!
    if weightQ == 1 %This indicates whether to weight by depth or not
        inputData       = OneRunData(1,:);
        weights         = OneRunData(2,:);
        inputDataClean  = inputData(~isnan(inputData));
        weightsClean    = weights(~isnan(weights));
        data            = makeWeightedReplicates(inputDataClean, weightsClean, weightRepDP, weightRepInflator); %Weight by replicating data according to weighting
    else
        %Set up data as unweighted input data
        data = OneRunData(1,:);
    end
    %Fit a mix log normal to the data
    [SR_MixLogNorm1RunHOLDER, ~, ~] = fitMixLogNorm(data, x, 2, regularizationVal);

    % %Plot histogram and mix log norm for each run
    % if i < 6
    %     figure;
    %     hold on
    %     histogram(data, "Normalization", "pdf", "BinEdges",[0:0.1:100, 1000])
    %     plot(x, SR_MixLogNorm1RunHOLDER(:,2), '-k', 'LineWidth', 1)
    %     xlim([0 10])
    %     ylim([0 1.5])
    %     xlabel("nSR")
    % end

%Save mix log norm pdf values 
  SR_MixLogNorm1Run(:, i) = SR_MixLogNorm1RunHOLDER(:,2);
end