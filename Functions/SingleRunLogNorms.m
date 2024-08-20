function[SR_MixLogNorm1Run] = SingleRunLogNorms(nSRcounts, coreLog, numruns, x)

%Initialise vector
SR_MixLogNorm1Run = NaN(length(x), numruns);

%Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove = cellfun(@isempty,nSRcountsChosen);
nSRcountsClean = nSRcountsChosen(~cores2remove);

% Take 1 run's worth of data from each core and calculate mix log normals.
for i = 1:numruns
OneRunData = [];
    % For each core, go through and get the ith run's data
    for j = 1:length(nSRcountsClean)

      %Open the core's cell
      nSRcount_openfile = nSRcountsClean{j};

      %Find the NaN's that split runs
      split_index = find(isnan(nSRcount_openfile(1,:)));
      split_indexStart = split_index + 1;
      split_indexEnd = split_index - 1;

      % use the NaN indeces to choose the data of that run for each core
      % and concatenate it into OneRunData
      nSRcount_1core1run = nSRcount_openfile(1:2,split_indexStart(i):split_indexEnd(i+1));
      OneRunData = cat(2, OneRunData, nSRcount_1core1run);
    end

    %% MUST ADD THIS WEIGHTING COMPONENT!
    X = OneRunData(1,:)';
    Y = OneRunData(2,:);
    Xclean = X(~isnan(X)); %Remove NaNs that separate cores and runs
    Yclean = Y(~isnan(X)); %Remove NaNs that separate cores and runs
    [X_u, ~, IC] = unique(Xclean);
    Y_u = accumarray(IC,Yclean);
    Y_uR = round(Y_u.*100); %%%%%%%%%% NEED TO DO A SENSITIVITY TEST TO THIS! IF I DON'T MULTIPY Y_u BY 100 THEN MUCH DATA HAS ITS WEIGHTING ROUNDED DOWN TO 0 WEIGHT... TRY
    data = repelem(X_u,Y_uR);

    % %Fit a mixture Log Normal to the data (not weighted)
    % [SR_MixLogNorm1RunHOLDER, ~, ~] = fitMixLogNorm(OneRunData(1,:), 2);

    %Fit a mixture Log Normal to the data (weighted)
    [SR_MixLogNorm1RunHOLDERweighted, ~, ~] = fitMixLogNorm(data, 2);

    %Plot histogram and mix log norm for each run
    if i < 6
        figure(i);
        hold on
        histogram(data, "Normalization", "pdf", "BinEdges",[0:0.1:100, 1000])
        plot(SR_MixLogNorm1RunHOLDERweighted(:,1), SR_MixLogNorm1RunHOLDERweighted(:,2), '-k', 'LineWidth', 1)
        xlim([0 10])
    end

%Save mix log norm pdf values 
  SR_MixLogNorm1Run(:, i) = SR_MixLogNorm1RunHOLDERweighted(:,2);
end