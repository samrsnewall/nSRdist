function[SR_MixLogNorm1Run, c95up, c95down, mus,Sigmas] = SingleRunLogNorms(nSRcounts, coreLog, numruns, x, numComponents, weightQ, weightRepDP, weightRepInflator, regularizationVal, fitS)

%Initialise vector
SR_MixLogNorm1Run = NaN(length(x), numruns);

%Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove    = cellfun(@isempty,nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~cores2remove);
runN            = NaN(numruns, length(nSRcountsClean));

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

        if i ==1
            runN(:,j) = randperm(length(split_index)-1, numruns);  %%%This allows me to make the scenarios chosen random

        end

        % use the NaN indeces to choose the data of that run for each core
        % and concatenate it into OneRunData
        nSRcount_1core1run = nSRcount_opencell(1:4,split_indexStart(runN(i,j)):split_indexEnd(runN(i, j)+1));
        OneRunData         = cat(2, OneRunData, nSRcount_1core1run);
    end

    if fitS.Lin2014AgeFiltering
        L2014Log = OneRunData(4,:) < 4500 & OneRunData(4,:) > 500;
        OneRunData = OneRunData(:,L2014Log);
    end


    %% Weighting Data
    if weightQ == 1 %This indicates whether to weight by depth or not
        inputData       = OneRunData(1,:);
        weights         = OneRunData(3,:); %Weights by depth, not the mixed weighting of depth and scenarios
        inputDataClean  = inputData(~isnan(inputData));
        weightsClean    = weights(~isnan(weights));
        data            = makeWeightedReplicates(inputDataClean, weightsClean, weightRepDP, weightRepInflator); %Weight by replicating data according to weighting
    else
        %Set up data as unweighted input data
        data = OneRunData(1,:);
    end

    %     % Plot histogram and mix log norm for each run
    % if i < 10
    %     figure;
    %     hold on
    %     histogram(data, "Normalization", "pdf", "BinEdges",[0:0.1:100, 1000])
    %     xlim([0 10])
    %     ylim([0 1.5])
    %     xlabel("nSR")
    % end

    %% Fitting mix log normal
    %Check that the data is a row vector
    [a,b] = size(data);
    if a>b
        data = data';
    end
    

    %Fit a mix log normal to the data
    [SR_MixLogNorm1RunHOLDER, ~, gmfit] = fitMixLogNorm(data, x, numComponents, regularizationVal, 1); %If this is taking too long, try reducing the weightRepInflator value

    %Get mus and sigmas from 
    mus     = zeros(numruns,1);
    Sigmas  = zeros(numruns,1);
    if numComponents == 1
        mus(i)      = gmfit.mu;
        Sigmas(i)   = gmfit.Sigma;
    end

    % % Plot histogram and mix log norm for each run
    % if i < 10
    %     plot(x, SR_MixLogNorm1RunHOLDER(:,2), '-k', 'LineWidth', 1)
    % end

    %Save mix log norm pdf values
    SR_MixLogNorm1Run(:, i) = SR_MixLogNorm1RunHOLDER(:,2);
end

%Find 95% confidence intervals
sortedCols  = sort(SR_MixLogNorm1Run');
c95up       = sortedCols(numruns.*0.975, :);
c95down     = sortedCols(numruns.*0.025, :);

end