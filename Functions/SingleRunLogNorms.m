function[SR_MixLogNorm1Run, c95up, c95down, mus,Sigmas, outS] = SingleRunLogNorms(nSRcounts, coreLog, numruns, x, numComponents,  weightRepDP, weightRepInflator, regularizationVal, fitS)

%Initialise vector
SR_MixLogNorm1Run = NaN(length(x), numruns);
weightedC = cell(1,numruns);
numCpairs = NaN(1,numruns);
chi2stats = NaN(1,numruns);
skippedIterations = 0;

%Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove    = cellfun(@isempty,nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~cores2remove);
runN            = NaN(numruns*2, length(nSRcountsClean));

% Take 1 run's worth of data from each core and calculate mix log normals.
for i = 1:numruns
    OneRunData = [];
    % For each core, go through and get the ith run's data
    skipIteration = 1;
    while skipIteration == 1
        for j = 1:length(nSRcountsClean)
            %Open the core's cell
            nSRcount_opencell = nSRcountsClean{j};

            %Find the NaN's that split runs
            split_index       = find(isnan(nSRcount_opencell(1,:)));
            split_indexStart  = split_index + 1;
            split_indexEnd    = split_index - 1;
            split_indexEnd    = [split_indexEnd(2:end), size(nSRcount_opencell, 2)];

            if i ==1
                runN(:,j) = randperm(length(split_index), numruns*2);  %%%This allows me to make the scenarios chosen random. The length(split_index) finds out how many times each core was sampled (when sampled previously, the sampling randomly chose scenarios). This chooses random samples (and therefore random scenarios).
            end

            % use the NaN indeces to choose the data of that run for each core
            % and concatenate it into OneRunData
            try
            nSRcount_1core1run = nSRcount_opencell(1:4,split_indexStart(runN(i,j)):split_indexEnd(runN(i, j)));
            catch
                disp("b")
            end
            OneRunData         = cat(2, OneRunData, nSRcount_1core1run);
        end

        if fitS.Lin2014AgeFiltering
            if sum(OneRunData<0) ~= 0
                warning("There are negative sed rates!")
            end
            L2014Log = OneRunData(4,:) < 4500 & OneRunData(4,:) > 500;
            OneRunData = OneRunData(:,L2014Log);
        end


        %% Weighting Data
        if fitS.weighting == "depth" %This indicates whether to weight by depth or not
            inputData       = OneRunData(1,:);
            weights         = OneRunData(3,:); %Weights by depth, not the mixed weighting of depth and scenarios
            inputDataClean  = inputData(~isnan(inputData));
            weightsClean    = weights(~isnan(weights));
            data            = makeWeightedReplicates(inputDataClean, weightsClean, weightRepDP, weightRepInflator); %Weight by replicating data according to weighting
        elseif fitS.weighting == "age"
            inputData       = OneRunData(1,:);
            weights         = OneRunData(4,:); %Weights by depth, not the mixed weighting of depth and scenarios
            inputDataClean  = inputData(~isnan(inputData));
            weightsClean    = weights(~isnan(weights));
            data            = makeWeightedReplicates(inputDataClean, weightsClean, weightRepDP, weightRepInflator); %Weight by replicating data according to weighting
        elseif fitS.weighting == "none"
            %Set up data as unweighted input data
            data = OneRunData(1,:);
        end

        % % Plot histogram and mix log norm for each run
        % if i < 10
        % figure;
        % hold on
        % yyaxis("left")
        % histogram(data, "BinEdges",[0:0.1:100, 1000])
        % xlim([0 10])
        % xlabel("nSR")
        % end

        %% Fitting mix log normal
        %Check that the data is a row vector
        [a,b] = size(data);
        if a>b
            data = data';
        end

        %Fit a mix log normal to the data. If there is an ill conditioned
        %covariance and an error is called, try again (initial conditions are
        %random so new initial conditions may allow the code to avoid an ill
        %conditioned covariance).
        allow = 0;
        errorNum = 0;
        while allow == 0 && skipIteration == 1
            try
                [SR_MixLogNorm1RunHOLDER, ~, gmfit] = fitMixLogNorm(data, x, numComponents, fitS.mln1RunReps); %If this is taking too long, try reducing the weightRepInflator value
                allow = 1;
                skipIteration = 0;
            catch %If there is an error in the fitMixLogNorm, which commonly happens because of an ill conditioned covariance, which may not occur if different initial conditions are used. If it continues to happen, skip this scenario and choose another scenario.
                errorNum = errorNum+1;
                disp("number of failed fittings = " + num2str(errorNum))
                allow = 0;
                if errorNum >10
                    skipIteration = 2;
                end
            end
            if mod(i, 10) == 0
            disp("Fit " + num2str(i) + " samplings")
            end
        end
    end

    if skipIteration == 2
        skippedIterations = skippedIterations + 1;
        %runN(i,:) = runN(numruns+skippedIterations, :);
    else
    %Get mus and sigmas from
    mus     = zeros(numruns,1);
    Sigmas  = zeros(numruns,1);
    if numComponents == 1
        mus(i)      = gmfit.mu;
        Sigmas(i)   = gmfit.Sigma;
    end

    %Run chi2gof on the fit
    if fitS.weighting ~= "none"
        numSRcalcs = numel(inputDataClean);
    else
        numSRcalcs = numel(data);
    end
    [~, ~, chiStat] = chi2gof_vsMLN(gmfit, log(data), numSRcalcs, fitS);

    weightedC{i} = data;
    numCpairs(i) = numSRcalcs;
    chi2stats(i) = chiStat.chi2stat;

    % % % Plot histogram and mix log norm for each run
    % if i < 10
    %     yyaxis("right")
    %     plot(x, SR_MixLogNorm1RunHOLDER(:,2), '-k', 'LineWidth', 1)
    %     ylabel("PDF")
    % end

    %Save mix log norm pdf values
    SR_MixLogNorm1Run(:, i) = SR_MixLogNorm1RunHOLDER(:,2);
    end
end

%Find 95% confidence intervals
sortedCols  = sort(SR_MixLogNorm1Run'); %#ok<TRSRT>
c95up       = sortedCols(ceil(numruns.*0.975), :);
c95down     = sortedCols(floor(numruns.*0.025), :);

%Put some useful info into an output structure
outS.weightedC = weightedC;
outS.numCpairs = numCpairs;
outS.chi2stats = chi2stats;
end