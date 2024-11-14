function[invSR_Gamma1Run, phat95up, phat95down, phat, nSR_Gamma1Run] = SingleRunGammas(nSRcounts, coreLog, numruns, x, weightQ, weightRepDP, weightRepInflator, regularizationVal)

%Initialise vector
invSR_Gamma1Run = NaN(length(x), numruns);
nSR_Gamma1Run   = NaN(length(x), numruns);
phat            = NaN(2,numruns);

%Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove    = cellfun(@isempty,nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~cores2remove);

%Set up inverse nSR grid
invx = sort(1./x);


% Take 1 run's worth of data from each core and calculate inverse gammas
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

    %% MUST ADD THIS WEIGHTING COMPONENT!
    if weightQ == 1 %This indicates whether to weight by depth or not
        inputData       = OneRunData(1,:);
        weights         = OneRunData(3,:);
        inputDataClean  = inputData(~isnan(inputData));
        weightsClean    = weights(~isnan(weights));
        data            = makeWeightedReplicates(inputDataClean, weightsClean, weightRepDP, weightRepInflator); %Weight by replicating data according to weighting
    else
        %Set up data as unweighted input data
        data = OneRunData(1,:);
    end



    %convert nSR to inverse nSR
    dataInv = 1./data;

    % Estimate the gamma fit parameters
    phat(:,i) = gamfit(dataInv);

    % Create gamma on invx values
    invxGamProb = gampdf(invx, phat(1,i), phat(2,i));

    % Plot inverse data
    binEdges = 0:0.1:15;

    % if i < 6
    %     figure;
    %     hold on
    %     histogram(dataInv, "BinEdges", binEdges, 'Normalization', "pdf");
    %     plot(invx, invxGamProb, 'LineWidth', 2)
    %     xlim([0 5])
    %     xlabel("Inverse normalised SR")
    %     title("Shape value = " + num2str(phat(1,i)))
    % end

    % Convert back to nSR
    [nSR_Gamma, nSR_GammaProb] = gammaAccRate2nSR(phat(1,i), phat(2,i), invx);
    % [nSRtest, nSRtestProb] = invSRtoSR(invx, invxGamProb); %Gives the same result as the line above
    
    %interpolate onto x
    nSR_GammaProb = interp1(nSR_Gamma, nSR_GammaProb, x);

    %Save gamma pdf values
    invSR_Gamma1Run(:, i)   = invxGamProb;
    nSR_Gamma1Run(:,i)      = nSR_GammaProb;
end

% %Find 95% of pdfs
% sortedCols = sort(invSR_Gamma1Run');
% c95up = sortedCols(numruns.*0.975, :);
% c95down = sortedCols(numruns.*0.025, :);

%Find 95% confidence values of phat
sortedPhat_alpha = sort(phat(1,:), 'ascend');
phat95up        = sortedPhat_alpha(numruns.*0.975);
phat95down      = sortedPhat_alpha(numruns.*0.025);

end