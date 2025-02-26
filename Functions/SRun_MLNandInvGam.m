function[outS] = SRun_MLNandInvGam(nSRcounts, coreLog, numruns, x, MLNcomponents, fitS)

%Initialise vectors
MLNfits = cell(length(x), numruns);
MLNs = NaN(length(x), numruns);
IGfits = cell(length(x), numruns);
invGams = NaN(length(x), numruns);
%phats  = NaN(2,numruns);
weightedC = cell(1,numruns);
numCpairs = NaN(1,numruns);
%chi2stats = NaN(1,numruns);
h1 = NaN(1,numruns);
p1 = NaN(1,numruns);
h2 = NaN(1,numruns);
p2 = NaN(1,numruns);
skippedIterations = 0;

%Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove    = cellfun(@isempty,nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~cores2remove);
runN            = NaN(numruns*2, length(nSRcountsClean));

%Initialise chiStat tables
MLNchiStatT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
InvGamChiStatT = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);

% Take 1 run's worth of data from each core and fits distributions
for i = 1:numruns
    OneRunData = [];
    % For each core, go through and get the ith run's data
    allow1 = 0;
    skipIteration = 0;
    while skipIteration == 0 & allow1 == 0
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

            nSRcount_1core1run = nSRcount_opencell(1:4,split_indexStart(runN(i,j)):split_indexEnd(runN(i, j)));
            
            if sum(nSRcount_1core1run(1,:) <=0) ~=0
                disp("check this core")
            end

            OneRunData         = cat(2, OneRunData, nSRcount_1core1run);
        end


        %Get rid of negative sedimentation rates (ideally, shouldn't be
        %there)
        if sum(OneRunData(1,:)<0) ~= 0
            warning("There are negative sed rates!... being removed")
            OneRunData = OneRunData(:, OneRunData(1,:) > 0);
        end
        
        %Get rid of SRs of 0 (ideally, shouldn't be there)
        if sum(OneRunData(1,:) == 0) ~= 0
            warning("There are SRs of 0!... being removed.")
            OneRunData = OneRunData(:,OneRunData(1,:) ~= 0);
        end

        %Apply filtering done by Lin2014 if wanted
        if fitS.Lin2014AgeFiltering
            L2014Log = OneRunData(4,:) < 4500 & OneRunData(4,:) > 500;
            OneRunData = OneRunData(:,L2014Log);
        end


        %% Weighting Data
        weightRepDP = fitS.OneRun.weightRepDP;
        weightRepInflator = fitS.OneRun.weightRepInflator;

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

        %Find number desired number of datapoints for chi2gof later
        if fitS.weighting ~= "none"
            numSRcalcs = numel(inputDataClean);
        else
            numSRcalcs = numel(data)./weightRepInflator;
        end

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
        allow2 = 0;
        errorNum = 0;
        while allow2 == 0 && skipIteration == 0
            try
                
                [SR_MixLogNorm1RunHOLDER, ~, gmfitH] = fitMixLogNorm(data, x, MLNcomponents, fitS.OneRun.MLNReps); %If this is taking too long, try reducing the weightRepInflator value
                [invGamPDFHOLDER, GamPDFHOLDER,IGfitH] = fitInvGamma(data, x);

                %Remove any NaNs in HOLDERS (think this is a bug?)
                SR_MixLogNorm1RunHOLDER(isnan(SR_MixLogNorm1RunHOLDER)) = 0;
                invGamPDFHOLDER(isnan(invGamPDFHOLDER)) = 0;
                GamPDFHOLDER(isnan(GamPDFHOLDER)) = 0;

                allow1 = 1;
                allow2 = 1;
                skipIteration = 0;
            catch %If there is an error in the fitMixLogNorm, which commonly happens because of an ill conditioned covariance, which may not occur if different initial conditions are used. If it continues to happen, skip this scenario and choose another scenario.
                errorNum = errorNum+1;
                disp("number of failed fittings = " + num2str(errorNum))
                allow2 = 0;
                if errorNum >10
                    skipIteration = 1;
                    gmfitH = NaN;
                    SR_MixLogNorm1RunHOLDER = NaN(length(x), 2);
                    IGfitH = NaN;
                    invGamPDFHOLDER = NaN(length(x), 1);
                    GamPDFHOLDER = NaN(length(x), 1);
                    chiStat1 = table2struct(table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]));
                    chiStat1.chi2stat = NaN;
                    chiStat1.df = NaN;
                    chiStat2 = table2struct(table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]));
                    chiStat2.chi2stat = NaN;
                    chiStat2.df = NaN;
                end
            end
            if mod(i-skippedIterations, 10) == 0
            disp("Fit " + num2str(i) + " samplings")
            end
        end
    end

    if skipIteration == 1
        skippedIterations = skippedIterations + 1;
    else

    %Run chi2gof on the fit
    %chi2gof on MLN
    MLNVEC.x = x;
    MLNVEC.px = SR_MixLogNorm1RunHOLDER(:,2);
    MLNVEC.numParams = 6;

    [logMLNVEC.x, logMLNVEC.px] = px_to_pfx(MLNVEC.x, MLNVEC.px, @log);
    logMLNVEC.numParams = MLNVEC.numParams;
    %[~, ~, MLNchiStat] = chi2gof_vsMLN(gmfit, log(data), numSRcalcs, fitS);\
    %[~,~,MLNchiStat] = chi2_dataVSpdfVEC(log(data), numSRcalcs, 20, MLNVEC, fitS);

    %chi2gof on inverse gamma
    InvGamVEC.x = x;
    InvGamVEC.px = invGamPDFHOLDER;
    InvGamVEC.numParams = 2;

    [logInvGamVEC.x, logInvGamVEC.px] = px_to_pfx(InvGamVEC.x, InvGamVEC.px, @log);
    logInvGamVEC.numParams = 2;

    [h1(i), p1(i), chiStat1, h2(i), p2(i), chiStat2] = chi2_dataVStwopdfVECs(log(data), numSRcalcs, 20, logMLNVEC, logInvGamVEC, fitS);


    end
    %Save this round's data
    MLNfits{i} = gmfitH;
    MLNs(:, i) = SR_MixLogNorm1RunHOLDER(:,2);
    IGfits{i} = IGfitH;
    invGams(:,i) = invGamPDFHOLDER;
    %GamPDF(:,i) = GamPDFHOLDER;
    weightedC{i} = data;
    numCpairs(i) = numSRcalcs;
    MLNchiStatT = [MLNchiStatT; struct2table(chiStat1, 'AsArray',true)];%#ok<AGROW>
    InvGamChiStatT = [InvGamChiStatT; struct2table(chiStat2, 'AsArray', true)];%#ok<AGROW>
end
outS.weightedC = weightedC;
outS.numCpairs = numCpairs;
outS.MLNPDFs = MLNs;
outS.invGamPDFs = invGams;
outS.MLNfits = MLNfits;
%outS.GamPDFs = GamPDF;
MLNchiStatT1 = addvars(MLNchiStatT, h1', p1', 'Before',"chi2stat", 'NewVariableNames', ["h", "p"]);
outS.MLNchiStats = MLNchiStatT1;
InvGamChiStatT1 = addvars(InvGamChiStatT, h2', p2', 'Before',"chi2stat", 'NewVariableNames', ["h", "p"]);
outS.InvGamChiStats = InvGamChiStatT1;
