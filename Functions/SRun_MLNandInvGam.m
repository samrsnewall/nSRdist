function[outS] = SRun_MLNandInvGam(nSRcounts, coreLog, numruns, x, fitS)

%Initialise vectors
OneRunDatas = cell(1, numruns);
MLNfits = cell(1, numruns);
MLNs = NaN(length(x), numruns);
logMLNs = NaN(length(x), numruns);
LNfits = cell(1, numruns);
LNs = NaN(length(x), numruns);
logLNs = NaN(length(x), numruns);
IGfits = cell(1, numruns);
invGams = NaN(length(x), numruns);
logInvGams =  NaN(length(x), numruns);
Gfits = cell(1, numruns);
Gams = NaN(length(x), numruns);
logGams =  NaN(length(x), numruns);
%phats  = NaN(2,numruns);
weightedC = cell(1,numruns);
numCpairs = NaN(1,numruns);
sedLength = NaN(1, numruns);
sedTimeSpan = NaN(1,numruns);
%chi2stats = NaN(1,numruns);
h = NaN(4, numruns);
p = NaN(4, numruns);
%skippedIterations = 0;
agediffsBinEdges = 0:100:10000;
agediffsWbC = NaN(length(agediffsBinEdges)-1,numruns);

% Initialize cell arrays for chi-squared statistics
MLNchi2stat = cell(numruns,1);
MLNchiStatTedges = cell(numruns,1);
MLNchiStatTO = cell(numruns,1);
MLNchiStatTE = cell(numruns,1);
InvGamChi2stat = cell(numruns,1);
InvGamChiStatTedges = cell(numruns,1);
InvGamChiStatTO = cell(numruns,1);
InvGamChiStatTE = cell(numruns,1);
GamChi2stat = cell(numruns,1);
GamChiStatTedges = cell(numruns,1);
GamChiStatTO = cell(numruns,1);
GamChiStatTE = cell(numruns,1);
LNChi2stat = cell(numruns,1);
LNChiStatTedges = cell(numruns,1);
LNChiStatTO = cell(numruns,1);
LNChiStatTE = cell(numruns,1);

% Initialize numeric arrays for degrees of freedom (df)
MLNchiStatTdf = NaN(numruns,1);
InvGamChiStatTdf = NaN(numruns,1);
GamChiStatTdf = NaN(numruns,1);
LNChiStatTdf = NaN(numruns,1);

%Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove    = cellfun(@isempty,nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~cores2remove);
runN            = NaN(numruns, length(nSRcountsClean));

%Create log(x) vector
logx = sort(log(x));

%% Take 1 run's worth of data from each core
for i = 1:numruns
    OneRunData1 = [];
    for j = 1:length(nSRcountsClean)
        %Open the core's cell
        nSRcount_opencell = nSRcountsClean{j};

        %Find the NaN's that split runs
        split_index       = find(isnan(nSRcount_opencell(1,:)));
        split_indexStart  = split_index + 1;
        split_indexEnd    = split_index - 1;
        split_indexEnd    = [split_indexEnd(2:end), size(nSRcount_opencell, 2)];
        if i ==1
            runN(:,j) = randperm(length(split_index), numruns);  %%%This allows me to make the scenarios chosen random. The length(split_index) finds out how many times each core was sampled (when sampled previously, the sampling randomly chose scenarios). This chooses random samples (and therefore random scenarios).
        end

        % use the NaN indeces to choose the data of that run for each core
        % and concatenate it into OneRunData
        nSRcount_1core1run = nSRcount_opencell(1:4,split_indexStart(runN(i,j)):split_indexEnd(runN(i, j)));
        if sum(nSRcount_1core1run(1,:) <=0) ~=0
            disp("check this core")
        end
        OneRunData1 = cat(2, OneRunData1, nSRcount_1core1run);
    end
    OneRunDatas{i} = OneRunData1;
end

%% Set up parallel process, if wanted (useful to turn off for debugging)
if fitS.useParallelComp
    if sum(size(gcp("nocreate"))) == 0
        parpool('Processes', 4)
    end
else
    ps = parallel.Settings;
    psACset = ps.Pool.AutoCreate;
    ps.Pool.AutoCreate = false;
end

%% Run loop to fit distributions to nSR data
for i = 1:numruns
    %To quiet some warnings, initialize temporaries within parfor loop;
    totalSedLength = [];
    totalSedAge = [];
    data = [];
    numSRcalcs = [];
    agediffsWC1R = [];
    MixLogNormPDFHOLDER = [];
    MLNfitH = [];
    invGamPDFHOLDER = [];
    IGfitH = [];
    chiStati = [];
    logMLNVEC = [];
    logInvGamVEC = [];

    %Get onerundata for this run
    OneRunData = OneRunDatas{i};

    %set up while loop conditions
    allow1 = 0; %If
    skipIteration = 0;
    while skipIteration == 0 & allow1 == 0
        %Get rid of negative sedimentation rates (ideally, shouldn't be
        %there, but BMode allows them to be there)
        if sum(OneRunData(1,:)<0) ~= 0
            warning("There are negative sed rates!... being removed")
            OneRunData = OneRunData(:, OneRunData(1,:) >= 0);
        end

        %Get rid of SRs of 0 (ideally, shouldn't be there)
        if sum(OneRunData(1,:) == 0) ~= 0
            warning("There are SRs of 0!... being removed.")
            OneRunData = OneRunData(:,OneRunData(1,:) ~= 0);
        end

        OneRunData_before = OneRunData;
        if fitS.merge_small_dt
            [OneRunData, mergeLog] = merge_small_dt_nSR(OneRunData, 500);
        end

        %Apply filtering done by Lin2014 if wanted
        if fitS.Lin2014AgeFiltering
            L2014Log = (OneRunData(4,:) < max(fitS.Lin2014AgeFilter) & OneRunData(4,:) > min(fitS.Lin2014AgeFilter)) | isnan(OneRunData(1,:));
            OneRunData = OneRunData(:,L2014Log);
        end


        %% Weighting Data
        weightRepDP = fitS.OneRun.weightRepDP;
        weightRepInflator = fitS.OneRun.weightRepInflator;

        %Set up nSR counts and weighting vectors without NaNs
        NaN_logi        = ~isnan(OneRunData(1,:));
        inputData       = OneRunData(1,NaN_logi);
        depthDiffs      = OneRunData(3,NaN_logi);
        ageDiffs        = OneRunData(4,NaN_logi);

        %Prepare data according to weighting choice
        if fitS.weighting == "none"
            weights = ones(size(inputData));
            data = inputData;
        else
            if fitS.weighting == "depth" %This indicates whether to weight by depth or not
                weights         = depthDiffs; %Weights by depth, not the mixed weighting of depth and scenarios
            elseif fitS.weighting == "age"
                weights         = ageDiffs; %Weights by age, not the mixed weighting of age and scenarios
            else
                error("Weighting type not recognized, must be 'none', 'depth', or 'age'")
            end
            data            = makeWeightedReplicates(inputData, weights, weightRepDP, weightRepInflator); %Weight by replicating data according to weighting
        end

        %% Calculate some other useful information
        %Find out how much length and time the data spans
        totalSedLength = sum(OneRunData(3,:), "all", "omitnan");
        totalSedAge = sum(OneRunData(4,:), "all", "omitnan");

        %Find number of SR calcs actually involved
        numSRcalcs = numel(inputData);

        %Get weighted bin counts of agediffs
        agediffsWC1R = makeWeightedBinCounts(ageDiffs, weights, agediffsBinEdges);
        %figure; histogram("BinEdges", agediffsBinEdges, "BinCounts", agediffsWC1R)

        %% Resample dataset to include impact of weighting avoid inflation of dataset size
        if fitS.resampleData
            rng(1);
            data = randsample(data, numSRcalcs, false);
        end

        %% Fit distributions

        if fitS.fitDists == true
            %Check that the data is a row vector

            [a,b] = size(data);
            if a>b
                data = data';
            end

            %Fit a mix log normal, inverse gamma, gamma and lognormal to the data. If there is an ill conditioned
            %covariance and an error is called, try again (initial conditions are
            %random so new initial conditions may allow the code to avoid an ill
            %conditioned covariance).
            allow2 = 0;
            errorNum = 0;
            while allow2 == 0 && skipIteration == 0
                try
                    [MixLogNormPDFHOLDER, ~, MLNfitH] = fitMixLogNorm(data, x, 2, fitS.OneRun.MLNReps, numSRcalcs); %If this is taking too long, try reducing the weightRepInflator value
                    [invGamPDFHOLDER, ~,IGfitH] = fitInvGamma(data, x, numSRcalcs);
                    [GamPDFHOLDER, GfitH] = fitGamma(data, x, numSRcalcs);
                    [LNPDFHOLDER, ~, LNfitH] = fitMixLogNorm(data, x, 1, fitS.OneRun.MLNReps, numSRcalcs);

                    %Remove any NaNs in HOLDERS (think this is a bug?)
                    MixLogNormPDFHOLDER(isnan(MixLogNormPDFHOLDER)) = 0;
                    invGamPDFHOLDER(isnan(invGamPDFHOLDER)) = 0;
                    GamPDFHOLDER(isnan(GamPDFHOLDER)) = 0;
                    LNPDFHOLDER(isnan(LNPDFHOLDER)) = 0;

                    allow1 = 1;
                    allow2 = 1;
                    skipIteration = 0;
                catch %If there is an error in the fitMixLogNorm, which commonly happens because of an ill conditioned covariance, which may not occur if different initial conditions are used. If it continues to happen, skip this scenario and choose another scenario.
                    errorNum = errorNum+1;
                    disp("number of failed fittings = " + num2str(errorNum))
                    allow2 = 0;
                    if errorNum >10
                        skipIteration = 1;
                        MLNfitH = NaN;
                        MixLogNormPDFHOLDER = NaN(length(x), 2);
                        IGfitH = NaN;
                        invGamPDFHOLDER = NaN(length(x), 1);
                        %GamPDFHOLDER = NaN(length(x), 1);
                        hi(:) = NaN;
                        pri(:) = NaN;
                        chiStati = cell(1,size(pdfs,1));
                        for ji = 1:size(pdfs, 1)
                            chiStati{ji}.chi2stat = NaN;
                            chiStati{ji}.E = NaN;
                            chiStati{ji}.O = NaN;
                            chiStati{ji}.edges = NaN;
                            chiStati{ji}.df = NaN;
                        end
                        logMLNVEC.px = NaN;
                        logInvGamVEC.px = NaN;
                    end
                end
                if mod(i, 50) == 0
                    disp("Fit " + num2str(i) + " samplings")
                end

            end
        else

            % ===== Dummy outputs when skipping distribution fitting =====
            nx = numel(x);

            % PDF holders with correct shapes
            MixLogNormPDFHOLDER = NaN(nx, 2);  % later uses (:,2)
            invGamPDFHOLDER     = NaN(nx, 1);
            GamPDFHOLDER        = NaN(nx, 1);
            LNPDFHOLDER         = NaN(nx, 2);  % later uses (:,2)

            % Fit “handles” / objects
            MLNfitH = NaN;
            IGfitH = NaN;
            GfitH  = NaN;
            LNfitH = NaN;


            allow1 = 1;
            allow2 = 1;
            skipIteration = 0;
        end
    end


    if skipIteration == 1
        disp("Skipped an iteration in SRun_MLNandInvGam")
    else


        %% Run chi2gof on the fits
        %chi2gof on MLN
        MLNVEC = [];
        MLNVEC.x = x';
        MLNVEC.px = MixLogNormPDFHOLDER(:,2);
        MLNVEC.numParams = 5;
        MLNVEC.pdfName = "2 Component Mixed Log Normal";
        [MLNVEC.mu, MLNVEC.var] = muVarPDFVec(MLNVEC);

        logMLNVEC = [];
        [logMLNVEC.x, logMLNVEC.px] = px_to_pfx(MLNVEC.x, MLNVEC.px, @log);
        logMLNVEC.numParams = MLNVEC.numParams;
        logMLNVEC.pdfName = MLNVEC.pdfName;
        [logMLNVEC.mu, logMLNVEC.var] = muVarPDFVec(logMLNVEC);

        %chi2gof on inverse gamma
        InvGamVEC = [];
        InvGamVEC.x = x';
        InvGamVEC.px = invGamPDFHOLDER;
        InvGamVEC.numParams = 2;
        InvGamVEC.pdfName = "Inverse Gamma";
        [InvGamVEC.mu, InvGamVEC.var] = muVarPDFVec(InvGamVEC);

        logInvGamVEC = [];
        [logInvGamVEC.x, logInvGamVEC.px] = px_to_pfx(InvGamVEC.x, InvGamVEC.px, @log);
        logInvGamVEC.numParams = InvGamVEC.numParams;
        logInvGamVEC.pdfName = InvGamVEC.pdfName;
        [logInvGamVEC.mu, logInvGamVEC.var] = muVarPDFVec(logInvGamVEC);

        %chi2gof on gamma
        GamVEC = [];
        GamVEC.x = x';
        GamVEC.px = GamPDFHOLDER;
        GamVEC.numParams = 2;
        GamVEC.pdfName = "Gamma";
        [GamVEC.mu, GamVEC.var] = muVarPDFVec(GamVEC);

        logGamVEC = [];
        [logGamVEC.x, logGamVEC.px] = px_to_pfx(GamVEC.x, GamVEC.px, @log);
        logGamVEC.numParams = GamVEC.numParams;
        logGamVEC.pdfName = GamVEC.pdfName;
        [logGamVEC.mu, logGamVEC.var] = muVarPDFVec(logGamVEC);

        %chi2gof on log normal
        LNVEC = [];
        LNVEC.x = x';
        LNVEC.px = LNPDFHOLDER(:,2);
        LNVEC.numParams = 2;
        LNVEC.pdfName = "Log Normal";
        [LNVEC.mu, LNVEC.var] = muVarPDFVec(LNVEC);

        logLNVEC = [];
        [logLNVEC.x, logLNVEC.px] = px_to_pfx(LNVEC.x, LNVEC.px, @log);
        logLNVEC.numParams = LNVEC.numParams;
        logLNVEC.pdfName = LNVEC.pdfName;
        [logLNVEC.mu, logLNVEC.var] = muVarPDFVec(logLNVEC);

        if fitS.run_chi2gof

            pdfs = {logLNVEC; logMLNVEC; logGamVEC; logInvGamVEC};

            if i <4
                fitS.dispChi2 = true;
            else
                fitS.dispChi2 = false;
            end
            [hi, pri, chiStati] = chi2_dataVStwopdfVECs(log(data), numSRcalcs, 20, pdfs, fitS);

        else %If fitS.run_chi2gof is false, set up dummy variables.
            % Number of distributions
            nPdf = 4;

            % h(:,i) = hi'; p(:,i) = pri'; so hi, pri should be 1×nPdf
            hi  = nan(1, nPdf);
            pri = nan(1, nPdf);

            % chiStati is a 1×nPdf cell array of structs with the same fields
            % as returned by chi2_dataVStwopdfVECs, but all NaN.
            chiStati = cell(1, nPdf);
            for k = 1:nPdf
                chiStati{k} = struct( ...
                    'chi2stat', NaN, ...
                    'df',       NaN, ...
                    'edges',    NaN, ...
                    'O',        NaN, ...
                    'E',        NaN );
            end
        end



    end
    %Save this round's data
    agediffsWbC(:,i) = agediffsWC1R';
    MLNfits{i} = MLNfitH;
    MLNs(:, i) = MLNVEC.px;
    MLNmu(i) = MLNVEC.mu;
    MLNvar(i) = MLNVEC.var;
    logMLNs(:,i) = logMLNVEC.px;
    logMLNmu(i) = logMLNVEC.mu;
    logMLNvar(i) = logMLNVEC.var;
    LNfits{i} = LNfitH;
    LNs(:, i) = LNVEC.px;
    LNmu(i) = LNVEC.mu;
    LNvar(i) = LNVEC.var;
    logLNs(:,i) = logLNVEC.px;
    logLNmu(i) = logLNVEC.mu;
    logLNvar(i) = logLNVEC.var;
    IGfits{i} = IGfitH;
    invGams(:,i) = InvGamVEC.px;
    invGammu(i) = InvGamVEC.mu;
    invGamvar(i) = InvGamVEC.var;
    logInvGams(:,i) = logInvGamVEC.px;
    loginvGammu(i) = logInvGamVEC.mu;
    loginvGamvar(i) = logInvGamVEC.var;
    Gfits{i} = GfitH;
    Gams(:,i) = GamVEC.px;
    Gammu(i) = GamVEC.mu;
    Gamvar(i) = GamVEC.var;
    logGams(:,i) = logGamVEC.px;
    logGammu(i) = logGamVEC.mu;
    logGamvar(i) = logGamVEC.var;
    weightedC{i} = data;
    numCpairs(i) = numSRcalcs;
    sedLength(i) = totalSedLength;
    sedTimeSpan(i) = totalSedAge;
    h(:,i) = hi';
    p(:,i) = pri';
    MLNchi2stat{i} = chiStati{1}.chi2stat;
    MLNchiStatTdf(i) = chiStati{1}.df;
    MLNchiStatTedges{i} = chiStati{1}.edges;
    MLNchiStatTO{i} = chiStati{1}.O;
    MLNchiStatTE{i} = chiStati{1}.E;
    InvGamChi2stat{i} = chiStati{2}.chi2stat;
    InvGamChiStatTdf(i) = chiStati{2}.df;
    InvGamChiStatTedges{i} = chiStati{2}.edges;
    InvGamChiStatTO{i} = chiStati{2}.O;
    InvGamChiStatTE{i} = chiStati{2}.E;
    GamChi2stat{i} = chiStati{3}.chi2stat;
    GamChiStatTdf(i) = chiStati{3}.df;
    GamChiStatTedges{i} = chiStati{3}.edges;
    GamChiStatTO{i} = chiStati{3}.O;
    GamChiStatTE{i} = chiStati{3}.E;
    LNChi2stat{i} = chiStati{4}.chi2stat;
    LNChiStatTdf(i) = chiStati{4}.df;
    LNChiStatTedges{i} = chiStati{4}.edges;
    LNChiStatTO{i} = chiStati{4}.O;
    LNChiStatTE{i} = chiStati{4}.E;
end

% Return ps.Pool.AutoCreate settings to original settings
if fitS.useParallelComp
else
    ps.Pool.AutoCreate = psACset;
end

outS.OneRunDatas = OneRunDatas;
outS.agediffsWbC = agediffsWbC;
outS.weightedC = weightedC;
outS.numCpairs = numCpairs;
outS.sedLength = sedLength;
outS.sedTimeSpan = sedTimeSpan;
outS.MLN.nSR.px = MLNs;
outS.MLN.nSR.x = x;
outS.MLN.nSR.mu = MLNmu;
outS.MLN.nSR.var = MLNvar;
outS.MLN.lnSR.px = logMLNs;
outS.MLN.lnSR.x = logx;
outS.MLN.lnSR.mu = logMLNmu;
outS.MLN.lnSR.var = logMLNvar;
outS.MLN.fits = MLNfits;
outS.invGam.nSR.px = invGams;
outS.invGam.nSR.x = x;
outS.invGam.nSR.mu = invGammu;
outS.invGam.nSR.var = invGamvar;
outS.invGam.lnSR.px = logInvGams;
outS.invGam.lnSR.x = logx;
outS.invGam.lnSR.mu = loginvGammu;
outS.invGam.lnSR.var = loginvGamvar;
outS.invGam.fits = IGfits;
outS.Gam.nSR.px = Gams;
outS.Gam.nSR.x = x;
outS.Gam.nSR.mu = Gammu;
outS.Gam.nSR.var = Gamvar;
outS.Gam.lnSR.px = logGams;
outS.Gam.lnSR.x = logx;
outS.Gam.lnSR.mu = logGammu;
outS.Gam.lnSR.var = logGamvar;
outS.Gam.fits = Gfits;
outS.LN.nSR.px = LNs;
outS.LN.nSR.x = x;
outS.LN.nSR.mu = LNmu;
outS.LN.nSR.var = LNvar;
outS.LN.lnSR.px = logLNs;
outS.LN.lnSR.x = logx;
outS.LN.lnSR.mu = logLNmu;
outS.LN.lnSR.var = logLNvar;
outS.LN.fits = LNfits;

MLNchiStatT1 = table(MLNchi2stat, h(1,:)', p(1,:)', MLNchiStatTdf, MLNchiStatTedges, MLNchiStatTO, MLNchiStatTE,'VariableNames', ["chi2stat","h","p","df","edges","O", "E"]);
outS.MLN.chiStats = MLNchiStatT1;
InvGamChiStatT1 = table(InvGamChi2stat, h(2,:)', p(2,:)', InvGamChiStatTdf, InvGamChiStatTedges, InvGamChiStatTO, InvGamChiStatTE,'VariableNames', ["chi2stat","h","p","df","edges","O", "E"]);
outS.invGam.chiStats = InvGamChiStatT1;
GamChiStatT1 = table(GamChi2stat, h(3,:)', p(3,:)', GamChiStatTdf, GamChiStatTedges, GamChiStatTO, GamChiStatTE,'VariableNames', ["chi2stat","h","p","df","edges","O", "E"]);
outS.Gam.chiStats = GamChiStatT1;
LNChiStatT1 = table(LNChi2stat, h(4,:)', p(4,:)', LNChiStatTdf, LNChiStatTedges, LNChiStatTO, LNChiStatTE,'VariableNames', ["chi2stat","h","p","df","edges","O", "E"]);
outS.LN.chiStats = LNChiStatT1;