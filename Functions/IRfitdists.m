function[outS] = IRfitdists(nSRcounts, coreLog, numruns, x, fitS)
% IRfitdists  Fit probability distributions to nSR data across
%             multiple Monte Carlo runs (Individual Runs).
%
% For each of numruns iterations, selects one run's worth of nSR data from
% each core, applies optional filtering and weighting, then fits four
% candidate distributions: Mixed Log-Normal (MLN), Log-Normal (LN),
% Inverse Gamma, and Gamma. Optionally runs chi-squared goodness-of-fit
% tests on each fitted distribution.
%
% Contrast with ARfitdists, which pools all runs from all cores into a
% single dataset before fitting. IRfitdists instead samples one run per
% core on each iteration, giving a distribution of fitted parameters
% that reflects sampling uncertainty across Monte Carlo runs.
%
% INPUTS
%   nSRcounts - (cell array) One cell per core, each containing a 3-row
%               nSR matrix (see README "Internal Data Formats")
%   coreLog   - (logical vector) Selects which cores from nSRcounts to use
%   numruns   - (scalar) Number of Monte Carlo fitting iterations
%   x         - (numeric vector) nSR evaluation points for PDF output
%   fitS      - (struct) Fitting settings. Key fields:
%                 .useParallelComp       (logical) Use parallel pool
%                 .fitDists              (logical) Fit distributions;
%                                        false = skip fitting, store NaNs
%                 .weighting             "none", "depth", or "age"
%                 .non_normalized_SR     (logical) Use raw SR instead of nSR
%                 .merge_small_dt        (logical) Merge short age intervals
%                 .Lin2014AgeFiltering   (logical) Apply Lin2014 age filter
%                 .resampleData          (logical) Resample to original size
%                 .run_chi2gof           (logical) Run chi-squared test
%                 .OneRun.weightRepDP         Weighting decimal places
%                 .OneRun.weightRepInflator   Weighting inflator
%                 .OneRun.MLNReps             MLN fitting repetitions
%
% OUTPUT
%   outS - Struct with per-run fitted distributions, PDF values, fit
%          handles, chi-squared statistics, and summary info. Fields:
%            .MLN, .invGam, .Gam, .LN  per-distribution results
%            .OneRunDatas, .weightedC  per-run input data
%            .numCpairs, .sedLength, .sedTimeSpan, .agediffsWbC
%
% See also: fitData, fitMixLogNorm, fitInvGamma, fitGamma, ARfitdists

%% --- Setup

logx             = sort(log(x));
nx               = numel(x);
nPdf             = 4;
agediffsBinEdges = 0:100:10000;

% Clean out any cores that returned empty nSR vectors
nSRcountsChosen = nSRcounts(coreLog);
cores2remove    = cellfun(@isempty, nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~cores2remove);
runN            = NaN(numruns, length(nSRcountsClean));

%% --- Sample one run per core, for each of numruns iterations

OneRunDatas = cell(1, numruns);
for i = 1:numruns
    OneRunData1 = [];
    for j = 1:length(nSRcountsClean)
        nSRcount_opencell = nSRcountsClean{j};

        % Find the NaN delimiters that separate individual runs
        split_index      = find(isnan(nSRcount_opencell(1,:)));
        split_indexStart = split_index + 1;
        split_indexEnd   = [split_index(2:end) - 1, size(nSRcount_opencell, 2)];

        if i == 1
            % Randomise which run is drawn from each core (done once,
            % reused across all numruns iterations)
            runN(:,j) = randperm(length(split_index), numruns);
        end

        nSRcount_1core1run = nSRcount_opencell(1:3, split_indexStart(runN(i,j)):split_indexEnd(runN(i,j)));
        if sum(nSRcount_1core1run(1,:) <= 0) ~= 0
            disp("check this core")
        end
        OneRunData1 = cat(2, OneRunData1, nSRcount_1core1run);
    end
    OneRunDatas{i} = OneRunData1;
end

%% --- Parallel pool setup

if fitS.useParallelComp
    if sum(size(gcp("nocreate"))) == 0
        parpool('Processes', 4)
    end
else
    ps = parallel.Settings;
    psACset = ps.Pool.AutoCreate;
    ps.Pool.AutoCreate = false;
end

%% --- Pre-allocate output struct

outS.OneRunDatas  = OneRunDatas;
outS.agediffsWbC  = NaN(length(agediffsBinEdges)-1, numruns);
outS.weightedC    = cell(1, numruns);
outS.numCpairs    = NaN(1, numruns);
outS.sedLength    = NaN(1, numruns);
outS.sedTimeSpan  = NaN(1, numruns);

outS.MLN.nSR.x    = x;    outS.MLN.nSR.px   = NaN(nx, numruns);
outS.MLN.nSR.mu   = NaN(1, numruns);  outS.MLN.nSR.var  = NaN(1, numruns);
outS.MLN.lnSR.x   = logx; outS.MLN.lnSR.px  = NaN(nx, numruns);
outS.MLN.lnSR.mu  = NaN(1, numruns);  outS.MLN.lnSR.var = NaN(1, numruns);
outS.MLN.fits     = cell(1, numruns);

outS.invGam.nSR.x    = x;    outS.invGam.nSR.px   = NaN(nx, numruns);
outS.invGam.nSR.mu   = NaN(1, numruns);  outS.invGam.nSR.var  = NaN(1, numruns);
outS.invGam.lnSR.x   = logx; outS.invGam.lnSR.px  = NaN(nx, numruns);
outS.invGam.lnSR.mu  = NaN(1, numruns);  outS.invGam.lnSR.var = NaN(1, numruns);
outS.invGam.fits     = cell(1, numruns);

outS.Gam.nSR.x    = x;    outS.Gam.nSR.px   = NaN(nx, numruns);
outS.Gam.nSR.mu   = NaN(1, numruns);  outS.Gam.nSR.var  = NaN(1, numruns);
outS.Gam.lnSR.x   = logx; outS.Gam.lnSR.px  = NaN(nx, numruns);
outS.Gam.lnSR.mu  = NaN(1, numruns);  outS.Gam.lnSR.var = NaN(1, numruns);
outS.Gam.fits     = cell(1, numruns);

outS.LN.nSR.x    = x;    outS.LN.nSR.px   = NaN(nx, numruns);
outS.LN.nSR.mu   = NaN(1, numruns);  outS.LN.nSR.var  = NaN(1, numruns);
outS.LN.lnSR.x   = logx; outS.LN.lnSR.px  = NaN(nx, numruns);
outS.LN.lnSR.mu  = NaN(1, numruns);  outS.LN.lnSR.var = NaN(1, numruns);
outS.LN.fits     = cell(1, numruns);

h           = NaN(4, numruns);
p           = NaN(4, numruns);
allChiStati = cell(numruns, 1);

%% --- Main fitting loop

for i = 1:numruns

    %--- 1. Prepare data

    OneRunData = OneRunDatas{i};

    % Remove negative SRs (ideally shouldn't exist, but BMode can produce them)
    if sum(OneRunData(1,:) < 0) ~= 0
        warning("There are negative sed rates!... being removed")
        OneRunData = OneRunData(:, OneRunData(1,:) >= 0);
    end

    % Remove zero SRs
    if sum(OneRunData(1,:) == 0) ~= 0
        warning("There are SRs of 0!... being removed.")
        OneRunData = OneRunData(:, OneRunData(1,:) ~= 0);
    end

    % Optionally use absolute SR instead of normalised SR
    if fitS.non_normalized_SR
        SRs = OneRunData(2,:) ./ OneRunData(3,:);
        SRs(isnan(OneRunData(1,:))) = NaN;
        OneRunData(1,:) = SRs;
    end

    % Merge short age intervals if wanted
    if fitS.merge_small_dt
        [OneRunData, ~] = merge_small_dt_nSR(OneRunData, 500);
    end

    % Apply Lin2014 age filter if wanted
    if fitS.Lin2014AgeFiltering
        L2014Log = (OneRunData(3,:) < max(fitS.Lin2014AgeFilter) & ...
                    OneRunData(3,:) > min(fitS.Lin2014AgeFilter)) | ...
                    isnan(OneRunData(1,:));
        OneRunData = OneRunData(:, L2014Log);
    end

    % Extract non-NaN data and set up weights
    NaN_logi   = ~isnan(OneRunData(1,:));
    inputData  = OneRunData(1, NaN_logi);
    depthDiffs = OneRunData(2, NaN_logi);
    ageDiffs   = OneRunData(3, NaN_logi);

    weightRepDP       = fitS.OneRun.weightRepDP;
    weightRepInflator = fitS.OneRun.weightRepInflator;

    if fitS.weighting == "none"
        weights = ones(size(inputData));
        data    = inputData;
    elseif fitS.weighting == "depth"
        weights = depthDiffs;
        data    = makeWeightedReplicates(inputData, weights, weightRepDP, weightRepInflator);
    elseif fitS.weighting == "age"
        weights = ageDiffs;
        data    = makeWeightedReplicates(inputData, weights, weightRepDP, weightRepInflator);
    else
        error("Weighting type not recognized, must be 'none', 'depth', or 'age'")
    end

    % Compute run-level summary statistics
    totalSedLength = sum(OneRunData(2,:), "all", "omitnan");
    totalSedAge    = sum(OneRunData(3,:), "all", "omitnan");
    numSRcalcs     = numel(inputData);
    agediffsWC1R   = makeWeightedBinCounts(ageDiffs, weights, agediffsBinEdges);

    % Optionally resample to original size (removes weighting inflation)
    if fitS.resampleData
        rng(1);
        data = randsample(data, numSRcalcs, false);
    end

    %--- 2. Fit distributions (with retry on ill-conditioned covariance)

    skipIteration = false; %This will change to true if 10 iterations fail, exiting while loop

    if fitS.fitDists
        % Ensure data is a row vector for fitting functions
        if size(data, 1) > size(data, 2)
            data = data';
        end

        allow2   = false;
        errorNum = 0;
        while ~allow2 && ~skipIteration
            try
                [MixLogNormPDFHOLDER, ~, MLNfitH] = fitMixLogNorm(data, x, 2, fitS.OneRun.MLNReps, numSRcalcs);
                [invGamPDFHOLDER,     ~, IGfitH]  = fitInvGamma(data, x, numSRcalcs);
                [GamPDFHOLDER,           GfitH]   = fitGamma(data, x, numSRcalcs);
                [LNPDFHOLDER,         ~, LNfitH]  = fitMixLogNorm(data, x, 1, fitS.OneRun.MLNReps, numSRcalcs);

                % Replace any NaNs in PDF outputs (think this is a bug in fitters?)
                MixLogNormPDFHOLDER(isnan(MixLogNormPDFHOLDER)) = 0;
                invGamPDFHOLDER(isnan(invGamPDFHOLDER))         = 0;
                GamPDFHOLDER(isnan(GamPDFHOLDER))               = 0;
                LNPDFHOLDER(isnan(LNPDFHOLDER))                 = 0;

                allow2 = true;
            catch
                % fitMixLogNorm can fail with an ill-conditioned covariance;
                % retrying with new random initial conditions often succeeds.
                % If it fails more than 10 times, skip this iteration.
                errorNum = errorNum + 1;
                disp("number of failed fittings = " + num2str(errorNum))
                if errorNum > 10
                    skipIteration = true;
                    [MLNfitH, IGfitH, GfitH, LNfitH, ...
                     MixLogNormPDFHOLDER, invGamPDFHOLDER, GamPDFHOLDER, LNPDFHOLDER] = makeDummyFitResult(nx);
                end
            end
            if mod(i, 50) == 0
                disp("Fit " + num2str(i) + " samplings")
            end
        end
    else
        % Distribution fitting disabled — fill with NaN placeholders
        [MLNfitH, IGfitH, GfitH, LNfitH, ...
         MixLogNormPDFHOLDER, invGamPDFHOLDER, GamPDFHOLDER, LNPDFHOLDER] = makeDummyFitResult(nx);
    end

    if skipIteration
        disp("Skipped an iteration in IRfitdists")
    end

    %--- 3. Build PDF VEC structs
    % buildVEC computes mu and var from the PDF; works correctly with NaN
    % PDFs (from skipped/disabled fitting), producing NaN mu and var.
    MLNVEC       = buildVEC(x, MixLogNormPDFHOLDER(:,2), 5, "2 Component Mixed Log Normal");
    logMLNVEC    = buildLogVEC(MLNVEC);
    InvGamVEC    = buildVEC(x, invGamPDFHOLDER,           2, "Inverse Gamma");
    logInvGamVEC = buildLogVEC(InvGamVEC);
    GamVEC       = buildVEC(x, GamPDFHOLDER,              2, "Gamma");
    logGamVEC    = buildLogVEC(GamVEC);
    LNVEC        = buildVEC(x, LNPDFHOLDER(:,2),          2, "Log Normal");
    logLNVEC     = buildLogVEC(LNVEC);

    %--- 4. Chi-squared goodness-of-fit test

    if ~skipIteration && fitS.run_chi2gof
        pdfs = {logLNVEC; logMLNVEC; logGamVEC; logInvGamVEC};
        fitS.dispChi2 = (i < 4);
        [hi, pri, chiStati_run] = chi2_dataVStwopdfVECs(log(data), numSRcalcs, 20, pdfs, fitS);
    else
        [hi, pri, chiStati_run] = makeDummyChi2(nPdf);
    end

    %--- 5. Store results directly in outS

    outS.agediffsWbC(:,i)    = agediffsWC1R';
    outS.weightedC{i}        = data;
    outS.numCpairs(i)        = numSRcalcs;
    outS.sedLength(i)        = totalSedLength;
    outS.sedTimeSpan(i)      = totalSedAge;

    outS.MLN.fits{i}         = MLNfitH;
    outS.MLN.nSR.px(:,i)     = MLNVEC.px;
    outS.MLN.nSR.mu(i)       = MLNVEC.mu;
    outS.MLN.nSR.var(i)      = MLNVEC.var;
    outS.MLN.lnSR.px(:,i)    = logMLNVEC.px;
    outS.MLN.lnSR.mu(i)      = logMLNVEC.mu;
    outS.MLN.lnSR.var(i)     = logMLNVEC.var;

    outS.invGam.fits{i}      = IGfitH;
    outS.invGam.nSR.px(:,i)  = InvGamVEC.px;
    outS.invGam.nSR.mu(i)    = InvGamVEC.mu;
    outS.invGam.nSR.var(i)   = InvGamVEC.var;
    outS.invGam.lnSR.px(:,i) = logInvGamVEC.px;
    outS.invGam.lnSR.mu(i)   = logInvGamVEC.mu;
    outS.invGam.lnSR.var(i)  = logInvGamVEC.var;

    outS.Gam.fits{i}         = GfitH;
    outS.Gam.nSR.px(:,i)     = GamVEC.px;
    outS.Gam.nSR.mu(i)       = GamVEC.mu;
    outS.Gam.nSR.var(i)      = GamVEC.var;
    outS.Gam.lnSR.px(:,i)    = logGamVEC.px;
    outS.Gam.lnSR.mu(i)      = logGamVEC.mu;
    outS.Gam.lnSR.var(i)     = logGamVEC.var;

    outS.LN.fits{i}          = LNfitH;
    outS.LN.nSR.px(:,i)      = LNVEC.px;
    outS.LN.nSR.mu(i)        = LNVEC.mu;
    outS.LN.nSR.var(i)       = LNVEC.var;
    outS.LN.lnSR.px(:,i)     = logLNVEC.px;
    outS.LN.lnSR.mu(i)       = logLNVEC.mu;
    outS.LN.lnSR.var(i)      = logLNVEC.var;

    h(:,i)         = hi';
    p(:,i)         = pri';
    allChiStati{i} = chiStati_run;

end

%% --- Restore parallel pool settings

if ~fitS.useParallelComp
    ps.Pool.AutoCreate = psACset;
end

%% --- Assemble chi-squared statistics tables

outS.MLN.chiStats    = buildChi2Table(allChiStati, h(1,:)', p(1,:)', 1);
outS.invGam.chiStats = buildChi2Table(allChiStati, h(2,:)', p(2,:)', 2);
outS.Gam.chiStats    = buildChi2Table(allChiStati, h(3,:)', p(3,:)', 3);
outS.LN.chiStats     = buildChi2Table(allChiStati, h(4,:)', p(4,:)', 4);

end

%% =========================================================================
%  Local helper functions
%% =========================================================================

function [MLNfitH, IGfitH, GfitH, LNfitH, MixLogNormPDF, invGamPDF, GamPDF, LNPDF] = makeDummyFitResult(nx)
% Returns NaN-filled fit handles and PDF arrays for use when fitting is
% skipped (fitDists = false) or has permanently failed (errorNum > 10).
MLNfitH       = struct('NumParams', NaN, 'mu', NaN, 'Sigma', NaN, ...
                       'ComponentProportion', NaN, 'NegativeLogLikelihood', NaN, ...
                       'BIC', NaN, 'BICtaeheefix', NaN);
IGfitH        = struct('NumParams', NaN, 'alpha', NaN, 'beta', NaN, ...
                       'NegativeLogLikelihood', NaN, 'BIC', NaN, 'BICtaeheefix', NaN);
GfitH         = struct('NumParams', NaN, 'alpha', NaN, 'beta', NaN, ...
                       'NegativeLogLikelihood', NaN, 'BIC', NaN, 'BICtaeheefix', NaN);
LNfitH        = struct('NumParams', NaN, 'mu', NaN, 'Sigma', NaN, ...
                       'NegativeLogLikelihood', NaN, 'BIC', NaN, 'BICtaeheefix', NaN);
MixLogNormPDF = NaN(nx, 2);
invGamPDF     = NaN(nx, 1);
GamPDF        = NaN(nx, 1);
LNPDF         = NaN(nx, 2);
end

function vec = buildVEC(x, px, numParams, pdfName)
% Constructs a PDF VEC struct and computes its mean and variance.
% Works with NaN px values (from skipped fitting), returning NaN mu/var.
vec.x         = x';
vec.px        = px;
vec.numParams = numParams;
vec.pdfName   = pdfName;
[vec.mu, vec.var] = muVarPDFVec(vec);
end

function logvec = buildLogVEC(vec)
% Transforms a VEC struct from nSR space into log(nSR) space.
[logvec.x, logvec.px] = px_to_pfx(vec.x, vec.px, @log);
logvec.numParams      = vec.numParams;
logvec.pdfName        = vec.pdfName;
[logvec.mu, logvec.var] = muVarPDFVec(logvec);
end

function [hi, pri, chiStati] = makeDummyChi2(nPdf)
% Returns NaN-filled chi-squared results for skipped or disabled tests.
hi       = nan(1, nPdf);
pri      = nan(1, nPdf);
chiStati = cell(1, nPdf);
for k = 1:nPdf
    chiStati{k} = struct('chi2stat', NaN, 'df', NaN, 'edges', NaN, 'O', NaN, 'E', NaN);
end
end

function T = buildChi2Table(allChiStati, h, p, distIdx)
% Assembles a chi-squared statistics table for one distribution across
% all runs. distIdx selects which distribution's results to extract from
% each run's chiStati cell (order matches the pdfs cell in the main loop).
numruns  = length(allChiStati);
chi2stat = cell(numruns, 1);
df       = NaN(numruns, 1);
edges    = cell(numruns, 1);
O        = cell(numruns, 1);
E        = cell(numruns, 1);
for i = 1:numruns
    chi2stat{i} = allChiStati{i}{distIdx}.chi2stat;
    df(i)       = allChiStati{i}{distIdx}.df;
    edges{i}    = allChiStati{i}{distIdx}.edges;
    O{i}        = allChiStati{i}{distIdx}.O;
    E{i}        = allChiStati{i}{distIdx}.E;
end
T = table(chi2stat, h, p, df, edges, O, E, ...
    'VariableNames', ["chi2stat","h","p","df","edges","O","E"]);
end
