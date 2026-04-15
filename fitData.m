%% fitData — Stage 2 of the nSRdist pipeline
%
% PURPOSE
%   Loads a calcData output file and fits four candidate probability
%   distributions to the nSR data from each NSR estimation method
%   (BMedian, BSampAR, RSR0, RSR500, RSR1000, and optionally individual-run
%   variants BSampIR, RSR0IR, etc.). Saves all fit results for downstream
%   plotting and analysis.
%
% PIPELINE CONTEXT
%   calcData  →  Results/dataT_*.mat
%   fitData   →  loads that file, fits distributions, saves d_*.mat
%
% DISTRIBUTIONS FIT
%   Results from ARfitdists (pooled-run / "AR" fits):
%     MLN  — 2-component Mixed Log-Normal, fitted by EM (fitgmdist)
%     LN   — Single-component Log-Normal (special case of MLN with 1 component)
%     Gam  — Gamma distribution fitted to nSR values
%     invGam — Inverse Gamma: a Gamma fitted to 1/nSR, then transformed back
%
%   Results from IRfitdists (individual-run / "IR" fits, if fitS.fitIRs = true):
%     Same four distributions fitted to each of S.OneRun.numruns individual
%     random samples (one run per core), giving a distribution of fit parameters
%     whose spread reflects run-to-run variability.
%
% NSR ESTIMATION METHODS ANALYSED
%   BMedian   — Bchron median age at each dated depth (analogous to Lin2014)
%   BSampAR   — All BSamp runs pooled and fitted together (ARfitdists)
%   RSR0AR    — All RSR0 runs pooled (no minimum dt restriction)
%   RSR500AR  — All RSR500 runs pooled (dt >= 500 yr required)
%   RSR1000AR — All RSR1000 runs pooled (dt >= 1000 yr required)
%   BSampIR   — Individual BSamp runs fitted separately (if fitS.fitIRs)
%   RSR0IR    — Individual RSR0 runs (if fitS.fitIRs)
%   RSR500IR  — Individual RSR500 runs (if fitS.fitIRs)
%   RSR1000IR — Individual RSR1000 runs (if fitS.fitIRs)
%
% CORE SUBSET
%   d.S1.chooseLog selects cores that pass all of:
%     meanSR  > fitS.minMeanSR   (excludes very low accumulation rate cores)
%     depths  > fitS.minDepth    (excludes shallow-water cores)
%     |lat|   < fitS.maxLat      (excludes high-latitude cores)
%
% OUTPUT (saved to filename2save)
%   d           — Struct containing all fit results. Fields:
%     d.dataT   — The loaded core data table (from calcData)
%     d.S       — The calcData settings structure
%     d.x       — nSR grid (exp(d.logx)) on which PDFs are evaluated
%     d.logx    — Log-space grid for nSR
%     d.S1      — Results struct for the primary core subset:
%       d.S1.fitS      — Copy of the fitS settings used
%       d.S1.chooseLog — Logical index of cores included in this subset
%       d.S1.numCores  — Number of cores in this subset
%       d.S1.BMedian   — ARfitdists output for Bchron Median
%       d.S1.BSampAR   — ARfitdists output for pooled BSamp
%       d.S1.RSR0AR    — ARfitdists output for pooled RSR0
%       d.S1.RSR500AR  — ARfitdists output for pooled RSR500
%       d.S1.RSR1000AR — ARfitdists output for pooled RSR1000
%       d.S1.BSampIR   — IRfitdists output for individual BSamp runs (if fitS.fitIRs)
%       d.S1.RSR0IR    — IRfitdists output for individual RSR0 runs (if fitS.fitIRs)
%       d.S1.RSR500IR  — IRfitdists output for individual RSR500 runs(if fitS.fitIRs)
%       d.S1.RSR1000IR — IRfitdists output for individual RSR1000 runs(if fitS.fitIRs)
%   BM          — BIGMACS reference distribution (fitted MLN and histogram)
%
% DEPENDENCIES
%   Functions/  — all helper functions (added to path at startup)
%   Results/    — must contain the calcData output file (added to path)

%Add important paths
addpath('Functions')
addpath('Results')

%% Load the calcData results
%filepath =    "Results/dataT_All1_RLGtrue_BchronJun2_3Apr26";
filepath =    "Results/dataT_All1_RLGtrue_Bchron8Apr26";
load(filepath)
d.dataT = dataT;
d.S = S;
d.label = "AllCores, R200M20";
split_filepath = split(filepath, ".mat");
filename2save =  split_filepath(1)+ "_fit9Apr26_ageweight.mat";

%% Set up fitting settings structure (fitS)
% fitS controls all aspects of how distributions are fitted to the nSR data.

% --- Core subset filters (applied to d.dataT) ---
fitS.minMeanSR = 8;                 % Minimum mean SR (cm/kyr); cores below this are excluded
fitS.minDepth = 1000;               % Minimum water depth (m)
fitS.maxLat = 40;                   % Maximum absolute latitude (degrees)
%Note, filters (that may be different) are also applied in calcData

% --- Distribution fitting options ---
fitS.fitDists = true;               % If true, fits all four distributions (MLN, LN, Gam, invGam)
fitS.weighting = "depth";           % How nSR estimates are weighted when fitting: "depth", "age", or "none"
fitS.non_normalized_SR = false;     % If true, fits absolute SR (cm/kyr) rather than nSR
fitS.mlnReps = 3;                   % Number of EM random restarts when fitting the Mixed Log-Normal;
                                    %   more restarts reduce the chance of converging to a local optimum

% --- Age filtering ---
fitS.Lin2014AgeFilter = [500 4000]; % [min max] age difference (yr) allowed between dated pairs
                                    % when fitS.Lin2014AgeFiltering = true; mirrors Lin2014 criterion

% --- Data pre-processing ---
fitS.merge_small_dt = false;         % If true, adjacent nSR bins with dt < 500 yr are merged before fitting
                                    % (via merge_small_dt_nSR) to reduce resolution bias

% --- Pooled-run (AR) fitting parameters ---
fitS.invXbinEdges = 0:0.1:15;                  % Bin edges for inverse-SR histogram (used in invGam fitting)
fitS.resampleData = false;                     % If true, after weighting resamples data back to%   original N (removes inflation from replication)
if fitS.weighting == "depth"
fitS.PooledRuns.weightRepDP = 3;               % Decimal places to round data to in replication process
fitS.PooledRuns.weightRepInflator = 0.5;         % Multiplier on weights before replication (increase or reduce size of weights; 1 = no inflation)
elseif fitS.weighting == "age"        
    fitS.PooledRuns.weightRepDP = 3;
    fitS.PooledRuns.weightRepInflator = 0.001;
else
    fitS.PooledRuns.weightRepDP = 3;
    fitS.OneRun.weightRepInflator = 1;
end


% --- Individual-run (IR) fitting (BSampIR, RSR0IR, etc.) ---
fitS.fitIRs = true; % If true, fits distributions to each run separately (slower; gives variability estimate)
if fitS.fitIRs
    fitS.OneRun.numruns = 1000;  % Number of individual runs to sample and fit
            fitS.OneRun.MLNReps = 3;
    if fitS.weighting == "depth"
        fitS.OneRun.weightRepDP = 3;        % Decimal places for weight rounding in IR fits
        fitS.OneRun.weightRepInflator = 1;  % Weight inflator for IR fits
            % MLN random restarts for IR fits
    elseif fitS.weighting == "age"
        fitS.OneRun.weightRepDP = 3;
        fitS.OneRun.weightRepInflator = 0.001;
    else
        fitS.OneRun.weightRepDP = 1;
        fitS.OneRun.weightRepInflator = 1;
    end
end

% --- Chi-squared goodness-of-fit test (optional diagnostic) ---
fitS.run_chi2gof = false;           % If true, runs chi-squared GOF test on each fit
fitS.chi2binN = 5;                  % Initial number of bins for chi-squared test
fitS.dispChi2 = false;              % If true, prints chi-squared results to the console
fitS.chi2MinCountNum = 5;           % Minimum allowed count per bin (bins are merged if below this)
fitS.chi2MaxCountNum = 50;          % Maximum allowed count per bin (bins are split if above this)
fitS.enforceBinSizeLimits = true;   % If true, enforces the min/max bin count limits above

d.S1.fitS = fitS;

fitS.useParallelComp = false;       % If true, uses parfor in IRfitdists (requires Parallel Computing Toolbox)
%% Load BIGMACS reference distribution
% Reads the BIGMACS (Lin2014-based) nSR log-normal and transition matrix
% for comparison with the n5ew fits. The MLN fitted here to BM.hist is
% used later as a reference when plotting.
BMlognorm = readtable("BIGMACSdata/lognormal.txt");                  % (use x values currently used in BIGMACS)
BM.nSR.x = BMlognorm.Var1';
BM.nSR.px = BMlognorm.Var2;
BM.hist = readmatrix("BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BM.hist = BM.hist(:,4);
BM.TM = readmatrix("BIGMACSdata/transition_parameter.txt");
[MLN_BIGMACS, ~, BM.gmfit] = fitMixLogNorm(BM.hist, BM.nSR.x, 2, 3, 343);
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);

%% Select core subset for fitting
d.S1.chooseLog = d.dataT.meanSR > fitS.minMeanSR & d.dataT.depths > fitS.minDepth & abs(d.dataT.lats) < fitS.maxLat;
%Find number of cores
d.S1.numCores = sum(d.S1.chooseLog);

%% Fit distributions to all NSR estimation methods
% ARfitdists pools all runs across all cores in the chosen subset,
% applies weighting (fitS.weighting), and fits MLN, LN, Gam, and invGam.
% The Lin2014 age filter is applied for Bchron-based methods (BMedian,
% BSampAR) but not for RSRx methods (whose age differences are already
% bounded by the dt_min restriction chosen during calcData).

% nSR grid on which all fitted PDFs are evaluated
d.logx = -6:0.01:6;
d.x = exp(d.logx);

% %Lin2014 Method:
fitS.Lin2014AgeFiltering = true;
d.S1.LinNSR = ARfitdists(d.dataT.LinNSR, d.x, d.S1.chooseLog, fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, 1, fitS);

%%
% BMedian: 
fitS.Lin2014AgeFiltering = true;
fitS.merge_small_dt = false;
[d.S1.BMedian] = ARfitdists(d.dataT.BMedian, d.x, d.S1.chooseLog, fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, 1, fitS);

%%
% BSampAR: 
fitS.Lin2014AgeFiltering = true;
fitS.merge_small_dt = false;
[d.S1.BSampAR] = ARfitdists(d.dataT.BSamp, d.x, d.S1.chooseLog, fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, 1000, fitS);

%%
% RSRx pooled fits: no Lin2014 age filter (dt_min restriction serves the
% same purpose); countDivisor set to number of RSR runs per core
fitS.Lin2014AgeFiltering = false;
[d.S1.RSR0AR]    = ARfitdists(d.dataT.RSR0,     d.x, d.S1.chooseLog, fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, S.numruns, fitS);
[d.S1.RSR500AR]  = ARfitdists(d.dataT.RSR500,  d.x, d.S1.chooseLog, fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, S.numruns, fitS);
[d.S1.RSR1000AR] = ARfitdists(d.dataT.RSR1000, d.x, d.S1.chooseLog, fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, S.numruns, fitS);

%%
if fitS.fitIRs
    %% Individual-run (IR) fits — captures run-to-run variability
    % IRfitdists randomly selects one run per core per iteration and fits
    % distributions to that single-run dataset. Repeating this fitS.OneRun.numruns
    % times gives a distribution of fit parameters whose spread quantifies
    % the uncertainty from using any single age model realisation.
    numruns = fitS.OneRun.numruns;
    rng(2) % Fixed seed for reproducibility
    fitS.Lin2014AgeFiltering = true;
    [d.S1.BSampIR] = IRfitdists(d.dataT.BSamp, d.S1.chooseLog, numruns, d.x, fitS);

    fitS.Lin2014AgeFiltering = false;
    fitS.dispChi2 = false;
    rng(2)
    [d.S1.RSR0IR]    = IRfitdists(d.dataT.RSR0,     d.S1.chooseLog, numruns, d.x, fitS);
    [d.S1.RSR500IR]  = IRfitdists(d.dataT.RSR500,  d.S1.chooseLog, numruns, d.x, fitS);
    [d.S1.RSR1000IR] = IRfitdists(d.dataT.RSR1000, d.S1.chooseLog, numruns, d.x, fitS);

    %% Log-space histogram counts for each individual run (diagnostic)
    % Stored as numruns x numBins matrices; used to visualise the spread
    % of empirical distributions across runs alongside the fitted PDFs.
    bw           = 0.1;        % Bin width in log(nSR) space
    logBinEdges  = -5:bw:5;   % Log-space bin edges
    d.S1.BSampIR.lnSRHistCounts   = NaN(numruns, length(logBinEdges)-1);
    d.S1.RSR0IR.lnSRHistCounts   = NaN(numruns, length(logBinEdges)-1);
    d.S1.RSR500IR.lnSRHistCounts = NaN(numruns, length(logBinEdges)-1);
    d.S1.RSR1000IR.lnSRHistCounts= NaN(numruns, length(logBinEdges)-1);

    for i = 1:numruns
        d.S1.BSampIR.lnSRHistCounts(i,:)   = histcounts(log(d.S1.BSampIR.weightedC{i}),   'BinEdges', logBinEdges) ./ d.S1.fitS.OneRun.weightRepInflator;
        d.S1.RSR0IR.lnSRHistCounts(i,:)    = histcounts(log(d.S1.RSR0IR.weightedC{i}),    'BinEdges', logBinEdges) ./ d.S1.fitS.OneRun.weightRepInflator;
        d.S1.RSR500IR.lnSRHistCounts(i,:)  = histcounts(log(d.S1.RSR500IR.weightedC{i}),  'BinEdges', logBinEdges) ./ d.S1.fitS.OneRun.weightRepInflator;
        d.S1.RSR1000IR.lnSRHistCounts(i,:) = histcounts(log(d.S1.RSR1000IR.weightedC{i}), 'BinEdges', logBinEdges) ./ d.S1.fitS.OneRun.weightRepInflator;
    end
end

%% Save all fit results
% Saves the d struct (containing all fit outputs, dataT, and settings) to
% filename2save. The file is not overwritten if it already exists — rename
% filename2save if you want to rerun and save fresh results.
if isfile(filename2save)
    disp("Warning, not saving file because it already exists. Please rename file")
else
    save(filename2save, "d", "-v7.3")
end