
%Add important paths
addpath('Functions')
addpath('Results')

filepath =    "dataT_All1_RLGtrue_DS0p05_Dec9_fitDec9_depthweight_400R.mat";
load(filepath)


%% Chi2gof testing of BMode data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
fitS.dispChi2 = true;

[h, p, chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.BMode.weightedC{1}), d.S1.BMode.numCpairs, fitS);
chiStat1RunT_BM = struct2table(chiStat1Run_BM);

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.BMode.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of BMedian data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
fitS.dispChi2 = true;

[h, p, chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.BMedian.weightedC{1}), d.S1.BMedian.numCpairs, fitS);
chiStat1RunT_BM = struct2table(chiStat1Run_BM);

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.BMedian.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of Bchron Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.BSampIR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.BSampIR.weightedC{i}), d.S1.BSampIR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.BSampIR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of RSR0 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New0IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New0IR.weightedC{i}), d.S1.New0IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New0IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New500 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New500IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New500IR.weightedC{i}), d.S1.New500IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New500IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New1000 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New1000IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New1000IR.weightedC{i}), d.S1.New1000IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New1000IR.chiStatTvsBM = chiStat1RunT_BM;

%% Chi2gof testing of New1500 Sampling data to BIGMACS fits
% Test the chi2gof of each fit to BIGMACS fit
h = NaN(numruns,1);
p = NaN(numruns,1);
fitS.dispChi2 = false;
chiStat1RunT_BM = table('Size', [0,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
for i = 1:numruns
    if ~isempty(d.S1.New1500IR.weightedC{i})
        [h(i), p(i), chiStat1Run_BM] = chi2gof_vsMLN(BM.gmfit, log(d.S1.New1500IR.weightedC{i}), d.S1.New1500IR.numCpairs(i), fitS);
        chiStat1RunT_BM = [chiStat1RunT_BM; struct2table(chiStat1Run_BM)]; %#ok<AGROW>
    else
        emptyRow = table('Size', [1,5], 'VariableTypes', ["double", "double", "cell", "cell", "cell"], 'VariableNames',["chi2stat","df","edges","O","E"]);
        chiStat1RunT_BM = [chiStat1RunT_BM; emptyRow]; %#ok<AGROW>
    end
end

%Save results to table and add to structure
chiStat1RunT_BM = addvars(chiStat1RunT_BM, h, p, 'Before', "chi2stat");
d.S1.New1500IR.chiStatTvsBM = chiStat1RunT_BM;