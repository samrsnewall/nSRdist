function [core_transnums, core_CSE2x, nSRcounts, agediffs] = oneCoreTM(corename, scenarios, LabIDs, incDepths, excLabIDs, excDepths)
%% Return empties if core not used
%Return Nans if core not used
if isempty(scenarios) 
    core_transnums = nan(3,3);
    core_CSE2x = nan(3,1);
    nSRcounts = [];
    agediffs = [];
    return
end

%% Read in Radiocarbon Data
%Read in some radiocarbon data from a net cdf file
WA_path = "/Applications/PaleoDataView/WA_Foraminiferal_Isotopes_2022";
fnm = fullfile(WA_path, "Age/", corename + ".age");
%Read in radiocarbon data from the core
depth_m = ncread(fnm, "Depth"); %(meters)
depth_cm = depth_m.*100; %convert to cm
age = ncread(fnm, "Age dated"); %(14C kyrs BP)
error = ncread (fnm, "Age +Error"); %(14C kyrs BP)
label = ncread(fnm, "Label"); %(Lab ID)
label = string(label);

%% Filter ages
[age, depth_cm, error, label, emptybreak1, emptybreak2] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths);

if emptybreak1 == 1 || emptybreak2 == 1
    core_transnums = nan(3,3);
    core_CSE2x = nan(3,1);
    nSRcounts = [];
    agediffs = [];
    return
end
    
%% Create fake labels if needed (think this will not be needed in this function anymore?)
% If there are no LabIDs, create temporary LabIDs so that scenarios can be made
if sum(contains(label, "NaN"))~=0
    disp("Creating Fake Labels in OneCoreTM for " + corename)
    for i = find(contains(label, "NaN"))
        label(i) = "FakeLabel" + num2str(i);
    end
end

%% Check which scenarios to analyse
%Set up empty vectors
validScenariosBool = ones(length(scenarios),1);

% Check which scenarios are invalid (due to reversals etc).
for i_sce = 1:length(scenarios)
    if scenarios{i_sce} == "invalid"
        validScenariosBool(i_sce) = 0;
    end
end
%% Loop through scenarios using random sampling SR calculations

%Initiate empty variables
scenarioTM = zeros(3,3,length(scenarios));
sce_transnums = zeros(3,3,length(scenarios));
sce_CSE2x = zeros(3,1,length(scenarios));

% Loop through scenarios to calculate TM
for i_sce = 1:length(scenarios)

    %Create empties if scenario is invalid
    if scenarios{i_sce} == "invalid"
        scenarioTM(:,:,i_sce) = NaN;
        continue
    end
    
    %Find the dates in the given scenario
    date_bool = contains(string(label), scenarios{i_sce});
    date_is = find(date_bool == 1);

    %% Perform Calibrations
    %Use MatCal to calibrate each age, storing the probabilities in vector
    %ageprob (note the AGE that each prob is relating to can be found by using
    %the index of that probability -1).
    dep_is = depth_cm(date_is); %Depths of each age
    ageprob = zeros(55001, length(date_is));
    for i = 1:length(date_is)
        [~,~,holder,~] = matcal(age(date_is(i))*1000, error(date_is(i))*1000,  'Marine20', 'CalBP','reserr', 0, 'plot', 0);
        ageprob(:,i) = holder(:,2);
    end
    clear holder %Gets rid of variable holder

    %% Set up years vector and reduce size of calibrated ages (by making NaN)
    m20_years = 0:55000;
    %m20_kyrs = m20_years./1000;

    %% Run a number of random samples, sampling each age with positivity rule
    numruns = 1500;
    transnums_allruns = zeros(3,3,numruns);
    CSE2x_allruns = zeros(3,1,numruns);
    for runN = 1:numruns
        %% Create a potential run of ages
        run_age = NaN(1,length(date_is));
        for i = 1:length(date_is)
            %Set up age probabilities so that every age is older than the
            %age of sample shallower than it
            if i == 1
                poss_ages = m20_years;
                age_probs = ageprob(:,i);
            else 
                min_age = run_age(i-1);
                poss_ages = m20_years(m20_years>min_age);
                age_probs = ageprob(:,i);
                age_probs = age_probs(m20_years>min_age);
                age_probs = age_probs./sum(age_probs);
            end
            run_age(i) = randsample(poss_ages, 1, true, age_probs);
        end

        %% Calculate mean SR for that potential run
        meanSR_run = (dep_is(end)-dep_is(1))./(run_age(end)-run_age(1));

        %% Calculate the sedrates for each pair of ages
        age_diffs = diff(run_age);
        dep_diffs = diff(dep_is)';
        SRs = dep_diffs./age_diffs;
        normSRs = SRs./meanSR_run;
        %Add normSRs to vector to count them (with their weighting)
        if runN == 1 && i_sce == 1
            nSRcounts = [normSRs; dep_diffs./(sum(validScenariosBool).*numruns)];
            agediffs = age_diffs;
        else
            nSRcounts = cat(2, nSRcounts,[normSRs; dep_diffs./(sum(validScenariosBool).*numruns)]);
            agediffs = cat(2, age_diffs);
        end

        %% Categorise normalised sed rates
        %Categories are Steady, Expansion, Contraction (S, E, C)
        char_categ = '';
        for i = 1:(length(date_is)-1)
            if normSRs(i) >= 0.922 && normSRs(i) <1.085
                char_categ(i) = 'S';
            else
                if normSRs(i)>=1.085 && normSRs(i) < inf
                    char_categ(i) = 'E';
                else
                    char_categ(i) = 'C';
                end
            end
        end

        %% Convert this to a transition matrix
        %Find the number of times each transition occurs
        transnums = NaN(3,3,1);
        for jj = 1:1
            transnums(1,1,jj) = length(strfind(char_categ, "CC")); %Probability of Contraction to Contraction
            transnums(1,2,jj) = length(strfind(char_categ, "CS")); %Probability of Contraction to Steady
            transnums(1,3,jj) = length(strfind(char_categ, "CE")); %Probability of Contraction to Expansion
            transnums(2,1,jj) = length(strfind(char_categ, "SC")); %Probability of Steady to Contraction
            transnums(2,2,jj) = length(strfind(char_categ, "SS")); %Probability of Steady to Steady
            transnums(2,3,jj) = length(strfind(char_categ, "SE")); %Probability of Steady to Expansion
            transnums(3,1,jj) = length(strfind(char_categ, "EC")); %Probability of Expansion to Contraction
            transnums(3,2,jj) = length(strfind(char_categ, "ES")); %Probability of Expansion to Steady
            transnums(3,3,jj) = length(strfind(char_categ, "EE")); %Probability of Expansion to Expansion
        end

        %Find the number of times transitions occur from each state
        numCSE2x = zeros(3,1);
        numCSE2x(1) = length(strfind(char_categ(1:end-1), "C"));
        numCSE2x(2) = length(strfind(char_categ(1:end-1), "S"));
        numCSE2x(3) = length(strfind(char_categ(1:end-1), "E"));
        
        %Store these numbers for all runs
        transnums_allruns(:,:,runN) = transnums;
        CSE2x_allruns(:,:,runN) = numCSE2x;
    end
    %Sum over all runs  
    sce_transnums(:,:,i_sce) = sum(transnums_allruns, 3);
    sce_CSE2x(:,:,i_sce) = sum(CSE2x_allruns, 3);
    scenarioTM(:,:,i_sce) = sum(transnums_allruns, 3)./sum(CSE2x_allruns, 3);
end
core_transnums = sum(sce_transnums, 3)./length(scenarios);
core_CSE2x = sum(sce_CSE2x, 3)./length(scenarios);
end