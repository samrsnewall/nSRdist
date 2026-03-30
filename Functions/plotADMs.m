%function[xxx] = plotADMs(nSRcounts,coresUsed,cores)

%
numCores = length(cores);
testCoreData = nSRcounts{1};
testseparators  = find([isnan(testCoreData(1,:)), 1]); %Find indeces to separate runs
numRuns     = length(separators)-1; %Find number of runs

allDepths = cell(1,numCores);
allAges = cell(1,numCores);
allSRs = cell(1,numCores);


for i = 1:numCores
    coreData        = nSRcounts{i}; %Pick out core
    
    if ~isempty(coreData)
    separators  = find([isnan(coreData(1,:)), 1]); %Find indeces to separate runs
    runDepths   = NaN(60, numRuns); %Initialise matrices
    runAges     = NaN(60, numRuns);
    runSRs      = NaN(60, numRuns);

    for j = 1:(numRuns)
        depthInfo       = coreData(2,separators(j):separators(j+1)-1); %Get the depths from each run
        runDepths1Run   = cumsum(depthInfo);
        ageInfo         = coreData(3,separators(j):separators(j+1)-1); %Get the ages from each run
        runAges1Run     = cumsum(ageInfo);
        runSRs1Run      = [coreData(1, (separators(j)+1):(separators(j+1)-1)), coreData(1, separators(j+1)-1)];

        runDepths(1:length(runDepths1Run), j) = runDepths1Run; %Put the depths into the initialised matrices
        runAges(1:length(runDepths1Run), j) = runAges1Run; %Put the ages...
        runSRs(1:length(runDepths1Run), j) = runSRs1Run;
    end

    %Check, for every row that isn't NaN, what is the mean nSR
    meanNSR{i} = mode(runSRs, 2);

    allDepths{i} = runDepths;
    allAges{i} = runAges;
    allSRs{i} = runSRs;

    % figure(fADMS)
    % subplot(5,11,i)
    % plot(runDepths', runAges', 'Color', [0.2 0.2 0.2 0.1])
    % xlabel("Depth (cm)")
    % ylabel("Age (cal kyr BP)")
    % title(cores{i})

    end


end


%% 
figure()
for i = 1:10     
    subplot(1,11,i)
    stairs(allDepths{i}, allSRs{i}, 'Color', [0.2 0.2 0.2 0.1])
    hold on
    stairs(allDepths{i}(:,1), meanNSR{i}, 'r')
    xlabel("Depth (cm)")
    ylabel("nSR")
    title(cores{i})
    stopper = 1;
end



%end