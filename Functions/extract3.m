function[LabIDs, incDepths, excLabIDs, excDepths, dataLoc] = extract3(rawdataMSPF, chosenCoresLog, S)

%Initialise cells
numCores    = sum(chosenCoresLog);
LabIDs      = cell(numCores, 1);
incDepths   = cell(numCores, 1);
excLabIDs   = cell(numCores, 1);
excDepths   = cell(numCores, 1);

%Set up all LabIDs to be "all"
for i = 1:numCores
    LabIDs{i}       = "all";
    incDepths{i}    = "";
end

data = rawdataMSPF(chosenCoresLog, :);

%How to access the data
dataLoc = strings(numCores,1);

%If we are only using Lin 2014 data, check whether we want to modify,
%get modifications and leave function
if S.useLin && ~S.usePF
    dataLoc(data.Lin2014 == 1) = "Lin2014";
    if S.modifyLin2014Data
        reversalIDs  = table2cell(data(:, "Lin2014ManualExclude"));
        for i = 1:numCores
            IDvector = [string(reversalIDs(i))];
            IDvector = IDvector(IDvector ~= "");
            IDstring = "";
            for j = 1:length(IDvector)
                if j == 1
                    IDstring = IDvector(j);
                else
                    IDstring = IDstring+ ", " + IDvector(j);
                end
            end
            %IDstring = [IDvector(1) + ", " + IDvector(2) + ", " + IDvector(3) + ", " + IDvector(4)];
            excLabIDs{i} = char(IDstring);
        end
    end
    return
end

%Add info about how to access
WAUseLog = data.WAuse == 1;
dataLoc(WAUseLog) = "WA";

LinKeepLog = data.Lin2014Keep == 1;
dataLoc(LinKeepLog) = "Lin2014";



%Set up excLabIDs to be a combination of all reasons to exclude LabIDs
%(depending on decision to remove large age gaps or not)
reversalIDs         = table2cell(data(:, "ReversalIDs"));
AgeGapIDs           = table2cell(data(:, "AgeGapIDs"));
NonPlanktonicIDs    = table2cell(data(:, "NonPlanktonicIDs"));
MiscRemovalIDs      = table2cell(data(:, "MiscRemovalIDs"));
Lin2014Removals     = table2cell(data(:, "Lin2014ManualExclude"));

%Set up excDepths to be a combination of all reasons to exclude Depths
ReversalDepths = table2cell(data(:, "ReversalDepths"));
AgeGapDepths = table2cell(data(:, "AgeGapDepths"));

for i = 1:numCores
    if WAUseLog(i) == 1 && LinKeepLog(i) ~= 1
        if S.removeLargeGaps
            IDvector = [string(reversalIDs{i}), string(AgeGapIDs{i}), string(NonPlanktonicIDs{i}), string(MiscRemovalIDs{i})];
        else
            IDvector = [string(reversalIDs{i}), string(NonPlanktonicIDs{i}), string(MiscRemovalIDs{i})];
        end

        IDvector = IDvector(IDvector ~= "");
        IDstring = "";
        for j = 1:length(IDvector)
            if j == 1
                IDstring = IDvector(j);
            else
                IDstring = IDstring+ ", " + IDvector(j);
            end
        end
        %IDstring = [IDvector(1) + ", " + IDvector(2) + ", " + IDvector(3) + ", " + IDvector(4)];
        excLabIDs{i} = char(IDstring);
        if isnan(ReversalDepths{i})
            ReversalDepths{i} = "";
        end

        %Check settings for choice about removing large gaps or not
        if S.removeLargeGaps
            DepthVector = [string(ReversalDepths{i}), string(AgeGapDepths{i})];
        else
            DepthVector = [string(ReversalDepths{i})];
        end

        DepthVector = DepthVector(DepthVector ~= "");
        DepthString = "";
        for j = 1:length(DepthVector)
            if j == 1
                DepthString = DepthVector(j);
            else
                DepthString = DepthString+ ", " + DepthVector(j);
            end
        end
        excDepths{i} = char(DepthString);

    elseif WAUseLog(i) == 1 && LinKeepLog(i) == 1
        if S.modifyLin2014Data
            IDvector = [string(Lin2014Removals{i})];
            IDvector = IDvector(IDvector ~= "");
            IDstring = "";
            for j = 1:length(IDvector)
                if j == 1
                    IDstring = IDvector(j);
                else
                    IDstring = IDstring+ ", " + IDvector(j);
                end
                excLabIDs{i} = char(IDstring);
            end
        end
    end
end