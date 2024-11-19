function[LabIDs, incDepths, excLabIDs, excDepths] = extract2(rawdataMSPF, chosenCoresLog, S)

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

%Set up excLabIDs to be a combination of all reasons to exclude LabIDs
%(depending on decision to remove large age gaps or not)
reversalIDs         = table2cell(rawdataMSPF(chosenCoresLog, "ReversalIDs"));
AgeGapIDs           = table2cell(rawdataMSPF(chosenCoresLog, "AgeGapIDs"));
NonPlanktonicIDs    = table2cell(rawdataMSPF(chosenCoresLog, "NonPlanktonicIDs"));
MiscRemovalIDs      = table2cell(rawdataMSPF(chosenCoresLog, "MiscRemovalIDs"));
for i = 1:numCores
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
end
%Set up excDepths to be a combination of all reasons to exclude Depths
ReversalDepths = table2cell(rawdataMSPF(chosenCoresLog, "ReversalDepths"));
AgeGapDepths = table2cell(rawdataMSPF(chosenCoresLog, "AgeGapDepths"));
for i = 1:numCores
    
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
end

end