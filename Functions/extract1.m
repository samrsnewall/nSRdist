function[LabIDs, incDepths, excLabIDs, excDepths] = extract1(rawdataMSPF, chosenCoresLog)

LabIDs      = table2cell(rawdataMSPF(chosenCoresLog, "LabIDs"));            %take list of LabIDs relating to MSPF dates of each core
incDepths   = table2cell(rawdataMSPF(chosenCoresLog, "IncludeDepths"));     % take list of depths (useful if no labels)
excLabIDs   = table2cell(rawdataMSPF(chosenCoresLog, "excludeLabIDs"));     %take list of manually removed dates for each core
excDepths   = table2cell(rawdataMSPF(chosenCoresLog, "excludeDepth"));      %take list of manually removed dates for each core (useful if no labels)
end
