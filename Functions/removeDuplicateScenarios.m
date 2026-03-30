function SS = removeDuplicateScenarios(SS)
% removeDuplicateScenarios  Remove duplicate scenarios from a scenarioStruct.
%
% Two scenarios are considered duplicates if they contain exactly the same
% lab IDs in the same order. The function groups scenarios by length, finds
% unique rows within each group, and returns only the unique ones, with all
% associated fields of SS pruned to match.
%
% INPUT / OUTPUT
%   SS - scenarioStruct with fields: scenarios, CFR, chosenLabels,
%        invSRvals, invSRprobs, meanSR, numdatepairs, ageModes,
%        lengthSed, MSI_byage, MSI_bydepth
%
% See also: oneCoreScenarios, scenariosDealWithReversals

if isempty(SS.scenarios)
    return
end

scenlengths  = cellfun(@length, SS.scenarios);
Uscenlengths = unique(scenlengths);
scenarios2keep = [];
shifter = 0;

for i = 1:length(Uscenlengths)
    scenLog         = scenlengths == Uscenlengths(i);
    sameLengthScens = SS.scenarios(scenLog);
    scenArray       = strings(Uscenlengths(i), sum(scenLog));
    for j = 1:sum(scenLog)
        scenArray(:,j) = sameLengthScens{j};
    end
    [~, scenarios2keepI, ~] = unique(scenArray', 'rows');
    scenarios2keepIs = sort(scenarios2keepI);
    scenarios2keep   = [scenarios2keep; scenarios2keepIs + shifter]; %#ok<AGROW>
    shifter          = shifter + sum(scenLog);
end

SS.scenarios    = SS.scenarios(scenarios2keep);
SS.CFR          = SS.CFR(scenarios2keep);
SS.chosenLabels = SS.chosenLabels(scenarios2keep);
SS.invSRvals    = SS.invSRvals(scenarios2keep);
SS.invSRprobs   = SS.invSRprobs(scenarios2keep);
SS.meanSR       = SS.meanSR(scenarios2keep);
SS.numdatepairs = SS.numdatepairs(scenarios2keep);
SS.ageModes     = SS.ageModes(scenarios2keep);
SS.lengthSed    = SS.lengthSed(scenarios2keep);
SS.MSI_byage    = SS.MSI_byage(scenarios2keep);
SS.MSI_bydepth  = SS.MSI_bydepth(scenarios2keep);

end
