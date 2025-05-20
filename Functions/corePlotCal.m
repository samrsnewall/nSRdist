function[] = corePlotCal(corename, LabIDs, incDepths, excLabIDs, excDepths, dataLoc, S)
%% Read in Radiocarbon Data
%%Read in radiocarbon data from the core
if dataLoc == "WA"
    [age, depth_cm, error, label] = getDataWA(corename, S);
elseif dataLoc == "Lin2014"
    [age, depth_cm, error, label] = getDatatxt(corename, S);
end

%% Filter for plotting
[age, depth, error, ~, ~,~, EM, ManE, NotCCR, ageGapTooHigh] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);

%% Plot all radiocarbon data
% figure
% hold on
% errorbar(depth, age, error, "vertical", 'o', "color", 'k', 'DisplayName', "Accepted")
% errorbar(EM.depth, EM.age, EM.error, "vertical", '.', 'color', 'r', 'DisplayName', 'Material Not Chosen')
% errorbar(ManE.depth, ManE.age, ManE.error, "vertical", 'o', 'color', 'r', 'DisplayName', "Manually Excluded")
% errorbar(NotCCR.depth, NotCCR.age, NotCCR.error, "vertical", 'x', 'color', 'r', 'DisplayName', "Outside CC Range")
% plot(ageGapTooHigh.depth, ageGapTooHigh.age, 'r-')
% %ylim([0 age(end)*1.1])
% set(gca, 'YTickLabel',get(gca,'YTick'))
% xlabel("Depth (cm)")
% ylabel(["Radiocarbon Age","(14C kyr BP)"])
% legend("location", "southeast")
% title(corename)

resAge = 0;
resError = 200; %Same as in multiMatcal currently
medCal = zeros(size(age))';
errCal = zeros(size(medCal));
for i = 1:length(age)
    [calAgeHolder, ~,~,medCal(i)] = matcal(age(i)*1000, error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    midCal(i) = sum(p95(1:2, i)/2);
    errCal(i) = p95(1,i) - p95(2,i);
end

%Find where the difference of the median cal ages are > 5000 yrs;
gap5kyrfinder = abs(diff(medCal)) > 5000;
gap5kyrLog = [gap5kyrfinder, 0] | [0, gap5kyrfinder];
gap5kyrAges = medCal(gap5kyrLog);
gap5kyrDepths = depth(gap5kyrLog);



EM.medCal = [];
EM.errCal = [];
for i = 1:length(EM.age)
    calAgeHolder = matcal(EM.age(i)*1000, EM.error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    EM.p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    EM.medCal(i) = sum(EM.p95(1:2, i)/2);
    EM.errCal(i) = EM.p95(1,i) - EM.p95(2,i);
end

ManE.medCal = [];
ManE.errCal = [];
for i = 1:length(ManE.age)
    calAgeHolder = matcal(ManE.age(i)*1000, ManE.error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    ManE.p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    ManE.medCal(i) = sum(ManE.p95(1:2, i)/2);
    ManE.errCal(i) = ManE.p95(1,i) - ManE.p95(2,i);
end

NotCCR.medCal = [];
NotCCR.errCal = [];
for i = 1:length(NotCCR.age)
    calAgeHolder = matcal(NotCCR.age(i)*1000, NotCCR.error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    NotCCR.p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    NotCCR.medCal(i) = sum(NotCCR.p95(1:2, i)/2);
    NotCCR.errCal(i) = NotCCR.p95(1,i) - NotCCR.p95(2,i);
end

if ~isempty(gap5kyrDepths)
    %There is an age pair that has a large gap here

    figure;
    hold on
    errorbar(depth, medCal, errCal, "vertical", "Marker","none" , "color", 'k', 'DisplayName', "95% CI", "LineStyle","none")
    plot(depth, medCal, "Marker", "o", "Color","k", "DisplayName", "Median Cal Age", "LineStyle","none")
    errorbar(EM.depth, EM.medCal, EM.errCal, "vertical", '.', "color", 'r', 'DisplayName', "Material Not Chosen")
    errorbar(ManE.depth, ManE.medCal, ManE.errCal, "vertical", 'o', "color", 'r', 'DisplayName', "Manually Excluded")
    errorbar(NotCCR.depth, NotCCR.medCal, NotCCR.errCal, "vertical", 'x', "color", 'r', 'DisplayName', "Outside CC Range")
    %upperlim = max([medCal(end), EM.medCal(end), ManE.medCal(end), NotCCR.medCal(end)])
    plot(gap5kyrDepths, ones(length(gap5kyrDepths)), 'r-')
    ylim([0 50000])
    %set(gca, 'YTickLabel',get(gca,'YTick'))
    ytickformat('%,6.4g')
    xlabel("Depth (cm)")
    ylabel(["Calibrated Age","(cal yr BP)"])
    legend("location", "southeast")
    title(corename)
end

end