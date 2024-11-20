function [interp_invSR, iProbs_wdN, meanSR, reversalpairs, numpairs, ageMode, depdiff, agediff, MSI_byage, MSI_bydepth, IDpairs, agediffV] = scenariopdfNorm(depth,date_is,label, ageprob, calAge, IDpairs, agediffV, S, plotfigs)
%% Perform Calibrations
%Use MatCal to calibrate each age, storing the probabilities in vector
%ageprob (note the AGE that each prob is relating to can be found by using
%the index of that probability -1).

%% Set up years vector and reduce size of calibrated ages (by making NaN)
calAgekyrs = calAge./1000;

%% Calculate mean sedimentation rate (oldest - youngest) in cm/y

%data are already sorted by depth
%Take mode age of each date
ageMode = zeros(length(date_is), 1);
for i = 1:length(date_is)
    [~,modeInd] = max(ageprob(:,i));
    ageMode(i) = calAge(modeInd);
end
%Calculate time between these ages
agediff = ageMode(end) - ageMode(1);
%Calculate depth difference between these ages
dep_is = depth(date_is); %Depths of each age
depdiff = dep_is(end)-dep_is(1);
meanSR = depdiff./agediff; % cm/y
%Calculate mean sampling interval by age
MSI_byage = agediff./length(date_is); %y/date
%Calculate mean sampling interval by depth
MSI_bydepth = depdiff./length(date_is); %cm/date

%% Calculate age difference and probability for each pairing of ages - PDF METHOD FUNCTION
%By using a function
numpairs = length(date_is)-1;
label_is = label(date_is);
[agediff_vals, agediff_probsums, IDpairs, agediffV] = agediffcalc(ageprob, calAge, numpairs, label_is, IDpairs, agediffV, S);

%Calculate depth differences between each pairing of ages;
deldep = diff(dep_is);

%% Find which pairs create reversals (reversal if p(-ve age difference) > 0.75 - arbitrary)
%Initiate a vector to store indices of which pair of ages results
%in a reversal
reversalpairs = zeros(1,numpairs);
%For each pairing, if the probability of negative age difference is
%greater than 0.75, classify as reversal
for m = 1:numpairs
    if sum(agediff_probsums{m}(agediff_vals{m} <= 0)) >=S.reversalCriteria
        %disp("Likely Reversal!")
        reversalpairs(m) = 1;
    end
end

%If there are reversal pairs, should break out of function
if sum(reversalpairs) > 0
    interp_invSR = 0;
    iProbs_wdN = 0;
    return
end

%% Impose the restriction that all age differences should be positive
agediff_probsumsPos = cell(1,numpairs);
agediff_valsPos = cell(1,numpairs);
for m = 1:numpairs
    gt0log = agediff_vals{m} > 0;
    agediff_probsumsPos{m} = agediff_probsums{m}(gt0log);
    agediff_valsPos{m} = agediff_vals{m}(gt0log);
    agediff_probsumsPos{m} = agediff_probsumsPos{m}./sum(agediff_probsumsPos{m}, 'all');

end
%% Calculate invSRvalues from agediff_vals by dividing by depth diff
invSR_vals = cell(1,numpairs);
invSR_probsums = cell(1,numpairs);
invSRAUC = NaN(1,numpairs);
for m = 1:numpairs
    invSR_vals{m} = agediff_valsPos{m}./deldep(m); % invSR with units of y/cm
    invSRAUC(m) = trapz(invSR_vals{m}, agediff_probsumsPos{m});
    invSR_probsums{m} = agediff_probsumsPos{m}./invSRAUC(m);
end

%invSR_probsums = agediff_probsumsPos;
%% Plot outputs of pairwise pdfs
if plotfigs ==1
    %Plot all the inverse sed rate pdfs
    pdfCreation = figure("Name", "pdfCreation");
    subplot(4,1,1)
    for n = 1:(numpairs)
        plot(invSR_vals{1,n}./1000, invSR_probsums{1,n}, 'LineWidth', 1)
        hold on
    end
    ylabel("Probability")
    xlabel("Inverse of Sed Rate ky/cm")
end
%% Calculate the mean sedimentation rates

%Find mode ages of top and bottom calibrated ages
[~,topage_ind] = max(ageprob(:,1));
[~,bottomage_ind] = max(ageprob(:,end));
topage_mode = calAgekyrs(topage_ind);
bottomage_mode = calAgekyrs(bottomage_ind);
total_agediff = bottomage_mode - topage_mode;
total_depthdiff = depth(end)-depth(1);
meanSR = total_depthdiff./total_agediff; % cm/kyr

%% Normalise wrt mean SR
%Multiply values of inverse SR by the mean SR (same as dividing inverse SR
%by inverse of mean SR) - *(1/1000) to convert to same units
invSR_normvals = cellfun(@(x) x*(meanSR.*(1/1000)), invSR_vals, 'un', 0);

if plotfigs == 1
    figure(pdfCreation)
    subplot(4,1,2)
    hold on
    for n = 1:(numpairs)
        plot(invSR_normvals{1,n}, invSR_probsums{1,n})
    end
    ylabel("Probability")
    xlabel("Normalised Inverse Sed Rate")
end
%% Interpolate normalised pdfs

if S.pdfMethod

    %First, interpolate each pdf so that they are compatible among the invSR values,
    %i.e each invSR_vals vector has a constant spacing and goes through the
    %integers (i.e. 15.0, 15.1, 15.2, instead of 15.032, 15.132, 15.232 etc).
    %Initiate Cells
    interp_invSR_indy = cell(1,numpairs);
    interp_invSRp_indy = cell(1,numpairs);
    %Define Spacing
    spacing = 0.00005;
    %Set up interpolated invSR vector
    vector_invSRvals = vertcat(invSR_normvals{:});
    lowerbound = floor(min(vector_invSRvals*(1/spacing)))*spacing;
    upperbound = max(vector_invSRvals);
    interp_invSR = lowerbound:spacing:(upperbound + spacing*0.9);

    %For each pair of dates:
    for nd = 1:numpairs
        if isempty(invSR_normvals{nd})
        else
            %Find where on the new invSR vector the invSR result of
            %this pair fits
            min_val = min(invSR_normvals{nd});
            max_val = max(invSR_normvals{nd});
            ind1 = find(interp_invSR>min_val,1,'first');
            ind2 = find(interp_invSR<max_val,1,'last');
            %Interpolate within this space
            interp_invSR_indy{nd} = interp_invSR(ind1:ind2);
            interp_invSRp_indy_hold = interp1(invSR_normvals{nd},invSR_probsums{nd},interp_invSR(ind1:ind2));
            interp_invSRp_indy{nd} = interp_invSRp_indy_hold./(sum(interp_invSRp_indy_hold)); %normalise the interpolated pdfs
        end
    end

    % if plotfigs ==1
    %     %Plot all of these interpolated pdfs
    %     figure(pdfCreation)
    %     subplot(4,1,3)
    %     for ii = 1:numpairs
    %         plot(interp_invSR_indy{ii}(:), interp_invSRp_indy{ii}(:))
    %         hold on
    %     end
    %     xlabel('Inverse Sed Rate Interpolations')
    %     ylabel("probability")
    % end

    %%

    %Find out where each pairs interpolated values fit onto the master vector and
    %align the probabilities so that they can be summed
    mat_iProbs = zeros(length(interp_invSR),numpairs);
    for i = 1:numpairs
        if isempty(interp_invSR_indy{i})
        else
            idx2 = find(abs(interp_invSR - interp_invSR_indy{i}(1)) <=0.00001);
            %idx2 = knnsearch(interp_invSR, interp_invSR_indy{i}(1));
            mat_iProbs(idx2:idx2+length(interp_invSR_indy{i}(:))-1,i) = interp_invSRp_indy{i}(:);
        end
    end

    %% Average with weighting by depth between the radiocarbon ages

    %Apply a weighting by multiplying the probability of the sed rate from each
    %pair by the depth between the two dates
    if S.weighting == true
        mat_iProbs_wd = mat_iProbs.*deldep';
    else
        mat_iProbs_wd = mat_iProbs;
    end

    if plotfigs == 1
        figure(pdfCreation)
        subplot(4,1,3)
        hold on
        plot(interp_invSR, mat_iProbs_wd)
        xlabel("Individual invSR Ratios weighted")
        ylabel("Probability")
    end

    %Sum all the probabilities from each pair pdf
    iProbs_wd = sum(mat_iProbs_wd, 2);

    %Calculate area under the curve and normalise again
    AUC = trapz(interp_invSR, iProbs_wd);
    iProbs_wdN = iProbs_wd./AUC;

    if plotfigs ==1
        figure(pdfCreation)
        subplot(4,1,4)
        hold on
        plot(interp_invSR, iProbs_wdN)
        ylabel("Probability")
        xlabel("Inverse of Sed Rate Ratio (Weighted by depth)")
    end
else
    interp_invSR = [];
    iProbs_wdN = [];
end

end