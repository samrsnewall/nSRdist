function [nSR_out, mergeLog] = merge_small_dt_nSR(nSR, dt_min)
%MERGE_SMALL_DT_NSR Merge adjacent bins in a 3xN nSR matrix until all dt>=dt_min,
% respecting core breaks marked by NaN in row 1.
%
% nSR rows:
%   row1: ratio normalized dx/dt (NOT additive) -> recomputed as (row2/row3)/normalizer after merges
%   row2: dx        (additive)
%   row3: dt        (additive)  [threshold applied here]
%
% Core breaks:
%   isnan(nSR(1,i)) == true marks a core header column (not a measurement bin).
%   Header columns are preserved and never merged with neighbors.
%   Measurement columns are those with ~isnan(nSR(1,i)).
%
% Output:
%   nSR_out: 4xM matrix with the same structure, headers preserved, row1 updated
%   mergeLog: struct array describing merges (core index and details)

    if nargin < 2
        error('Usage: [nSR_out, mergeLog] = merge_small_dt_nSR(nSR, dt_min)');
    end
    if size(nSR,1) ~= 3
        error('nSR must be a 3xN matrix.');
    end
    if ~isscalar(dt_min) || dt_min <= 0
        error('dt_min must be a positive scalar.');
    end

    header = isnan(nSR(1,:));
    if ~any(header)
        error('No core breaks found: expected at least one NaN in row 1.');
    end

    % Identify core blocks using header columns as starts.
    coreStarts = find(header);

    nSR_out = zeros(3,0);
    mergeLog = struct('core', {}, 'j_in_core', {}, 'direction', {}, 'dt_before', {}, 'dt_after', {});
    logk = 0;

    for c = 1:numel(coreStarts)
        s = coreStarts(c);
        if c < numel(coreStarts)
            e = coreStarts(c+1) - 1;
        else
            e = size(nSR,2);
        end

        % Append header column unchanged (preserve all 4 rows as provided)
        nSR_out(:,end+1) = nSR(:,s);

        % Measurement columns for this core are s+1:e excluding any other headers
        measCols = (s+1):e;
        measCols = measCols(~isnan(nSR(1,measCols)));

        if isempty(measCols)
            continue; % core has no measurement bins
        end

        seg = nSR(:, measCols);

        % Merge within this core segment
        [segMerged, segLog] = merge_segment(seg, dt_min);

        % Append merged segment
        nSR_out = [nSR_out, segMerged];

        % Add logs with core number
        for k = 1:numel(segLog)
            logk = logk + 1;
            mergeLog(logk).core      = c;
            mergeLog(logk).j_in_core = segLog(k).j_in_core;
            mergeLog(logk).direction = segLog(k).direction;
            mergeLog(logk).dt_before = segLog(k).dt_before;
            mergeLog(logk).dt_after  = segLog(k).dt_after;
        end
    end
end


function [seg, segLog] = merge_segment(seg, dt_min)
% seg is 4xK, all columns are measurements (row1 not NaN)
% Merging rule: for any dt<dt_min, merge with left or right neighbor to
% make new dt as small as possible while >= dt_min (if possible in one step).
% If neither side reaches dt_min in one step, merge with the smaller dt-sum and repeat.

    segLog = struct('j_in_core', {}, 'direction', {}, 'dt_before', {}, 'dt_after', {});
    logk = 0;

    %put dts into their own vector
    dt = seg(3,:);

    %Get normalizer to convert SR estimates to nSR
    normalizer = (seg(2,1)./seg(3,1))./seg(1,1);

    %Check if any merging needs to happen
    while any(dt < dt_min)
        if size(seg,2) == 1
            warning('Only one bin left in a core and dt < dt_min; cannot merge further.');
            break;
        end

        %Find first bin to merge
        j = find(dt < dt_min, 1, 'first'); % policy: first offending bin

        %Check if it has bins on left and right
        hasL = (j > 1);
        hasR = (j < size(seg,2));

        %Calculate dt of merging with each possible bin
        dtL = inf; dtR = inf;
        if hasL, dtL = dt(j) + dt(j-1); end
        if hasR, dtR = dt(j) + dt(j+1); end

        % See if the merge with each bin brings it above dt_min
        canL = hasL && (dtL >= dt_min);
        canR = hasR && (dtR >= dt_min);

        % Decide which bin to merge with
        if canL && canR
            dir = ternary(dtL <= dtR, 'L', 'R');
        elseif canL
            dir = 'L';
        elseif canR
            dir = 'R';
        else
            dir = ternary(dtL <= dtR, 'L', 'R');
        end

        %Save the original dt
        dt_before = dt;

        % Do the merging
        if dir == 'L'
            % Merge j into j-1: add rows 2:3, recompute row1, delete column j
            seg(2:3, j-1) = seg(2:3, j-1) + seg(2:3, j);
            seg(1,   j-1) = (seg(2, j-1) / seg(3, j-1))/normalizer;
            seg(:, j) = [];
        else
            % Merge j into j+1: add rows 2:3, recompute row1, delete column j
            seg(2:3, j+1) = seg(2:3, j+1) + seg(2:3, j);
            seg(1,   j+1) = (seg(2, j+1) / seg(3, j+1))/normalizer;
            seg(:, j) = [];
        end

        % Get the new dt
        dt = seg(3,:);

        %Store all changes in the log
        logk = logk + 1;
        segLog(logk).j_in_core = j;
        segLog(logk).direction = dir;
        segLog(logk).dt_before = dt_before;
        segLog(logk).dt_after  = dt;
    end
end


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
