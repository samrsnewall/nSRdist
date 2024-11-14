% Original matrix A
d2coast = readmatrix("dist2coast5pt.csv", "NumHeaderLines", 2);
A = d2coast;

%%% Convert A into a format where we can plot it as a matrix
% Get the unique x and y values
x_vals = unique(A(:,1));  % Unique x values
y_vals = unique(A(:,2));  % Unique y values

% Initialize the matrix B of size [length(x_vals), length(y_vals)]
B = zeros(length(x_vals), length(y_vals));

% Use sub2ind to map x and y values to linear indices of matrix B
% Get the row indices corresponding to x values
[~, x_idx] = ismember(A(:,1), x_vals);

% Get the column indices corresponding to y values
[~, y_idx] = ismember(A(:,2), y_vals);

% Use sub2ind to compute the linear indices for matrix B
linear_idx = sub2ind(size(B), x_idx, y_idx);

% Assign z values to their respective positions in B
B(linear_idx) = A(:,3);

%%% Remove information from the land
Bocean = B;
Bocean(B <=0) = 0;

%%%

%Create a reasonable depth colormap map
depthMap = [0 0 0; linspace(1,0,10000)',linspace(1,0,10000)',ones(10000,1)];

%Find lat and long limits
x = [min(x_vals) max(x_vals)];
y = [min(y_vals) max(y_vals)];

%Plot as a map
figure; imagesc(y, x, Bocean); colormap(depthMap); colorbar; set(gca, 'YDir', 'normal')
hold on
plot(longs(highSRCoresLog), lats(highSRCoresLog), 'rs', "LineWidth", 1)
plot(longs(lowSRCoresLog), lats(lowSRCoresLog), 'bs', "LineWidth", 1)