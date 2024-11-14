function[Y, pY] = inverseProbDist(X, pX)
% Given:
% X  - a vector of x values (unevenly spaced)
% pX - the pdf values corresponding to X

% Step 1: Transform X to Y
Y = 1 ./ X;

% Step 2: Compute the new PDF pY
pY = pX ./ (Y.^2);  % Use the change of variables formula

% Optional: Sort Y and pY because 1/X inverts the order of values
[Y, sortIdx] = sort(Y);  % Sort Y values
pY = pY(sortIdx);        % Rearrange pY to match sorted Y

%Compare to answer from invSRtoSR
[invX, pinvX] = invSRtoSR(X, pX);

% Plot original and transformed distributions for comparison
figure;
subplot(2,1,1);
plot(X, pX);
title('Original Gamma Distribution pX vs X');
xlabel('X');
ylabel('pX');

subplot(2,1,2);
hold on
plot(Y, pY, "k-");
plot(invX, pinvX, "r--");
title('Transformed Distribution pY vs Y (Y = 1/X)');
xlabel('Y');
ylabel('pY');
end