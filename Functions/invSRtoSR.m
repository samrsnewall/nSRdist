function[SRvals_interp, SRvalsprob_norm] = invSRtoSR(invSRvals, invSRprobs)

% Concatenate the invSRvals and probabilities at just above 0
invSRvalsPOS = invSRvals(invSRvals>=0.000001);
invSRprobsPOS = invSRprobs(invSRvals>=0.000001);

%Test whether the values are on equal spacing
if length(unique(diff(invSRvals))) ~= 1
    %Convert them to equal spacing
    invSRvals_equalspc = min(invSRvalsPOS):0.001:max(invSRvalsPOS);
    invSRprobs_equalspc = interp1(invSRvalsPOS, invSRprobsPOS, invSRvals_equalspc);
    
    %Assign these to the variables used in the next part of the code.
    invSRvalsPOS = invSRvals_equalspc;
    invSRprobsPOS = invSRprobs_equalspc;
end

%Apply correction for differing bin sizes after inversion
%find bin width in invSR space
invSRbinwidth = invSRvalsPOS(2) - invSRvalsPOS(1);
%use bin width in invSR space to find bin width in SR space
binsizes = 1./(invSRvalsPOS - invSRbinwidth./2) - 1./(invSRvalsPOS + invSRbinwidth./2);
%divide prob by binsize in SR space
SRvalsprob = invSRprobsPOS./binsizes;
%Find the SR values that apply to the bins
newSRvals = (1./(invSRvalsPOS - invSRbinwidth./2) + 1./(invSRvalsPOS + invSRbinwidth./2))./2;

%Interpolate to equal x-spacing from minimum SR value to 50 (gives plenty
%of space to capture high SR events)
SRvals_interp = newSRvals(end):0.001:50; %create new spacing
SRvalsprob_interp = interp1(newSRvals, SRvalsprob, SRvals_interp); %interpolate to it
%normalise prob vector
AUC = trapz(SRvals_interp, SRvalsprob_interp);
SRvalsprob_norm = SRvalsprob_interp./AUC;
end