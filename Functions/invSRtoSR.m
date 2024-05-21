function[SRvals_interp, SRvalsprob_norm] = invSRtoSR(invSRvals, invSRprobs)

% Concatenate the invSRvals and probabilities at just above 0
invSRvalsPOS = invSRvals(invSRvals>=0.000001);
invSRprobsPOS = invSRprobs(invSRvals>=0.000001);

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
SRvalsprob_norm = SRvalsprob_interp./sum(SRvalsprob_interp);
end