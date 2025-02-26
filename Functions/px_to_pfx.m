function[fx, pfx] = px_to_pfx(x, px, f)

% Concatenate x and px above x = 0
xPOS = x(x>=1e-5);
pxPOS = px(x>=1e-5);

%Test whether the x is on equal spacing
dxU = unique(diff(x));                                                      %Find unique values of dx
dxU_residual = abs(dxU - dxU(1));                                           %Find how different any unique values are from each other
dxU_big = dxU_residual(1:end)>1e-5; dxU_big(1) = true;                                                         
dxCounts = sum(dxU_big);                                                %Find out how many extra dx values are actually significant

if dxCounts ~= 1
    %Convert the x to equal spacing
    xbinwidth = 5e-4;
    x_equalspc = min(xPOS):xbinwidth:max(xPOS);
    px_equalspc = interp1(xPOS, pxPOS, x_equalspc);
    
    %Assign these to the variables used in the next part of the code.
    xPOS = x_equalspc;
    pxPOS = px_equalspc;
else
    xbinwidth = dxU(1);
end

%%%% Apply correction for differing bin sizes after applying function

%find bin width in fx space
fx_halfstepdown = f(xPOS - xbinwidth./2);
fx_halfstepup = f(xPOS + xbinwidth./2);
fx_binsizes = abs(fx_halfstepdown - fx_halfstepup);
%divide prob by binsize in fx space
fxvalsprob = pxPOS./fx_binsizes;
%Find the fx values that apply to the bins
fxvals = (fx_halfstepdown + fx_halfstepup)./2;

%Interpolate to x spacing equivalent to original x spacing but after
%function applied (make sure no part of it is outside of fxvals, to avoid
%NaNs)
fxvals_interp = sort(f(x)); %create new spacing
fx_toobig = fxvals_interp > max(fxvals);
fx_toosmall = fxvals_interp < min(fxvals);
fxvalsprob_interp = interp1(fxvals, fxvalsprob, fxvals_interp); %interpolate to it
fxvalsprob_interp(fx_toobig) = fxvalsprob_interp(find(fx_toobig, 1, 'first') -1);
fxvalsprob_interp(fx_toosmall) = fxvalsprob_interp(find(fx_toosmall, 1, 'last') +1);
%normalise prob vector
AUC = trapz(fxvals_interp, fxvalsprob_interp);
fxvalsprob_norm = fxvalsprob_interp./AUC;

%Assign to values to be output
fx = fxvals_interp;
pfx = fxvalsprob_norm;
end