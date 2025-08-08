function[mu, var] = muVarPDFVec(pdf)
mu = trapz(pdf.x, pdf.x.*pdf.px);
var = trapz(pdf.x, pdf.x.*pdf.x.*pdf.px) - mu^2;
end