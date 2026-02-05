function out = lorentzian(x, x0, wL)
out = (wL / pi) ./ (wL^2 + (x - x0).^2);
