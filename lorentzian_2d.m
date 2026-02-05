function out = lorentzian_2d(x, y, x0, y0, w_l, varargin)

out = (w_l / pi)^2 ./ (w_l^2 + (x - x0).^2) ./ (w_l^2 + (y - y0).^2);
