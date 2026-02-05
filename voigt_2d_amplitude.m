function out = voigt_2d_amplitude(x_in, y_in, x0, y0, w_g, w_l)

n = 2^8;
method = 'cubic';

ll = min(x_in(1), y_in(1));
ul = max(x_in(end), y_in(end));
n_x = n;
n_y = n;
x = linspace(ll, ul, n);
dx = (ul - ll) / (n - 1);
dy = dx;
[X, Y] = meshgrid(x, x);

% lorentzian
L = 1 ./ (1 + ((X - x0) / w_l).^2) ./ (1 + ((Y - y0) ./ w_l).^2);

% setup axes in fft order (zero is first element)
xx = [(0:floor(n_x / 2) - 1), (-floor(n_x / 2):-1)] * dx;
yy = [(0:floor(n_y / 2) - 1), (-floor(n_y / 2):-1)] * dy;
[XX, YY] = meshgrid(xx, yy);

% diagonal gaussian starting at 0
D = (XX + YY) ./ sqrt(2);
A = (XX - YY) ./ sqrt(2);
da = sqrt(dx^2 + dy^2);
s_a = da / 2.355;  % about 1 point wide ~ delta fxn
G = exp(-D.^2 ./ (2 * w_g^2)) .* exp(-A.^2 ./ (2 * s_a^2));

% convolution
G_hat = fft2(G);
G_hat = G_hat ./ G_hat(1);
L_hat = fft2(L);
F = real(ifft2(G_hat .* L_hat));

amplitude = interp2(X, Y, F, x0, y0, method);
F = F ./ amplitude;
out = interp2(X, Y, F, x_in, y_in, method);
