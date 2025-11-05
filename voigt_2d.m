function out = voigt_2d(x, y, x0, y0, w_g, w_l)
L = (w_l / pi)^2 ./ (w_l^2 + (x - x0).^2) ./ (w_l^2 + (y - y0).^2);

n_x = length(x);
dx = x(1, 2) - x(1, 1);
if mod(n_x, 2) == 0
    xx = (-n_x / 2: n_x / 2 - 1) * dx;
else
    xx = (-floor(n_x / 2): floor(n_x / 2)) * dx;
end
n_y = length(y);
dy = y(2) - y(1);
if mod(n_y, 2) == 0
    yy = (-n_y / 2: n_y / 2 - 1) * dy;
else
    yy = (-floor(n_y / 2): floor(n_y / 2)) * dy;
end
[XX, YY] = meshgrid(xx, yy);

D = (XX + YY) ./ sqrt(2);
A = (XX - YY) ./ sqrt(2);
da = sqrt(dx^2 + dy^2);
s_a = da / 2.355;
G = 1 / (2*pi*s_a*w_g).*exp(-D.^2 ./ (2 * w_g^2)).*exp(-A.^2 ./ (2 * s_a^2));
G = ifftshift(G);

F = real(ifft2(fft2(G).*fft2(L)));
norm = trapz(trapz(F));
out = F ./ norm;
