function [h, p] = errorbar_fill(x, y, err, options)
%ERRORBAR_FILL Confidence interval as filled band
%   Plot the data points and confidence interval
%   Right now use a weighted average for 2 or less points and use the
%   observed std of the points for 3 or more points.
arguments (Input)
    x (1,:)
    y (1,:)
    err (1,:)
    options.ax = gca;
    options.Marker = ".";
    options.MarkerSize = 1;
    options.Color = [1, 0, 0];
    options.FaceColor = [1, 0, 0];
    options.FaceAlpha = 0.1;
    options.EdgeColor = "none";
    options.DefaultErr = 0.05;
    options.MinErr = 0.025;
    options.XScale = "log";
    options.YScale = "linear";
end

arguments (Output)
    h
    p
end

ax = options.ax;

[xbar, ybar] = determine_fill_border(x, y, err, options);

p = draw_patch(ax, xbar, ybar, options);

hold on
h = plot(ax, x, y);
set(ax, ...
    "XScale", options.XScale, ...
    "YScale", options.YScale)
set(h, ...
    "Marker", options.Marker, ...
    "Color", options.Color, ...
    "LineStyle", "none")
hold off
end


function [xconf,yconf] = determine_fill_border(x, y, err, options)
[xu, ~, ic] = unique(x);
yp = zeros(size(xu));
yn = zeros(size(xu));
for ii = 1:length(xu)
    % y_mean = weighted_average(y(ic == ii), err(ic == ii));
    % delta = rms(err(ic == ii));
    ind = ic == ii;
    n = sum(ind);
    if n > 2
        % if several points just use them to est error
        [delta, y_mean] = std(y(ind));
        delta = delta * 2;
    else
        % for few points propagate error
        [y_mean, delta] = weighted_average(y(ind), err(ind), options);
    end
    yp(ii) =  y_mean + delta; % Calculate upper border
    yn(ii) =  y_mean - delta; % Calculate lower border
end

xconf = [xu, xu(end:-1:1)];
yconf = [yp, yn(end:-1:1)];
end


function p = draw_patch(ax, xconf, yconf, options)
hold on
p = fill(ax, xconf, yconf, 'red');
p.FaceColor = options.FaceColor;
p.FaceAlpha = options.FaceAlpha;
p.EdgeColor = options.EdgeColor;
end


function out = rms(in)
out = sqrt(mean(in.^2, "omitnan")); % Calculate the root mean square of the input
end


function [ym, e] = weighted_average(y, err, options)
err(isnan(err)) = options.DefaultErr;
err(err<options.MinErr) = options.MinErr;
inv_var = 1./(err.^2);
% inv_var(isnan(inv_var)) = mean(inv_var, "omitnan");
sum_inv_var = sum(inv_var, "omitnan");
ym = sum(y .* inv_var) / sum_inv_var; % Calculate the weighted average
e = 2/sqrt(sum_inv_var);  % rms *2 for 95%
end
