function [ax] = plot_fit_data_errorbar_fill(...
    x_fit, y_fit, x_data, y_data, err, options...
    )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    x_fit (1,:)
    y_fit (1,:)
    x_data (1,:)
    y_data (1,:)
    err (1,:)
    options.ax = gca;
    options.LineWidth = 1.5;
    options.Marker = ".";
    options.MarkerSize = 4;
    options.MarkerFaceColor = "none";
    options.Color = [1, 0, 0];
    options.FaceColor = [1, 0, 0];
    options.FaceAlpha = 0.2;
    options.EdgeColor = "none";
    options.DefaultErr = 0.05;
    options.MinErr = 0.025;
    options.XScale = "log";
    options.YScale = "linear";
end

arguments (Output)
    ax
end
ax = options.ax;

hold on
plot(ax,x_fit, y_fit, '-', 'LineWidth', options.LineWidth, 'Color', options.Color);

options = rmfield(options,"LineWidth");
nv_pairs = namedargs2cell(options);
errorbar_fill(x_data, y_data, err, nv_pairs{:});

end