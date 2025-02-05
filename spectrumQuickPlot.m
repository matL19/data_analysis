function [fig, ax] = spectrumQuickPlot(freq, data, varargin)
fig = figure();
ax = axes(fig);
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "xlim"
            x_limits = val;
        case "ylim"
            y_limits = val;
        otherwise
            error("Invalid name/value pair")
    end
    varargin = varargin(3:end);
end
plot(ax,freq,data)
xlabel('Frequency (cm^{-1})')
ylabel('Absorbance (AU)')
if exist("x_limits", 'var')
    xlim(x_limits)
end
if exist("y_limits", 'var')
    ylim(y_limits)
end
end