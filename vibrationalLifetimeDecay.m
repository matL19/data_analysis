function [VL_traces, t2s] = vibrationalLifetimeDecay(dataIn, peak_regions, peak_directions)

% This function takes in an SGR 2D-IR data structure of isotropic data,
% and calculates the vibrational lifetime of certain peaks.
% 
% ----- Input arguments: 
%
% dataIn - SGR 2D-IR data structure, assumed to already be isotropic and
% sorted by t2
%
% peak_regions - a cell array defining the bounds of each peak to be
% considered, where each element of the cell array is a 2x2 double with the
% w1 and w3 regions of the peak in question. ex:
%       peak_regions = {[
%           [ w1 range ];
%           [ w3 range ];
%       ]};
%
% peak_directions - a double array the same length as peak_regions
% specifiying whether each region represents a positive or negative peak.
% Input 1 for positive and -1 for negative.
%
% ----- Output arguments:
%
% VL_traces - a double array where the columns represent the vibrational
% lifetime traces. Column n is the nth peak requested in peak_regions.
%
% t2s - a double array where the columns represent the t2 axis for the
% vibrational lifetime traces. Column n is the nth peak requested in
% peak_regions.

% argument validation
if numel(peak_regions) ~= numel(peak_directions)
    error("peak_regions and peak_directions must be the same length.")
end
if any(~(peak_directions == 1 | peak_directions == -1))
    error("elements of peak_directions can only be 1 or -1.")
end

n_peaks = numel(peak_regions);
n_t2s = numel([dataIn.t2]);
VL_traces = zeros(n_t2s,n_peaks);
for ii = 1:n_peaks
    
    w1 = peak_regions{ii}(1,:);
    w3 = peak_regions{ii}(2,:);
    VL = [];
    
    d = cropData(dataIn,w1,w3);
    for jj = 1:numel(d)
        if peak_directions(ii) == 1
            y_new = max(max(d(jj).R));
        elseif peak_directions(ii) == -1
            y_new = min(min(d(jj).R));
        end
        VL = [VL y_new];
    end
    VL_traces(:,ii) = VL;
end
t2s = [dataIn.t2];

end