function out = calibrateILConc(data,freq,ref1,freq1,ref2,freq2,varargin)
%     Returns a 2x1 matrix M = [A;B] of the relative concentrations of IL and PEGDA
%     in your sample where A = fraction of reference 1 and B = fraction of
%     reference 2.
%     Function call of the form:
%     M = calibrateILConc(data,freq,ref1,freq1,ref2,freq2)
% Optional inputs include
% calibrateILConc(...,"manual",X)
% This calibrates the input spectrum according to the values provided in X,
% instead of using the GUI click-input method.
% where X = [ref_freq1 baseline_freq1 ref_freq2 baseline_freq2]
%
%     arguments
%         data double
%         freq (:,1) double
%         ref1 double
%         freq1 (:,1) double
%         ref2 double
%         freq2 (:,1) double
%         varargin
%     end

%     if freq ~= freq1 & freq ~= freq2
%         error('All three frequency axes must be the same.')
%         return
%     end

while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "manual"
            inputmethod = "manual";
            calibration_pts = val;
        otherwise
            error("Invalid name/value pair")
    end
    varargin = varargin(3:end);
end
if ~exist('inputmethod')
    inputmethod = "click";
end
if inputmethod == "click"
    str = {'Select your peaks IN THIS ORDER:',...
        '1. Top of reference peak 1',...
        '2. Baseline of reference peak 1',...
        '3. Top of reference peak 2',...
        '4. Baseline of reference peak 2'};
    
    figure(8412);clf
    set(figure(8412),'Units','normalized')
    set(figure(8412),'Position',[0.3333 0.0162 0.3315 0.8848])
    subplot(2,1,1)
    plot(freq1,ref1)
    annotation('textbox',[0.55 0.4 0.5 0.5],'String',str,'FitBoxToText','on')
    title('Your first reference spectrum')
    subplot(2,1,2)
    plot(freq2,ref2)
    title('Your second reference spectrum')
    
    pause
    
    [X,Y] = ginput(4);
elseif inputmethod == "manual"
    X = calibration_pts;
    dw = abs(freq1(2)-freq1(1));
    Y(1) = ref1(freq1 >= X(1)-dw/2 & freq1 <= X(1)+dw/2);
    Y(2) = ref1(freq1 >= X(2)-dw/2 & freq1 <= X(2)+dw/2);
    dw = abs(freq2(2)-freq2(1));
    Y(3) = ref2(freq2 >= X(3)-dw/2 & freq2 <= X(3)+dw/2);
    Y(4) = ref2(freq2 >= X(4)-dw/2 & freq2 <= X(4)+dw/2);
end
% find the coefficients from the pure spectra
epsilon = [0 0;0 0];
%range = emimntf2_pure(freq1 > 1550 & freq1 < 1600);
epsilon(1,1) = Y(1)-Y(2);
dw = abs(freq1(2)-freq1(1));
idx1 = find(freq1 >= X(1)-dw/2 & freq1 <= X(1)+dw/2);
if numel(idx1) > 1
    idx1 = idx1(end/2);
end
idx2 = find(freq1 >= X(2)-dw/2 & freq1 <= X(2)+dw/2);
epsilon(1,2) = ref2(idx1) - ref2(idx2);
%range = pegda_pure(freq2 > 1700 & freq2 < 1750);
epsilon(2,2) = Y(3)-Y(4);
dw = abs(freq2(2)-freq2(1));
idx1 = find(freq2 >= X(3)-dw/2 & freq2 <= X(3)+dw/2);
if numel(idx1) > 1
    idx1 = idx1(end/2);
end
idx2 = find(freq2 >= X(4)-dw/2 & freq2 <= X(4)+dw/2);
epsilon(2,1) = ref1(idx1) - ref1(idx2);


% solve for the concentrations in the input spectrum
absorbance = [0;0];
%range = data(freq > 1550 & freq < 1600);
dw = abs(freq(2)-freq(1));
idx1 = find(freq >= X(1)-dw/2 & freq <= X(1)+dw/2);
idx2 = find(freq >= X(2)-dw/2 & freq <= X(2)+dw/2);
absorbance(1) = data(idx1) - data(idx2);

%range = data(freq > 1700 & freq < 1750);
idx1 = find(freq >= X(3)-dw/2 & freq <= X(3)+dw/2);
idx2 = find(freq >= X(4)-dw/2 & freq <= X(4)+dw/2);
absorbance(2) = data(idx1) - data(idx2);

concs = epsilon\absorbance;
fractions = 1/sum(concs)*concs;

fprintf(1,'Concentration determination / relative ratios \n')
for ii = 1:length(fractions)
    fprintf(1,'    component %2i',ii);
end
fprintf(1,'\n')
for ii = 1:length(fractions)
    fprintf(1,'%16.6f',fractions(ii));
end
fprintf(1,'\n')

out = fractions;

end