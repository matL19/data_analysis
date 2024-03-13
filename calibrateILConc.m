function out = calibrateILConc(data,freq,ref1,freq1,ref2,freq2)
    % Returns a 2x1 matrix M = [I;P] of the relative concentrations of IL and PEGDA
    % in your sample where I = IL fraction and P = PEGDA fraction. 

    arguments
        data double
        freq (:,1) double
        ref1 double
        freq1 (:,1) double
        ref2 double
        freq2 (:,1) double
    end
    
    if freq ~= freq1 & freq ~= freq2
        error('All three frequency axes must be the same.')
        return
    end

    str = {'Select your peaks IN THIS ORDER:',...
        '1. Top of reference peak 1',...
        '2. Baseline of reference peak 1',...
        '3. Top of reference peak 2',...
        '4. Baseline of reference peak 2'};
    
    figure(8412);clf
    subplot(2,1,1)
    plot(freq1,ref1)
    annotation('textbox',[0.55 0.4 0.5 0.5],'String',str,'FitBoxToText','on')
    title('Your first reference spectrum')
    subplot(2,1,2)
    plot(freq2,ref2)
    title('Your second reference spectrum')
    
    pause
    
    [X,Y] = ginput(4);
    
    % find the coefficients from the pure spectra
    epsilon = [0 0;0 0];
    %range = emimntf2_pure(freq1 > 1550 & freq1 < 1600);
    epsilon(1,1) = Y(1)-Y(2);
    dw = abs(freq1(2)-freq1(1));
    idx1 = find(freq1 >= X(1)-dw/2 & freq1 <= X(1)+dw/2);
    epsilon(1,2) = ref2(idx1);
    %range = pegda_pure(freq2 > 1700 & freq2 < 1750);
    epsilon(2,2) = Y(3)-Y(4);
    dw = abs(freq2(2)-freq2(1));
    idx2 = find(freq2 >= X(3)-dw/2 & freq2 <= X(3)+dw/2);
    epsilon(2,1) = ref1(idx2);
    
    % solve for the concentrations in the input spectrum
    absorbance = [0;0];
    %range = data(freq > 1550 & freq < 1600);
    idx1 = find(freq >= X(1)-dw/2 & freq <= X(1)+dw/2);
    absorbance(1) = data(idx1);
    %range = data(freq > 1700 & freq < 1750);
    idx2 = find(freq >= X(3)-dw/2 & freq <= X(3)+dw/2);
    absorbance(2) = data(idx2);
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