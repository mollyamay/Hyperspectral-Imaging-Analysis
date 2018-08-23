% df = '032718/0327182';
 refSpec = importdata('/Users/mollyamay/Documents/MATLAB/2018/June/062518_hs/062018_9.csv');
% %x_in_nm = 4.13566*10^-6*299792458./(refSpec(:,1));
x_in_nm = refSpec(:,1);
% figure
% plot(x_in_nm,refSpec(:,3))
% Region of interest in cm^-1 (WHAT PEAKS ARE YOU LOOKING AT?)
rois = 640;
roie = 900;

x_index = 1;
y_index = 3;

guess = [10000 660 50 50 ...
    10000 810 80 80 ...
    0 -2.5e-3 3000];    

guess_delta = [10000 20 50 50 ...
    10000 50 80 80 ...
     10000 10 10000];
     
free_parameters = [1 2 3 4 ...
    5 6 7 8 ...
    9 10 11];
  
y_intensity = transpose(Spectrum(41, :));  %refSpec(:,3);

        % Region of interest
        roi_start = find(x_in_nm >= rois, 1);
        roi_end = find(x_in_nm >= roie, 1);
        roi = [roi_start:roi_end];
        roguess = find(x_in_nm >= 0.5, 1);
            hf=figure(1);
            clf;
    

        % Set baseline param
        guess(length(guess) - 2) = y_intensity(roguess); 
        guess_delta(length(guess) - 2) = y_intensity(roguess)/2;
        
        guess(length(guess)) = max(y_intensity);
        guess_delta(length(guess)) = max(y_intensity)/2;
        
        %Setup guesses
        high_guess = guess + guess_delta;
        low_guess = guess - guess_delta;

        % Run the fitting procedure
        
        [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, y_intensity(roi), x_in_nm(roi), 1);
        [f, G, fit, out] = fitvoigt(answer, y_intensity(roi), x_in_nm(roi), 1);
   
            % -------- Plot Fit -------- %
            subplot(2,1,1)
            plot(out{1}, out{2}-out{5}, out{1}, out{3}-out{5});
            title('')
            ylabel('Intensity (arb. u.)')
            xlabel('')
  

        % -------- Process data -------- %
   
       results = transpose(answer);

        % No need to show baseline fit, hence -3 parameters
        table_rows = (length(answer) - 3) / 4;
        table_data = zeros(table_rows, 5);
        
        for j=1:table_rows
            index = (j - 1) * 4 + 1;
            gauss_fwhm = answer(index + 3);
            lorentz_fwhm = answer(index + 2);
            table_data(j, 1) = answer(index); % Intensity
            table_data(j, 2) = answer(index + 1); % Peak Pos
            table_data(j, 3) = lorentz_fwhm;
            table_data(j, 4) = gauss_fwhm;
            %voigt peak width
            table_data(j, 5) = gauss_fwhm*(1-2.0056*1.0593+sqrt((lorentz_fwhm/gauss_fwhm)^2+2*1.0593*lorentz_fwhm/gauss_fwhm+2.0056^2*1.0593^2));

        end

       
            %-------- Plot data -------- %
%             MATLAB trickery, produce a subplot, get its position and delete it.
%             Then put the uitable into the subplot position

            sp = subplot(2, 1, 2);
            pos = get(sp, 'Position');
            un = get(sp, 'Units');
            delete(sp);
            cnames={'Amp', 'Position', 'Lorentzian FWHM', 'Gaussian FWHM', 'Voigt FWHM'};
            t = uitable(hf, 'Data', table_data, 'ColumnName', cnames, 'Units', un, 'Position', pos);
        
 
      


