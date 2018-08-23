refSpec = importdata('/Users/mollyamay/Documents/MATLAB/2018/March/032718/0327180_1.txt');
x_in_cm = 1/0.00005328-1./(refSpec(:,1)*10^-7);

% Region of interest in cm^-1 (WHAT PEAKS ARE YOU LOOKING AT?)
rois = 360;
roie = 410;

% Peak Height, Peak position, gaussian fwhm, lorentzian fwhm

guess = [4000 375 30 30 ... 
    4000 400 30 30 ...
    0 0 3000];    
%4000 441 50 50 ...
% Guess range to fit. (guess plus/minus this value)
guess_delta = [4000 5 30 30 ...
     4000 3 30 30 ...
     100000 100 80000];
     %4000 3 40 40 ...

% to fix peak parameters, just delete the corresponding number (the
% subsequent numbers are unchanged
free_parameters = [1 2 3 4 ...
    5 6 7 8 ... 
    9 10 11];
   % 9 10 11 12 ...
   
        y_intensity = transpose(Spectrum(4000,:));
        
        % Region of interest
        roi_start = find(x_in_cm >= rois, 1);
        roi_end = find(x_in_cm >= roie, 1);
        roi = [roi_start:roi_end];
        roguess = find(x_in_cm >= 0.5, 1);
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
        
        [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, y_intensity(roi), x_in_cm(roi), 1);
        [f, G, fit, out] = fitvoigt(answer, y_intensity(roi), x_in_cm(roi), 1);
   
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
        
 
      

%         % Region of interest
%         roi_start = find(x_in_cm >= rois, 1);
%         roi_end = find(x_in_cm >= roie, 1);
%         roi = [roi_start:roi_end];
%         roguess = find(x_in_cm >= 50, 1);
%             hf=figure(1);
%             clf;
%     
% 
%         % Set baseline param
%         guess(length(guess) - 2) = y_intensity(roguess); 
%         guess_delta(length(guess) - 2) = y_intensity(roguess)/2;
%         
%         guess(length(guess)) = max(y_intensity);
%         guess_delta(length(guess)) = max(y_intensity)/2;
%         
%         %Setup guesses
%         high_guess = guess + guess_delta;
%         low_guess = guess - guess_delta;
% 
%         % Run the fitting procedure
%         
%         [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, y_intensity(roi), x_in_cm(roi), 1);
%         [f, G, fit, out] = fitvoigt(answer, y_intensity(roi), x_in_cm(roi), 1);
%        
%             figure
%             plot(out{1}, out{2}-out{5}, out{1}, out{3}-out{5});
%             title('178,11')
%             ylabel('Intensity (arb. u.)')
%             xlabel('Raman Shift (cm^-^1)')
%       
% %        results(:,1) = answer; 
% 
% % figure
% % plot(x_in_cm,refSpec(:,3))
% % xlim([rois roie])
% % axis 'auto y'
