close all
clear all
clc

% addpath(genpath(strcat(fileparts(which(mfilename)), '/subtightplot')));
data = '112717_mos21';

delimiterIn = '\t';
headerlinesIn = 3;
A = importdata(data, delimiterIn, headerlinesIn);
num_steps_x = str2double(A.textdata{1});
num_steps_y = str2double(A.textdata{2});

matSize = num_steps_x * num_steps_y;

% Peak to integrate over (in pixels) and the title
peak = {80 250};
beta0 = [7000000, 51, 115];
% Create a matrix with the spectra for each point as rows

Spectrum = zeros(num_steps_x*num_steps_y, 1340);

for i = 0:num_steps_y-1
    for j=1:num_steps_x
        k = j+num_steps_x*i;
        l = (i*((num_steps_x+1)*2-1)+j*2)-1;
        Spectrum(k,:)=A.data(l,:);
    end
end

%{
This is the threshold for a single pixel to increase by before being
subtracted (set to the previous pixel value). This removes cosmic rays to a
degree; occasionally they will be two pixels wide.
%}
cosmicThreshold = 300;
rowsSpectrum = size(Spectrum, 1);
colsSpectrum = size(Spectrum, 2);
for j=1:colsSpectrum
    minVal = 999999;
    for i=1:rowsSpectrum
        if i ~= 1 && i ~= rowsSpectrum && j ~= 1 && j ~= colsSpectrum
            lastVal = Spectrum(i - 1, j - 1);
            thisVal = Spectrum(i, j);
            nextVal = Spectrum(i + 1, j + 1);
            if thisVal - lastVal > cosmicThreshold && thisVal - nextVal > cosmicThreshold
                Spectrum(i, j) = lastVal;
            end
        end
    end
end   %end %loop over colsSpectrum
    
    %  fit spectra to Lorentzian or Gaussian and export fit params to Signal
    length = peak{1, 2}-peak{1, 1};
    
    pixel = zeros(1,length);
    for i = peak{1, 1}:peak{1, 2}
        pixel(1,i-peak{1, 1}+1)=i;
    end
    
    model = @(b,pixel)(b(1)./((b(2)-pixel).^2+b(3)^2));
    
    Signal = zeros(matSize,3);
    fprintf('starting loop\n')
    figure()
    hold on;
     for i = 0:matSize-1
        scratch = Spectrum(i+1,[peak{1, 1}:peak{1, 2}]);
        beta = nlinfit(pixel,scratch,model,beta0);
        %      fit = (beta(1)+beta(2)*exp(-((wn1-655).^2/(2*121)^2)));  % Gaussian
        fit = (beta(1)./((beta(2)-pixel).^2+beta(3)^2));  %Lorentzian Fit
        Signal(i+1,1) = beta(1);
        Signal(i+1,2) = beta(2);
        Signal(i+1,3) = beta(3)^2;
        plot(pixel,fit)
        plot(Spectrum(i+1,:))
        clear scratch
    end % loop over spectrum data

    %Put spectral data on grid
    
    Int = zeros(num_steps_y, num_steps_x, 1);
    coord = zeros(num_steps_y, num_steps_x, 1);
    Shift = zeros(num_steps_y, num_steps_x, 1);
    Width = zeros(num_steps_y, num_steps_x, 1);
    

for y = 0:num_steps_y-1
    for x = 1:num_steps_x
        index = (num_steps_x * y) + x;
        
        if mod(y, 2) == 0
            x_pos = x;
        else
            x_pos = num_steps_x - x + 1;
        end
        Int(y+1, x_pos, 1) = sum(Spectrum(index, peak{1, 1}:peak{1, 2}));
        coord (y+1,x_pos,1) = index;
        Shift (y+1,x_pos,1) = Signal(index,2);
        Width (y+1,x_pos,1) = Signal(index,3);
    end
end


figure
imagesc(Int(:, :, 1))
pbaspect([1 num_steps_y/num_steps_x 1])
colorbar
title('Intensity')
int_name = strcat(data,'_intensity');
savefig(int_name)
figure
imagesc(Shift(:, :, 1))
pbaspect([1 num_steps_y/num_steps_x 1])
colorbar
title('Peak Position')
Shift_name = strcat(data,'_Shift');
savefig(Shift_name)
figure
imagesc(Width(:, :, 1))
pbaspect([1 num_steps_y/num_steps_x 1])
colorbar
title('Peak Width')
width_name = strcat(data,'_Width');
savefig(width_name)
figure
imagesc(coord(:, :, 1))
pbaspect([1 num_steps_y/num_steps_x 1])
colorbar
title('Mapping')
map_name = strcat(data,'_mapping');
savefig(int_map)