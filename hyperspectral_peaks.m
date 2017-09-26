mydata = 'dentin-enamel-40s-A';
delimiterIn = '\t'
headerlinesIn = 3;
A = importdata(mydata, delimiterIn, headerlinesIn)
num_stepsx = str2num(A.textdata{1})
num_stepsy = str2num(A.textdata{2})


% Peak range to start. Run the program first and look at the all plot data
% to get the proper range for the peak you want to fit.
peakStart = 895;
peakFinish = 930;

% Wether or not to create new figures on each run (otherwise output to 666
% and 667)
saveFig = false;

%{
This is the threshold for a single pixel to increase by before being 
subtracted (set to the previous pixel value). This removes cosmic rays to a 
degree; occasionally they will be two pixels wide.
%}
cosmicThreshold = 500;

matSize = num_stepsx*num_stepsy;
Signal = zeros(num_stepsx*num_stepsy,4);

% Create the x column 

for i = 0:num_stepsy-1
    for j=1:num_stepsx
        k = j+num_stepsx*i;
        l = i*((num_stepsx+1)*2-1)+j*2;
        Signal(k,3)=A.data(l,1);
    end
end

% Create the y column 

for i = 0:num_stepsy-1   
    for j=1:num_stepsx
        k = (num_stepsx)*((i+1)*2)+(i+1);
        Signal(j+num_stepsx*i,2)=A.data(k,1);
    end
end

% % Create the total intensity column 
% 
% for i = 0:num_steps-1
%     
%      for j=1:num_steps
%         k = j+num_steps*i;
%         l = i*((num_steps+1)*2-1)+j*2;
%         Signal(k,1)=A.data(l,2);
% end
% end

% Create a matrix with the spectra for each point as rows 

for i = 0:num_stepsy-1
    for j=1:num_stepsx
        k = j+num_stepsx*i;
        l = (i*((num_stepsx+1)*2-1)+j*2)-1;
        Spectrum(k,:)=A.data(l,:);
    end
end


%{
Remove cosmic ray values. This works by checking if adjacent pixels to the
current one change by the `cosmicThreshold` value in a peak like manner.  

e.g. pixels p1, p2, p3 have values of 0, 10000, and 50 respectivly. 
%}
rowsSpectrum = size(Spectrum, 1);
colsSpectrum = size(Spectrum, 2);
for j=1:colsSpectrum
    minVal = 999999;
    for i=1:rowsSpectrum
        if i ~= 1 && i ~= rowsSpectrum && j ~= 1 && j ~= colsSpectrum
            if thisVal < minVal
               minVal = thisVal 
            end
            lastVal = Spectrum(i - 1, j - 1);
            thisVal = Spectrum(i, j);
            nextVal = Spectrum(i + 1, j + 1);  
            if thisVal - lastVal > cosmicThreshold && thisVal - nextVal > cosmicThreshold
                Spectrum(i, j) = lastVal;
            end
        end
    end
    %Spectrum(i, :) -= (zeros(1, colsSpectrum) + minVal);
end

% Sum each spectra to create intensities between start and finish pixels

avg = floor((peakStart+peakFinish) ./ 2);
avg

for i = 1:matSize
   Peak_1(i,1) = sum(Spectrum(i, peakStart:peakFinish));% - Spectrum(i, avg);% * (peakFinish - peakStart);
end

    
% Fit each spectrum to a gaussian 
%{
wn1 = zeros(1,colsSpectrum);
for i = 1:colsSpectrum
  wn1(1,i)=i;
end
 
 model = @(b,wn1)(b(1)+b(2)*exp(-((wn1-655).^2/(2*121)^2)));
 
 rng('default');
 
 for i = 0:matSize
     beta0 = [100,1000];
     scratch = Spectrum(i+1,:);
     beta = nlinfit(wn1,scratch,model,beta0);
     fit = (beta(1)+beta(2)*exp(-((wn1-655).^2/(2*121)^2)));
     Signal(i+1,4) = beta(1);
     plot(wn1,fit)
     hold on;
     plot(x1,spectrum2)
     clear scratch
 end

%}
% Make a contour plot of the intensities in z

z=Peak_1(:,1);  %plots peak
y=Signal(:,2);
x=Signal(:,3);


%z = z./max(z);  % normalize intensities

dx=0.001;
dy=0.001;

x_edge = [floor(min(x)*10)/10:dx:ceil(max(x*10)/10)];
y_edge = [floor(min(y)):dy:ceil(max(y))];
[X,Y] = meshgrid(x_edge,y_edge);

Z = griddata(x, y, z, X, Y);

if saveFig
    figure(666)
else
    figure
end


surf(X, Y, Z, 'EdgeColor', 'interp')
ylabel('y position (\mum)')
xlabel('x position (\mum)')
zlabel('Intensity (arb. units)')
title(mydata)

% Plot all spectra
if saveFig
    figure(667)
else
    figure
end

% TODO plot specific ranges
zz = transpose(Spectrum);
numRows = size(zz, 1);
numCols = size(zz, 2);
xx = zeros(numRows, numCols);
yy = zeros(numRows, numCols);
for i=1:numRows
    for j=1:numCols
        xx(i, j) = j;
        yy(i, j) = i;
    end
end
plot3(xx, yy, zz, 'Color', [0.2 0.2 0.8])

ylabel('pixel')
xlabel('sample #')
zlabel('Intensity (arb. units)')