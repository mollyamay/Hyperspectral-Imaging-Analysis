%refSpec = importdata('/Users/mollyamay/Documents/MATLAB/2018/May/052818/0425181_1.txt');
%  refSpec = importdata('/Users/mollyamay/Documents/MATLAB/2017/October/20171003_HS_453.csv');
% x_in_cm = 1/0.00005328-1./(refSpec(:,1)*10^-7);

%x_in_cm = refSpec(:,1);
%x_in_cm = 4.13567*10.^-6.*299792458./refSpec(:,1);
figure
% plot(refSpec.data(:,1),refSpec.data(:,2))
% plot(x_in_cm, refSpec(:,3))
% figure

for i = 1:50:size(Spectrum);
    plot(Spectrum(i,:))  %plot(x_in_cm, Spectrum(i,:))
    hold on
end
