varname= genvarname('0711175')
data2=load('0711175');
[m,n]=size(data2)

for i=1:1:m
    Laser(i,1)=0;
    for j=1:1:4
        Laser(i,1)=Laser(i,1)+data2(i,25+j);
    end
end

Laser_norm=Laser./max(Laser);

figure
% plot(frequency,Laser_norm)
plot(Laser_norm)
title('Laser')

for i=1:1:m
    Total(i,1)=0;
    for j=1:1:n
        Total(i,1)=Total(i,1)+data2(i,j);
    end
end

% Total_norm=Total9./max(Total9);

figure
% plot(frequency,Total_norm)
plot(Total)
title(varname)

% Laser_Total1=Total8_norm./Laser1_norm;
% ZPL_Total1=Total8_norm./ZPL1_norm;

% figure(4)
% plot(frequency,ZPL_Total1)
% plot(ZPL_Total1)
% title('ZPL/Total1')
% figure(5)
% plot(frequency,Laser_Total1)
% plot(Laser_Total1)
% title('Laser/Total1')

FileName=[varname]
save(FileName)

%% For laser line, use starting data point 26, da=1:1:5
%% For ZPL, use starting data point 422 and da=1:1:23 
%% For total intensity, use starting data point 1 and da=1:1:1340