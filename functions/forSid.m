clear; clc; close all
range =  5:2.5:100;    range=round(range,0); % size range*range
for iV = 1:length(range)
    resultsDir = ['P:\Abdo\EBSD Data\30 WS\30MPa_[mm]_Strain_'  num2str(range(iV))  'S_1G'];
    if ~exist([resultsDir '\Linear\Integrated Uxy.mat'],'file')
        load(['P:\Abdo\EBSD Data\30 WS\30MPa_[mm]_Strain_' num2str(range(iV)) '.mat'], 'alldata');
        tic;        xEBSD2DIS(alldata,1,[resultsDir '\Linear'],'Linear');
        Time = toc; save([resultsDir '\Linear\Integrated Uxy.mat'],'Time','-append'); 
        tic;        xEBSD2DIS(alldata,4,[resultsDir '\Quadratic'],'Quadratic');
        Time = toc; save([resultsDir '\Quadratic\Integrated Uxy.mat'],'Time','-append');
        clear alldata
    end
	load(['P:\Abdo\EBSD Data\30 WS\30MPa_[mm]_DISP_' num2str(range(iV)) '.mat'], 'alldata');
	[X,Y,Ux,Uy,~,~,~] = FE_Mesh_Generator([alldata alldata(:,3)],2,1,'Linear');
    load([resultsDir '\Linear\Integrated Uxy.mat'],'M4','Time');
    x(1,iV) = mean(abs(X-M4.X1));       y(1,iV)  =  mean(abs(Y-M4.Y1));
    UX(1,iV) =  mean(abs(Ux-M4.Ux));    UY(1,iV) =  mean(abs(Uy-M4.Uy));
    TT (1,iV) = Time;
    clear  M4 Time
    
    [X,Y,Ux,Uy,~,~,~] = FE_Mesh_Generator([alldata alldata(:,3)],2,1,'Linear');
    load([resultsDir '\Quadratic\Integrated Uxy.mat'],'M4','Time');
    x(2,iV) = mean(abs(X-M4.X1));       y(2,iV)  =  mean(abs(Y-M4.Y1));
    UX(2,iV) =  mean(abs(Ux-M4.Ux));    UY(2,iV) =  mean(abs(Uy-M4.Uy));
    TT (2,iV) = Time;
    clear alldata M4 Time
end
plot3(range.^2,UY(1,:),TT(1,:));    hold on
plot3(range.^2,UY(2,:),TT(12,:));   hold off
xlabel('El. No.');          ylabel('Error (abs.)');     zlabel('Time (s)');
set(gca,'ZScale', 'log');   zticks([10, 100,1000,10000]); 
legend('Linear','Quadratic');
set(gcf,'position',[10 50 1300 950]); 
saveas(gcf, [fileparts(resultsDir) '\Q vs L Error.fig']);  saveas(gcf, [InoF '\Mean Error.tif']); close