function [X,Y,Ux,Uy] = gatherGrains(scan_dir)
set(0,'defaultAxesFontSize',25); 	set(0,'DefaultLineMarkerSize',14) 
d = dir(fullfile(scan_dir)); % add _Pyxpb %scan all the files with .nxs format
count=0;
for iv = 1:length(d)
    if d(iv).isdir == 1 && ~strcmpi(d(iv).name, '.') && ~strcmpi(d(iv).name, '..')
        count = count+1;
        Dut{count} = load(fullfile(d(iv).folder,d(iv).name,'Intergated Uxy.mat'));
        if count ==1   
            Ux = zeros(size(Dut{count}.Ux));    
            Ux(Ux==0)=NaN;                      Uy = Ux; 
            X  = Dut{count}.X1;                 Y  = Dut{count}.Y1;
        elseif size(Dut{count}.Ux) ~= size(Dut{1}.Ux)
            Dut{count}.Ux = flip(griddata(Dut{count}.X1,Dut{count}.Y1,Dut{count}.Ux,X,Y),2);
            Dut{count}.Uy = flip(griddata(Dut{count}.X1,Dut{count}.Y1,Dut{count}.Uy,X,Y),2);
        end
        
%%
        for  iX=1:size(Dut{count}.Ux,1)
            for iY=1:size(Dut{count}.Ux,2)
                if ~isnan(Dut{count}.Ux(iX,iY))
                    Ux(iX,iY) = Dut{count}.Ux(iX,iY) ;
                    Uy(iX,iY) = Dut{count}.Uy(iX,iY) ;
                end
            end
s1=subplot(1,2,1);	imagesc(unique(X),unique(Y),Ux);        set(gca,'Ydir','normal'); 
s1.XDir='reverse';	s1.YDir='reverse'; axis image;axis xy;  colormap jet;%axis off;
xlabel('X[\mum]');	ylabel('Y[\mum]');   
c = colorbar;     	c.Label.String = 'Ux[\mum]';       v1 = caxis;

s1=subplot(1,2,2);	imagesc(unique(X),unique(Y),Uy);        set(gca,'Ydir','normal'); 
s1.XDir='reverse';	s1.YDir='reverse'; axis image;axis xy;  colormap jet;%axis off;
xlabel('X[\mum]');	ylabel('Y[\mum]');   
c = colorbar;    	c.Label.String = 'Uy[\mum]';       v2 = caxis;
set(gcf,'position',[10 50 1900 950]); 
pause(0.1)
        end
    end
end

%% plot
% figure
close all
s1=subplot(1,2,1);	imagesc(unique(X),unique(Y),Ux);        set(gca,'Ydir','normal'); 
s1.XDir='reverse';	s1.YDir='reverse'; axis image;axis xy;  colormap jet;%axis off;
xlabel('X[\mum]');	ylabel('Y[\mum]');   
c = colorbar;     	c.Label.String = 'Ux[\mum]';       v1 = caxis;

s1=subplot(1,2,2);	imagesc(unique(X),unique(Y),Uy);        set(gca,'Ydir','normal'); 
s1.XDir='reverse';	s1.YDir='reverse'; axis image;axis xy;  colormap jet;%axis off;
xlabel('X[\mum]');	ylabel('Y[\mum]');   
c = colorbar;    	c.Label.String = 'Uy[\mum]';       v2 = caxis;
set(gcf,'position',[10 50 1900 950]); 
saveas(gcf,[scan_dir '\All_Dis.tif'],'tiffn'); 
saveas(gcf,[scan_dir '\All_Dis.tif']); 
save([scan_dir '\All Intergated Uxy.mat'],'X','Y','Ux','Uy'); close