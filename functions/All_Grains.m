function [Maps]=All_Grains(fname,email_add,started,ended)
tic; set(0,'defaultAxesFontSize',35); 	set(0,'DefaultLineMarkerSize',10) 

%% load data
load(fname,'Map_RefID','C_rotated', 'Data_InputMap',...
    'Grain_Map_strain_sample', 'GrainData','Data');

%% plot to select the boundray
try
    close all;  	warning off;        s3=subplot(1,1,1);
    imagesc(Data_InputMap.X_axis,Data_InputMap.Y_axis,Map_RefID); hold on
    axis off;           axis image;         axis xy;          
    colormap jet;       s3.XDir='reverse';  s3.YDir='reverse';                 
if isempty(GrainData.RefPoint)==1
    [GrainData.RefPoint] = to_label(Data_InputMap, Map_RefID);
else
    for iv=1:length(GrainData.RefPoint.x)       
        GrainData.RefPoint.prop.labels{iv} = num2str(iv);     
    end
end
scatter(GrainData.RefPoint.x,GrainData.RefPoint.y,'k','filled');
scatter(GrainData.RefPoint.x,GrainData.RefPoint.y,'w');
labelpoints(GrainData.RefPoint.x,GrainData.RefPoint.y,GrainData.RefPoint.prop.labels);     
hold off;      set(gcf,'position',[30 50 1300 950])
[A,B,~] = fileparts(fname);
SavingD = fullfile(A,[B ' Integrated Grains']);       mkdir(SavingD);
saveas(gcf,[SavingD '\All Grains.png']);     
saveas(gcf,[SavingD '\All Grains.fig']);     close

%% save data
if ~exist('started','var')    
        GrinT = 1:length(GrainData.RefPoint.x);
elseif length(started)==1 && exist('ended','var')
        GrinT = started:1:ended;
else 
        GrinT = started;
end
       
for Spec = GrinT
    clc; fprintf('Started Solving for grain no. %d\n',Spec);
    E = squeeze(Grain_Map_strain_sample(:,:,Spec,:,:));         % STRAIN
    Maps.E11 = E(:,:,1,1);      Maps.E12 = E(:,:,1,2);      Maps.E22 = E(:,:,2,2);      
    
    if Spec == GrinT(1)      % Dim.
        Maps.Stiffness = C_rotated;                     
        Maps.Y   = Data.YSample;                Maps.X   = Data.XSample;
        Maps.stepsize  = (abs(Maps.X(1,1)-Maps.X(1,2))); 
        Maps.nu  =  Maps.Stiffness(1,2)/(Maps.Stiffness(1,1)+ Maps.Stiffness(1,2));
        Maps.E   =  Maps.Stiffness(1,1)*(1-2*Maps.nu)*(1+Maps.nu)/(1-Maps.nu);
        Maps.units.xy = 'um';                   Maps.units.S  = 'GPa';         
        Maps.units.E  = 'Abs.';                 Maps.units.W = 'rad';  
%         Maps.X = Maps.X';  Maps.Y = Maps.Y';
    end
    
    DoS = MoreOutNaN(Maps.E12,unique(Maps.X));
    Maps.E11(isnan(DoS))=NaN;      Maps.E12(isnan(DoS))=NaN;      Maps.E22(isnan(DoS))=NaN; 
%%
    percentTry = [100,85,60,45,30,15,7.5];    ctx = 1; % to reduce data density when data is huge
    ActP = 0; % data size is not reduced
        while ctx<length(percentTry)+1
            try    
                [Map,ActP] = smotherData(Maps,percentTry(ctx));
                savTo = fullfile(SavingD,[num2str(Spec) '_' num2str(ActP)]); 	
                        mkdir(savTo);
                inData = [Map.X(:),Map.Y(:),Map.E11(:),Map.E22(:),Map.E12(:)];
                [~,~,~,~,M4{Spec}] = xEBSD2DIS(inData,2,savTo);    
                ctx = length(percentTry)+1;
            catch err
                w = warning;    warning on;  warning(err.message);  warning(w.state);
                rmdir(savTo);
                ctx = ctx+1;
                if ctx == length(percentTry)+1            
                    rethrow(err);            
                end   
            end
        end
	M4Elements{Spec} = M4{Spec}.Elements;
end
if started ==1 &&  ended ==  length(GrainData.RefPoint.x)
    Inp4Grains(M4{started}.Nodes,M4Elements,SavingD);   
    gatherGrains(SavingD);
end

%%
done = {['All_Grains is done integrateing for ' num2str(GrinT)...
        ' Grain(s) in ' num2str(toc/3600) ' hours and saved in ' SavingD]};
gmail(email_add,['xEBSD2ABAQUS on ' getenv('COMPUTERNAME') ' for ' getenv('USERNAME')],done)
catch err
    emailErr(err,SavingD,email_add);                 
    rethrow(err)
end