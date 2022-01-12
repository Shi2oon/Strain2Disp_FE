%% 2D
restoredefaultpath;clc;clear;close all
addpath(genpath([pwd '\functions'])); PWD = pwd;
% create a 30 MPa mode-I field (for mode II input 'II' for mixed mode input 'fun'
% for 3D field use the 'Westergaard_3D' function with similar input. then
% decide number of elements LxL. here I choose L = 12, to get faster
% results
[files,Operation,unit,E,nu] = Westergaard_Modes(30,'plane_stress','I');
% load the strain field data with [um] units you can choose other fields 
% but be aware that figures will always assumes the input units are in um
% which will not alter in any thing when it comes to caclulated field
load(files{2,3})
percentTry = 100:-10:10;    ctx = 1; 
    while ctx<length(percentTry)+1
        try
            % reduce data size if needed
		[inData,MatP.ActP] = smotherData(alldata,percentTry(ctx)); 
[MatP.X1,MatP.Y1,MatP.Ux,MatP.Uy,MatP.M4] = ...
    FE_OOM(inData,'Linear',erase(files{2,1},'.mat'));
    ctx=length(percentTry)+1;
        catch err
            ctx = ctx+1;
            if ctx == length(percentTry)+1;     rethrow(err);       end
        end
    end
	
% to go to abaqus you need to declear some mechanical values
% and you need to download DIC2ABAQUS from https://github.com/Shi2oon/DIC2ABAQUS
MatP.input_unit = unit{2,3}; % units to adjust in abaqus which is unitless
MatP.Operation = Operation{2,3}; %for integreated data
MatP.stressstat = 'plane_stress'; % plane stress or strain
MatP.Mat = 'ferrite'; % materil name (any will do)
MatP.type = 'E'; % E for istropic linear elastic, A for anistropic with 
                 % stinfess tensor as input and R for elastic plastic mateirals
MatP.E = E;  % young's modulus
MatP.nu = nu;   % possiov's ratio
MatP.unique = 'calibration'; % unique name
MatP.results = erase(files{2,1},'.mat');

clc;cd('A:\OneDrive - Nexus365\GitHub\DIC2ABAQUS');
addpath(genpath(pwd));  
% Locate the crack 
[DATA,UnitOffset,MatP,msk,SaveD] = Locate_Crack(MatP,MatP.input_unit,MatP.results,MatP); 
% prepare and run abaqus cae
[MatP.Abaqus,Abaqus.CAE] = PrintRunCode(MatP, ...
        msk,SaveD,ceil(min(size(DATA.X1))*0.5-2),UnitOffset);
% Post Processing
[Abaqus.J,Abaqus.Keff,Abaqus.KI,Abaqus.KII] = PlotKorJ(MatP.Abaqus,MatP.E,UnitOffset);

cd(PWD)
save([fileparts(MatP.Abaqus) '\AbaqusOutput.mat'],'Abaqus','MatP');

%% 3D
restoredefaultpath;clc;clear;close all
addpath(genpath([pwd '\functions'])); PWD = pwd;
[files,Operation,unit,stif] = Westergaard_3D([20,10,50],'fun',8);
load(files{2,3})
percentTry = 100:-10:10;    ctx = 1; 
    while ctx<length(percentTry)+1
        try
            % reduce data size if needed
		[inData,MatP.ActP] = smotherData(alldata,percentTry(ctx)); 
[MatP.X1,MatP.Y1,MatP.Ux,MatP.Uy,MatP.M4,MatP.Z1,MatP.Uz] = ...
    FE_OOM(inData,'Linear',erase(files{2,3},'.mat'));
    ctx=length(percentTry)+1;
        catch err
            ctx = ctx+1;
            if ctx == length(percentTry)+1;     rethrow(err);       end
        end
    end

% to go to abaqus you need to declear some mechanical values
% and you need to download DVC2Abaqus from https://github.com/Shi2oon/DVC2ABAQUS
cd('A:\OneDrive - Nexus365\GitHub\DVC2ABAQUS')
addpath(genpath(pwd));  
MatP.Operation = 'DVC';
MatP.input_unit   = unit{2,3};        % meter (m) or milmeter (mm) or micrometer(um);
MatP.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
MatP.results = [erase(files{2,3},'.mat') '\3D_Integrated_Uxy'];
MatP.stressstat = 'DVC'; % plane stress or strain
MatP.Mat = 'ferrite'; % materil name (any will do)
MatP.type = 'E'; % E for istropic linear elastic, A for anistropic with
% stinfess tensor as input and R for elastic plastic mateirals
MatP.E = 210e9;  % young's modulus
MatP.nu = 0.3;   % possiov's ratio
MatP.unique = 'calibration'; % unique name
[offset,RadEulerAng,rotCentre,Abaqus,Len] = DVC2J(MatP);% input to abaqus
% Plot results
[J,KI,KII,KIII,Direction] = Plot3DKJ(erase(Abaqus,'.inp'),offset,Len(3,1));