clc;clear;close all
% create a 30 MPa mode-I field (for mode II input 'II' for mixed mode input 'fun'
% for 3D field use the 'Westergaard_3D' function with similar input. then
% decide number of elements LxL. here I choose L = 12, to get faster
% results
[files,Operation,unit,stif] = Westergaard_Modes(30,'plane_stress','I',12);
% [files,Operation,unit,stif] = Westergaard_3D(90,'fun',5);
% load the strain field data with [um] units you can choose other fields 
% but be aware that figures will always assumes the input units are in um
% which will not alter in any thing when it comes to caclulated field
load(files{2,3})
percentTry = 100:-10:10;    ctx = 1; 
    while ctx<length(percentTry)+1
        try
            % reduce data size if needed
		[inData,Dir.ActP] = smotherData(alldata,percentTry(ctx)); 
[X1,Y1,Ux,Uy,M4,Z1,Uz] = FE_OOM(alldata,'Linear',erase(files{2,1},'.mat'));
        catch err
            ctx = ctx+1;
            if ctx == length(percentTry)+1;     rethrow(err);       end
        end
    end
	

% to go to abaqus you need to declear some mechanical values
MatP.input_unit = unit{2,3}; % units to adjust in abaqus which is unitless
MatP.stressstat = Operation{2,3}; % plane stress or strain
MatP.Mat = 'ferrite'; % materil name (any will do)
MatP.type = 'E'; % E for istropic linear elastic, A for anistropic with 
                 % stinfess tensor as input and R for elastic plastic mateirals
MatP.E = 220e9;  % young's modulus
MatP.nu = 0.3;   % possiov's ratio
MatP.unique = 'calibration'; % unique name
[BCf, UnitOffset] = ...
    AbaqusBC(M4.Nodes,M4.Elements, [X1(:) Y1(:) Ux(:) Uy(:)], MatP, [], erase(files{2,1},'.mat'));

% you can then upload the .inp file to abaqus
% if you faced with any issues or you thinking about extending this code to
% other things please drop me a line at abdo.aog@gmail.com or leave a
% comment in GitHub