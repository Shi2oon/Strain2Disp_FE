function Maps = inandOut(fname,Spec)
    load([fname '_XEBSD.mat'],'C_voight', 'GND','Map_EBSD_MTEX','MicroscopeData',...
        'Grain_Map_A0_sample','Data', 'Grain_Map_rotation_sample',...
        'Grain_Map_strain_sample', 'Grain_Map_stress_sample','GrainData','Maps',...
        'Grain_Map_PH2','Grain_Map_MAE2');
SavingD = [fname '_results\full_map'];
%% save data
A = squeeze(Grain_Map_A0_sample(:,:,Spec,:,:));             % DEFORMAION
W = squeeze(Grain_Map_rotation_sample(:,:,Spec,:,:));       % ROTATION
E = squeeze(Grain_Map_strain_sample(:,:,Spec,:,:));         % STRAIN
S = squeeze(Grain_Map_stress_sample(:,:,Spec,:,:));         % STRESS
% save GNDS PH and MAE
[rowC,colZ] = find(S(:,:,1,1));  % linear indices for nonzero element
Maps.GND    = zeros(size(squeeze(S(:,:,1,1))));
try for i=1:length(colZ);   Maps.GND(rowC(i),colZ(i)) = GND.total(rowC(i),colZ(i)); end; end
try Maps.PH  = mean(squeeze(Grain_Map_PH2(:,:,Spec,:)),3);
Maps.MAE = mean(squeeze(Grain_Map_MAE2(:,:,Spec,:)),3); end
% Rotation
Maps.W11 = W(:,:,1,1);      Maps.W12 = W(:,:,1,2);      Maps.W13 = W(:,:,1,3);  
Maps.W21 = W(:,:,2,1);      Maps.W22 = W(:,:,2,2);      Maps.W23 = W(:,:,2,3);  
Maps.W31 = W(:,:,3,1);      Maps.W32 = W(:,:,3,2);      Maps.W33 = W(:,:,3,3);  
% Stress
Maps.S11 = S(:,:,1,1);      Maps.S12 = S(:,:,1,2);      Maps.S13 = S(:,:,1,3);  
Maps.S21 = S(:,:,2,1);      Maps.S22 = S(:,:,2,2);      Maps.S23 = S(:,:,2,3);  
Maps.S31 = S(:,:,3,1);      Maps.S32 = S(:,:,3,2);      Maps.S33 = S(:,:,3,3); 
% Strain
Maps.E11 = E(:,:,1,1);      Maps.E12 = E(:,:,1,2);      Maps.E13 = E(:,:,1,3);  
Maps.E21 = E(:,:,2,1);      Maps.E22 = E(:,:,2,2);      Maps.E23 = E(:,:,2,3);  
Maps.E31 = E(:,:,3,1);      Maps.E32 = E(:,:,3,2);      Maps.E33 = E(:,:,3,3); 
% Deformation gradient
Maps.A11 = A(:,:,1,1);      Maps.A12 = A(:,:,1,2);      Maps.A13 = A(:,:,1,3);  
Maps.A21 = A(:,:,2,1);      Maps.A22 = A(:,:,2,2);      Maps.A23 = A(:,:,2,3);  
Maps.A31 = A(:,:,3,1);      Maps.A32 = A(:,:,3,2);      Maps.A33 = A(:,:,3,3);

% stifness:  crystal orientation is defined as the rotation that transforms crystal 
% coordinates, i.e., a description of a vector or a tensor with respect to the crystal
% reference frame, into specimen coordinates, i.e., a desciption of the same object
% with respect to a specimen fixed reference frame.
Maps.Stiffness = StiffnessRot(Map_EBSD_MTEX(sub2ind([MicroscopeData.NROWS,MicroscopeData.NCOLS],...
                              GrainData.RefPoint.prop.yi(Spec),GrainData.RefPoint.prop.xi(Spec)))...
                               .orientations.matrix,C_voight(:,:,Spec));            
Maps.nu  =  Maps.Stiffness(1,2)/(Maps.Stiffness(1,1)+ Maps.Stiffness(1,2)); %meaningless,not used
Maps.E   =  Maps.Stiffness(1,1)*(1-2*Maps.nu)*(1+Maps.nu)/(1-Maps.nu);%meaningless,not used
% Dim
Maps.X   = Data.XSample;    Maps.Y   = Data.YSample;  
Maps.stepsize  =(abs(Maps.X(1,1)-Maps.X(1,2)));         
Maps.Wo  = (1/2).*(Maps.S11.*Maps.E11 + Maps.S12.*Maps.E12 + Maps.S13.*Maps.E13 +...
                   Maps.S21.*Maps.E21 + Maps.S22.*Maps.E22 + Maps.S23.*Maps.E23 +...
                   Maps.S31.*Maps.E31 + Maps.S32.*Maps.E32 + Maps.S33.*Maps.E33);
% units (defualt xEBSD units)
Maps.units.xy = 'um';       Maps.units.S  = 'GPa';         
Maps.units.E  = 'Abs.';     Maps.units.W = 'rad'; 
alldata = [Maps.X(:)   Maps.Y(:)   Maps.E11(:)   Maps.E22(:) Maps.E12(:)];
save(SavingD,'Maps','alldata'); % save
[Maps.X1,Maps.Y1,Maps.Ux,Maps.Uy,Maps.M4,Maps.ScaleYN] = xEBSD2DIS(alldata,'Linear',fileparts(SavingD));
save(SavingD,'Maps','alldata'); % save
end