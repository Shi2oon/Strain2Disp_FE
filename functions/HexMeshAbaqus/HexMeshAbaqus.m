function [Nodes, Elements,xEl,yEl,zEl,E11,E22,E33,E12,E13,E23,ScaleYN] ...
            = HexMeshAbaqus(dvcdata)
close all;
    dvcdata = sortrows(dvcdata,[3,1,2]);
% organize data and getting steps
X  = dvcdata(:,1);              x = unique(X);          Sx = abs(x(2)-x(1));
Y  = dvcdata(:,2);              y = unique(Y);          Sy = abs(y(2)-y(1));
Z  = dvcdata(:,3);              z = unique(Z);          
try; Sz = abs(z(2)-z(1)); catch; Sz = z(1); end
[~,dataum ] = reshapeStrainData(dvcdata);
Dis_3D(:,:,:,1) = dataum.X1;     	Dis_3D(:,:,:,2) = dataum.Y1;        
Dis_3D(:,:,:,3) = dataum.Z1;        Dis_3D(:,:,:,4) = dataum.Exx;       
Dis_3D(:,:,:,5) = dataum.Eyy;       Dis_3D(:,:,:,6) = dataum.Ezz;
Dis_3D(:,:,:,7) = dataum.Exy;       Dis_3D(:,:,:,8) = dataum.Exz;      
Dis_3D(:,:,:,9) = dataum.Eyz;
[a, b, c, ~ ] = size(Dis_3D);
BinaryMask_3D = ones(a, b, c);

% add Nans to mask
[a]=isnan(squeeze(Dis_3D(:,:,:,5)));
BinaryMask_3D(a)=0;

if 1-nansum(BinaryMask_3D(:))/(size(BinaryMask_3D,1)*size(BinaryMask_3D,2)) > 0.01
    ScaleYN = 'Y';
else
    ScaleYN = 'N';
end

% Make point clouds of voxelized data
[BinaryMask_PtCloud, XYZmskDis] = MakePtClouds(BinaryMask_3D,Dis_3D);

% Create Hexahedral mesh in ABAQUS format from voxelized data
[VertMskCloud, Faces_MaskPtCloud] = HexMesh3D_Abaqus8Node(BinaryMask_PtCloud);

% Display Hexahedral mesh
DisplayAbaqusMesh(VertMskCloud, Faces_MaskPtCloud); 

dcdataO = [XYZmskDis(:,1).*Sx  XYZmskDis(:,2).*Sy  XYZmskDis(:,3).*Sz ...
           XYZmskDis(:,7)      XYZmskDis(:,8)      XYZmskDis(:,9) ...
           XYZmskDis(:,10)     XYZmskDis(:,11)     XYZmskDis(:,12)];
VertMskC = [VertMskCloud(:,1).*Sx  VertMskCloud(:,2).*Sy  VertMskCloud(:,3).*Sz];
% saveas(gcf,[INPfolder '\Hex_Mesh.fig']);  
% dcdataO = XYZmskDis(:,4:end);
% VertMskC = XYZmskDis(:,4:6);
% Write data to file .. need to modfiy the writting to resemble abaqus
[~, Nodes, Elements, ~, ~] = CreateAbaqusInputData_unique(VertMskC, dcdataO);
x = unique(dcdataO(:,1));	y = unique(dcdataO(:,2));	z = unique(dcdataO(:,3));       
for iV=1:8
    for iO=1:size(Elements,1)
        [~,Indx] = min(abs(x-Nodes(Elements(iO,iV),2)));      xEl(iV,iO) = x(Indx);
        [~,Indx] = min(abs(y-Nodes(Elements(iO,iV),3)));      yEl(iV,iO) = y(Indx);
        [~,Indx] = min(abs(z-Nodes(Elements(iO,iV),4)));      zEl(iV,iO) = z(Indx);
        A = ismember(dcdataO(:,1:3),[xEl(iV,iO) yEl(iV,iO) zEl(iV,iO)],'row');
        if sum(A)==0
%             Nodes(Elements(iO,iV))=[];
%         	  E11(iV,iO) = NaN;     E22(iV,iO) = NaN;
%             E33(iV,iO) = NaN;     E12(iV,iO) = NaN;
%             E13(iV,iO) = NaN;     E23(iV,iO) = NaN;
        else
            E11(iV,iO) = dcdataO(A,4);     E22(iV,iO) = dcdataO(A,5);
            E33(iV,iO) = dcdataO(A,6);     E12(iV,iO) = dcdataO(A,7);
            E13(iV,iO) = dcdataO(A,8);     E23(iV,iO) = dcdataO(A,9);
        end
    end
end
end