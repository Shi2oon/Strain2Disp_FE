function [OBJMapPtCloud, XYZDispPtCloud] = MakePtClouds(ObjMap_crop,dispXYZ)

% object_bin_val is the value of the binary that corresponds to the
% segmentation of interest, e.g.

% Reduce Object Map to include only segmentation of interest (i.e.,
% corresponding to binary value)
ObjMap_crop(ObjMap_crop ~= 1) = 0;
ObjMap_crop(ObjMap_crop == 1) = 1;

% Determine dimensions of segmented volume
[objmap_rows, objmap_cols, objmap_slices] = size(ObjMap_crop);

% Create (Binary) Point Cloud of segmented volume
% disp('Segmenting of Displacement field')
a = 1;
for i=1:objmap_rows
    for j = 1:objmap_cols
        for k = 1:objmap_slices
            if ObjMap_crop(i,j,k)==1
                OBJMapPtCloud(a,1)=j;
                OBJMapPtCloud(a,2)=i;
                OBJMapPtCloud(a,3)=k;
                a = a+1;
            end
        end
    end
end

% Get Displacement (Transformation) Data from ANTs Registration
X_disp  = dispXYZ(:,:,:,1);
Y_disp  = dispXYZ(:,:,:,2);
Z_disp  = dispXYZ(:,:,:,3);
XX_disp = dispXYZ(:,:,:,4);
YY_disp = dispXYZ(:,:,:,5);
ZZ_disp = dispXYZ(:,:,:,6);
XY_disp = dispXYZ(:,:,:,7);
XZ_disp = dispXYZ(:,:,:,8);
YZ_disp = dispXYZ(:,:,:,9);

% Create Point Cloud of Displacements from Transformation Matrix
% 6 columns, 1:3 = X,Y,Z coordinates of points in segmented volume = MATCH
% ObjectMap Point Cloud
% columns 4:6 are displacements in the X,Y,Z directions

XYZDispPtCloud = zeros(length(OBJMapPtCloud),9);
XYZDispPtCloud(:,1:3) = OBJMapPtCloud;
b = 1;
for m=1:objmap_rows
    for n = 1:objmap_cols
        for p = 1:objmap_slices
            if ObjMap_crop(m,n,p)==1
                XYZDispPtCloud(b,4)  = X_disp(m,n,p);
                XYZDispPtCloud(b,5)  = Y_disp(m,n,p);
                XYZDispPtCloud(b,6)  = Z_disp(m,n,p);
                XYZDispPtCloud(b,7)  = XX_disp(m,n,p);
                XYZDispPtCloud(b,8)  = YY_disp(m,n,p);
                XYZDispPtCloud(b,9)  = ZZ_disp(m,n,p);
                XYZDispPtCloud(b,10) = XY_disp(m,n,p);
                XYZDispPtCloud(b,11) = XZ_disp(m,n,p);
                XYZDispPtCloud(b,12) = YZ_disp(m,n,p);
                b = b+1;
            end
        end
    end
end

end