function [Abaqus_DispBCs] = Abaqus8NodeDispBC_unique(DispPtCloud, UniqueABNodesMapped)

%| Inputs:
%|      UniqueABNodesMapped
%|              - element vertices derived from voxelized mesh
%|              - output by HexMesh3D_Abaqus8Node.m

%|      DispPtCloud
%|              - n x 6 matrix of point cloud data
%|                  Column 1-3: X, Y, Z vector position coordinates
%                           i.e., columns 1:3 MATCH Point Cloud coordinates
%|                  Column 4-6: X, Y, Z displacement magnitude

%| Outputs:
%|      Abaqus_DispBCs
%|              - n x 4 matrix of 3D Displacments for ABAQUS input file, as follows
%                       <node number>, 1, 1, <X displacement of node>
%                       <node number>, 2, 2, <Y displacement of node>
%                       <node number>, 3, 3, <Z displacement of node>

%**********************************************************************%
% disp('Creating Displacement BC');
dX=DispPtCloud(:,4);
dY=DispPtCloud(:,5);
dZ=DispPtCloud(:,6);

% PtCloudPts = number of coordinate points in segmented image
PtCloudPts = size(DispPtCloud,1);

% Sort ...Nodes... matrix to match order of coordinates in XYZdisp matrix
UniqueABNodesMapped_SORTED = sortrows(sortrows(sortrows(UniqueABNodesMapped, 4), 2), 3);

% Initialize Displacement Boundary condition
% Need to write the displacements of each voxel to the corresponding node
    % Tests if coordinates match, indicating node matches voxel
    % if so --> sets displacements from DispPtCloud matrix to the
    % corresponding node
DispBC_NodeNum = zeros(PtCloudPts,1);
Ptcloudrow = 1;
for n=1:length(UniqueABNodesMapped_SORTED)
    if Ptcloudrow <= PtCloudPts
        if UniqueABNodesMapped_SORTED(n,2:3) == DispPtCloud(Ptcloudrow,1:2)
            DispBC_NodeNum(Ptcloudrow,1) = UniqueABNodesMapped_SORTED(n,1);
            Noderow(Ptcloudrow,1) = n;
            Ptcloudrow = Ptcloudrow + 1;
        end
    end
end

% Double Check to make sure coordinates are equal
%{
for i=1:PtCloudPts
    nodematch(i) = isequal(UniqueABNodesMapped_SORTED(Noderow(i,1),2:3),DispPtCloud(i,1:2));
end
% i.e., nodesdontmatch = 0 if all coordinates in Nodes v. Displacement matrices dont match
nodesdontmatch = ismember(0, nodematch);    

% Display error message if coordinates in nodes matrix dont match those in point cloud 
if nodesdontmatch
    disp ('Nodes do not match');
end
%}

% initialize Displacement Boundary condition matrix
Abaqus_DispBCs = zeros(3*PtCloudPts, 4);

% ABAQUS BC format:     <node number>, 1, 1, <X displacement of node>
%                       <node number>, 2, 2, <Y displacement of node>
%                       <node number>, 3, 3, <Z displacement of node>

% Write Displacement Matrix
for i=1:PtCloudPts
    for j=1:3
        Abaqus_DispBCs(3*(i-1)+j, 1) = DispBC_NodeNum(i);
        Abaqus_DispBCs(3*(i-1)+j, 2) = j;
        Abaqus_DispBCs(3*(i-1)+j, 3) = j;
        if j == 1
            Abaqus_DispBCs(3*(i-1)+j, 4) = dX(i);
        else if j == 2
                Abaqus_DispBCs(3*(i-1)+j, 4) = dY(i);
            else if j == 3
                    Abaqus_DispBCs(3*(i-1)+j, 4) = dZ(i);
                end
            end
        end
    end
end
end