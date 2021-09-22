function [UniqueABNodes, UniqueABNodesMapped, UniqueABElements, UniqueAB_StepBCs, missingnodes] ...
    = CreateAbaqusInputData_unique(Vertices, XYZDispPtCloud)


%| Inputs:
%|      vertices_MaskPtCloud
%|              - element vertices derived from voxelized mesh
%|              - output by HexMesh3D_Abaqus8Node.m

%|      XYZDispPtCloud
%|              - displacement vector: n x 6 matrix
%|                  Column 1-3: X, Y, Z vector position coordinates
%|                  Column 4-6: X, Y, Z displacement magnitude
%|              - output by MakePtClouds.m

%| Outputs:
%|      UniqueABNodes
%|              - Unique Node coordinates for ABAQUS input file (*.inp)

%|      UniqueABElements
%|              - Unique Element definitions for ABAQUS input file (*.inp)

%|      UniqueAB_StepBCs
%|              - 3D Displacments for ABAQUS input file (*.inp)

%|      missingnodes
%|              - error handling
%|                  returns # missing nodes if all nodes are not included
%|                  in element definitions

[Vert_rows] = size(Vertices,1);

% Initialize Nodes Matrix & populate with vertices
AbaqusNodes = zeros(Vert_rows,4);
Indices = 1:1:Vert_rows;
Indices = Indices.';

AbaqusNodes(:,1) = Indices;
AbaqusNodes(:,2:4) = Vertices;

% Remove redundant nodes and keep first instance
[UniqueABNodes, ~, indexUnique] = unique(AbaqusNodes(:,2:4), 'rows', 'stable');
% First column = index of vertex corresponding to first instance of
% coordinates
ABNodesMapped(:,1) = indexUnique;
% Columns 2-4 = voxel coordinates, X,Y,Z
ABNodesMapped(:,2:4) = AbaqusNodes(:,2:4);

% Initialize Unique Nodes Matrix & populate with first instance of vertex coordinates
UniqueABNodesMapped=zeros(size(UniqueABNodes,1),4);
% First column = index of vertex corresponding to first instance of
% coordinates
UniqueABNodesMapped(:,1) = Indices(1:size(UniqueABNodes,1));
% Columns 2-4 = voxel coordinates, X,Y,Z
UniqueABNodesMapped(:,2:4) = UniqueABNodes(:,1:3);

% Initialize Unique Elements Matrix
UniqueABElements = zeros(size(XYZDispPtCloud,1), 9, 'int32');

% First column = Element number
% Columns 2-9 = node numbers corresponding to 8 vertices of element
    % NOTE order = IMPORTANT
for q=1:size(XYZDispPtCloud,1)
    UniqueABElements(q, 1) = q;
    UniqueABElements(q, 2) = ABNodesMapped(8*(q-1)+1,1);
    UniqueABElements(q, 3) = ABNodesMapped(8*(q-1)+2,1);
    UniqueABElements(q, 4) = ABNodesMapped(8*(q-1)+4,1);
    UniqueABElements(q, 5) = ABNodesMapped(8*(q-1)+3,1);
    UniqueABElements(q, 6) = ABNodesMapped(8*(q-1)+5,1);
    UniqueABElements(q, 7) = ABNodesMapped(8*(q-1)+6,1);
    UniqueABElements(q, 8) = ABNodesMapped(8*(q-1)+8,1);
    UniqueABElements(q, 9) = ABNodesMapped(8*(q-1)+7,1);
end

% disp('1. Missing Node Error Handling Disabled - Edit CreateAbaqusInputData_unique.m line 76 to re-enable.');
missingnodes = 'Disabled';
% Test to ensure no missing nodes in Element Matrix
%disp('Testing to ensure no missing nodes in Element Matrix');
%tic
%containnode = ones(size(UniqueABNodes,1),1);
%for n = 1:size(UniqueABNodes,1)
%    containnode(n) = ismember(n, UniqueABElements(:,2:9));
%    
%end
%missingnodes = ismember(0, containnode);    % i.e., missingnodes = 0 if all nodes included in Element matrix


% Create Matrix of Displacement magnitudes in X,Y,Z for voxels within segmented image for input to Abaqus FEM

% disp('2. Create Matrix of Displacement magnitudes in X,Y,Z for voxels within Point Cloud');
UniqueAB_StepBCs = Abaqus8NodeDispBC_unique(XYZDispPtCloud, UniqueABNodesMapped);
end