function [vert_ObjPtCloud, Faces_ObjPtCloud] = HexMesh3D_Abaqus8Node(PtCloudCoords)

% Function creates 8*vertices per voxel, which creates redundant nodes,
    % these are removed in subsequent functions (e.g., those with suffix
    % _unique.m
% PtCloudCoords is a 3-column matrix of point cloud data, 
    % e.g., representing voxelized segmentation of an image
        % X = column 1
        % Y = column 2
        % Z = column 3
        
disp('Generates vertices and faces for a 3D Hexahedral Mesh')
X=PtCloudCoords(:,1);
Y=PtCloudCoords(:,2);
Z=PtCloudCoords(:,3);


% PtCloudPts = number of coordinate points in segmented image
[PtCloudPts] = size(PtCloudCoords,1);

% initialize vertices and faces matrices
vert_ObjPtCloud = zeros(8*PtCloudPts, 3);
Faces_ObjPtCloud = zeros(6*PtCloudPts, 4);

% *** THE ORDER OF +/- IS IMPORTANT HERE in order for the 'patch' function to
    % properly connect vertices AND ALSO for ABAQUS to properly define
    % elements/faces
for q=1:PtCloudPts
    vert_ObjPtCloud(8*(q-1)+1, 1) = X(q);
    vert_ObjPtCloud(8*(q-1)+2, 1) = X(q)+1;
    vert_ObjPtCloud(8*(q-1)+3, 1) = X(q);
    vert_ObjPtCloud(8*(q-1)+4, 1) = X(q)+1;
    vert_ObjPtCloud(8*(q-1)+5, 1) = X(q);
    vert_ObjPtCloud(8*(q-1)+6, 1) = X(q)+1;
    vert_ObjPtCloud(8*(q-1)+7, 1) = X(q);
    vert_ObjPtCloud(8*(q-1)+8, 1) = X(q)+1;
    
    vert_ObjPtCloud(8*(q-1)+1, 2) = Y(q);
    vert_ObjPtCloud(8*(q-1)+2, 2) = Y(q);
    vert_ObjPtCloud(8*(q-1)+3, 2) = Y(q)+1;
    vert_ObjPtCloud(8*(q-1)+4, 2) = Y(q)+1;
    vert_ObjPtCloud(8*(q-1)+5, 2) = Y(q);
    vert_ObjPtCloud(8*(q-1)+6, 2) = Y(q);
    vert_ObjPtCloud(8*(q-1)+7, 2) = Y(q)+1;
    vert_ObjPtCloud(8*(q-1)+8, 2) = Y(q)+1;
    
    vert_ObjPtCloud(8*(q-1)+1, 3) = Z(q);
    vert_ObjPtCloud(8*(q-1)+2, 3) = Z(q);
    vert_ObjPtCloud(8*(q-1)+3, 3) = Z(q);
    vert_ObjPtCloud(8*(q-1)+4, 3) = Z(q);
    vert_ObjPtCloud(8*(q-1)+5, 3) = Z(q)+1;
    vert_ObjPtCloud(8*(q-1)+6, 3) = Z(q)+1;
    vert_ObjPtCloud(8*(q-1)+7, 3) = Z(q)+1;
    vert_ObjPtCloud(8*(q-1)+8, 3) = Z(q)+1;

% Define mesh element faces for each cubic element

    Faces_ObjPtCloud(6*(q-1)+1,:) = 8*(q-1)+[1 2 4 3];
    Faces_ObjPtCloud(6*(q-1)+2,:) = 8*(q-1)+[5 6 8 7];
    Faces_ObjPtCloud(6*(q-1)+3,:) = 8*(q-1)+[1 5 6 2];
    Faces_ObjPtCloud(6*(q-1)+4,:) = 8*(q-1)+[3 4 8 7]; 
    Faces_ObjPtCloud(6*(q-1)+5,:) = 8*(q-1)+[2 4 8 6];
    Faces_ObjPtCloud(6*(q-1)+6,:) = 8*(q-1)+[1 3 7 5];

end
end
