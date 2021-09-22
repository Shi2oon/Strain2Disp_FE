function [nodes,elem,X,Y,Z,E1,E2,E3,E4,E5,E6] = Mesh2and3D(rawData,DIM)
if DIM == 2
    xy     = rawData(:,1:2);
    Uxy    = rawData(:,3:5);
    [~,order] = sortrows(xy,[1,2]);
    Uxy = Uxy(order,:);
     
    % get the number of vox along x, y, z
    xVox = size(unique(xy(:,1)),1);
    yVox = size(unique(xy(:,2)),1);

    % get step size and boundaries for the mesh
    step   = zeros(1,3);
    boxmin = zeros(1,3);
    boxmax = zeros(1,3);
    for ii = 1:2
        step(ii) = min(diff(unique(xy(:,ii))));
        boxmin(ii) = xy(1,ii)-step(ii)/2;
        boxmax(ii) = xy(end,ii)+step(ii)/2;
    end

    % generate nodes 
    [x,y] = meshgrid(boxmin(1):step(1):boxmax(1),boxmin(2):step(2):boxmax(2));
    numNodes = numel(x);
    coord = [reshape(x,numNodes,1), reshape(y,numNodes,1)];
    nodes = [(1:numNodes)', sortrows(coord,[2,1])];

    % allocate array for elements
    elem = zeros(size(xy,1),5);
    count = 1;

    % start loop over voxel dimensions
    for ix = 1:xVox
            for iy = 1:yVox
                % get element label
                elem(count,1) = count;

                % nodes on the plane with higher x
                elem(count,2) = iy + (0-1)*(yVox+1) + ix*(yVox+1)*(0+1);
                elem(count,3) = elem(count,2) + 1;
                elem(count,4) = elem(count,2) + xVox + 2;
                elem(count,5) = elem(count,2) + xVox + 1;
                count = count+1;
            end
    end
    
    F1 = scatteredInterpolant(xy(:,1),xy(:,2),Uxy(:,1),'natural');
    F2 = scatteredInterpolant(xy(:,1),xy(:,2),Uxy(:,2),'natural');
    F3 = scatteredInterpolant(xy(:,1),xy(:,2),Uxy(:,3),'natural');
    nodes(:,4) = F1(nodes(:,2),nodes(:,3));
    nodes(:,5) = F2(nodes(:,2),nodes(:,3));
    nodes(:,6) = F3(nodes(:,2),nodes(:,3));
    
    X  = [nodes(elem(:,2),2) nodes(elem(:,3),2) nodes(elem(:,4),2) nodes(elem(:,5),2)]';
    Y  = [nodes(elem(:,2),3) nodes(elem(:,3),3) nodes(elem(:,4),3) nodes(elem(:,5),3)]';
    E1 = [nodes(elem(:,2),4) nodes(elem(:,3),4) nodes(elem(:,4),4) nodes(elem(:,5),4)]';
    E2 = [nodes(elem(:,2),5) nodes(elem(:,3),5) nodes(elem(:,4),5) nodes(elem(:,5),5)]';
    E3 = [nodes(elem(:,2),6) nodes(elem(:,3),6) nodes(elem(:,4),6) nodes(elem(:,5),6)]';
    nodes = nodes(:,1:3);
   
    Z = [];     E4 = [];    E5 = [];    E6 = [];
    
elseif DIM == 3 % 3d
    xyz     = rawData(:,1:3);
    Uxy    = rawData(:,4:9);
    [~,order] = sortrows(xyz,[1,3,2]);
    Uxy = Uxy(order,:);
     
    % get the number of vox along x, y, z
    xVox = size(unique(xyz(:,1)),1);
    yVox = size(unique(xyz(:,2)),1);
    zVox = size(unique(xyz(:,3)),1);

    % get step size and boundaries for the mesh
    step   = zeros(1,3);
    boxmin = zeros(1,3);
    boxmax = zeros(1,3);
    for ii = 1:3
        step(ii) = min(diff(unique(xyz(:,ii))));
        boxmin(ii) = xyz(1,ii)-step(ii)/2;
        boxmax(ii) = xyz(end,ii)+step(ii)/2;
    end

    %% Generate 3D mesh
    % generate nodes 
    [x,y,z] = meshgrid(boxmin(1):step(1):boxmax(1),boxmin(2):step(2):boxmax(2),boxmin(3):step(3):boxmax(3));
    numNodes = numel(x);
    coord = [reshape(x,numNodes,1), reshape(y,numNodes,1), reshape(z,numNodes,1)];
    nodes = [(1:numNodes)', sortrows(coord,[1,3,2])];

    % allocate array for elements
    elem = zeros(size(xyz,1),9);
    count = 1;

    % start loop over voxel dimensions
    for ix = 1:xVox
        for iz = 1:zVox
            for iy = 1:yVox

                % get element label
                elem(count,1) = count;

                % nodes on the plane with lower x
                elem(count,2) = iy + (iz-1)*(yVox+1) + (ix-1)*(yVox+1)*(zVox+1);
                elem(count,3) = elem(count,2) + 1;
                elem(count,4) = elem(count,3) + yVox + 1;
                elem(count,5) = elem(count,2) + yVox + 1;

                % nodes on the plane with higher x
                elem(count,6) = iy + (iz-1)*(yVox+1) + ix*(yVox+1)*(zVox+1);
                elem(count,7) = elem(count,6) + 1;
                elem(count,8) = elem(count,7) + yVox + 1;
                elem(count,9) = elem(count,6) + yVox + 1;

                count = count+1;
            end
        end
    end
    
    F1 = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),Uxy(:,1),'natural');
    nodes(:,5) = F1(nodes(:,2),nodes(:,3),nodes(:,4));
    F2 = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),Uxy(:,2),'natural');
    nodes(:,6) = F2(nodes(:,2),nodes(:,3),nodes(:,4));
    F3 = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),Uxy(:,3),'natural');
    nodes(:,7) = F3(nodes(:,2),nodes(:,3),nodes(:,4));
    F4 = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),Uxy(:,4),'natural');
    nodes(:,8) = F4(nodes(:,2),nodes(:,3),nodes(:,4));
    F5 = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),Uxy(:,5),'natural');
    nodes(:,9) = F5(nodes(:,2),nodes(:,3),nodes(:,4));
    F6 = scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),Uxy(:,6),'natural');
    nodes(:,10)= F6(nodes(:,2),nodes(:,3),nodes(:,4));
    
    X  = [nodes(elem(:,2),2) nodes(elem(:,3),2) nodes(elem(:,4),2) nodes(elem(:,5),2)...
          nodes(elem(:,6),2) nodes(elem(:,7),2) nodes(elem(:,8),2) nodes(elem(:,9),2)]';
    Y  = [nodes(elem(:,2),3) nodes(elem(:,3),3) nodes(elem(:,4),3) nodes(elem(:,5),3)...
          nodes(elem(:,6),3) nodes(elem(:,7),3) nodes(elem(:,8),3) nodes(elem(:,9),3)]';
    Z  = [nodes(elem(:,2),4) nodes(elem(:,3),4) nodes(elem(:,4),4) nodes(elem(:,5),4)...
          nodes(elem(:,6),4) nodes(elem(:,7),4) nodes(elem(:,8),4) nodes(elem(:,9),4)]';
    E1 = [nodes(elem(:,2),5) nodes(elem(:,3),5) nodes(elem(:,4),5) nodes(elem(:,5),5)...
          nodes(elem(:,6),5) nodes(elem(:,7),5) nodes(elem(:,8),5) nodes(elem(:,9),5)]';
    E2 = [nodes(elem(:,2),6) nodes(elem(:,3),6) nodes(elem(:,4),6) nodes(elem(:,5),6)...
          nodes(elem(:,6),6) nodes(elem(:,7),6) nodes(elem(:,8),6) nodes(elem(:,9),6)]';
    E3 = [nodes(elem(:,2),7) nodes(elem(:,3),7) nodes(elem(:,4),7) nodes(elem(:,5),7)...
          nodes(elem(:,6),7) nodes(elem(:,7),7) nodes(elem(:,8),7) nodes(elem(:,9),7)]';
    E4 = [nodes(elem(:,2),8) nodes(elem(:,3),8) nodes(elem(:,4),8) nodes(elem(:,5),8)...
          nodes(elem(:,6),8) nodes(elem(:,7),8) nodes(elem(:,8),8) nodes(elem(:,9),8)]';
    E5 = [nodes(elem(:,2),9) nodes(elem(:,3),9) nodes(elem(:,4),9) nodes(elem(:,5),9)...
          nodes(elem(:,6),9) nodes(elem(:,7),9) nodes(elem(:,8),9) nodes(elem(:,9),9)]';
    E6 = [nodes(elem(:,2),10) nodes(elem(:,3),10) nodes(elem(:,4),10) nodes(elem(:,5),10)...
          nodes(elem(:,6),10) nodes(elem(:,7),10) nodes(elem(:,8),10) nodes(elem(:,9),10)]';
    nodes = nodes(:,1:4);
end