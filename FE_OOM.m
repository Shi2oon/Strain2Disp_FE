function [X1,Y1,Ux,Uy,M4,Z1,Uz] = FE_OOM(alldata,ShapeFunOrder,resultsDir)
% last update 24/12/2020.

% type of shape function 'Linear' for Linear 'Quadratic' Elements

% include example data here
% this is matlab function that integrate the strain field to displacement
% either in 2D or 3D. the input should be the a matrix with each column
% vector containing a strain field components arranged as below
% alldata = [X(:) Y(:) Z(:) e11(:) e22(:) e33(:) e12(:) e13(:) e23(:)];
% then you need to decide the between a 

set(0,'defaultAxesFontSize',25);        set(0,'DefaultLineMarkerSize',14)
mkdir(resultsDir);

%% Create the mesh object:
fprintf ('Started Meshing ... ');
if size(alldata,2) == 9
    % alldata = [X(:) Y(:) Z(:) e11(:) e22(:) e33(:) e12(:) e13(:) e23(:)];
    NDIM  = 3;              % 3D
    NNODE = 8;
    NGP   = 8;
    [M4.Nodes,M4.Elements,M4.Xall,M4.Yall,M4.Zall,M4.E11,M4.E22,M4.E33,M4.E12,...
        M4.E13,M4.E23,M4.ScaleYN] = HexMeshAbaqus(alldata); % non unifrom data
    resultsDir = [resultsDir '\3D'];
    saveas(gcf,[resultsDir '_Meshed ' num2str(NGP) '.fig']);
    saveas(gcf,[resultsDir '_Meshed ' num2str(NGP) '.png']); close
else
    % alldata = [X(:) Y(:) e11(:) e22(:) e12(:)];
    NDIM          = 2;      % 2D
    alldata = sortrows(alldata,[1,2]);
    %     try
    %     catch
    %     fprintf('FE_Mesh_Generator failed ! .. Meshing ');
    [M4.Xall,M4.Yall,M4.E11,M4.E22,M4.E12,M4.Nodes,M4.Elements] = ...
        Meshing(alldata,NDIM,'Strain');  % support both unifrom and non unifrom data
        NGP = 4;        NNODE = 4;  % full integration, change NGP to 1 for reduced
        
    if strcmpi(ShapeFunOrder,'Quadratic') || strcmpi(ShapeFunOrder,'Linear')
        try% full integration, change NGP to 4 or 1 for reduced
        [NGP,NNODE,M4.Xall,M4.Yall,M4.E11,M4.E22,M4.E12,M4.Nodes,M4.Elements] = ...
            FE_Mesh_Generator(alldata,NDIM,ShapeFunOrder);
        catch err
            disp(err.message)
            disp('The current Quadratic meshing algorithm does not support non-rectnagular maps')
            disp('The shape function is now linear, 4 Gauss and nodes')
            ShapeFunOrder = 'Linear';
        end
    end
    resultsDir = [resultsDir '\' ShapeFunOrder];
    Z1 = []; Uz = [];
    [M4.Elements,M4.E11,M4.E22,M4.E12,~,~,M4.ScaleYN] = ... % remove outliers
        WhatNaN(M4.Elements,M4.E11,M4.E22,M4.E12,M4.Xall,M4.Yall);
end

%% Mesh
myMesh = Mesh_ME(NDIM,NNODE,NGP,ShapeFunOrder,M4.Nodes,M4.Elements); % Mesh
close all;          [M4.X,M4.Y] = myMesh.Plot_Mesh();   M4.IX = myMesh.IX';
if M4.ScaleYN == 'Y' && NDIM == 2
    addScale(1,[M4.X(:) M4.Y(:)]);    % Plot the mesh
end
% saveas(gcf,[resultsDir '_Meshed ' num2str(NGP) '.fig']);
% saveas(gcf,[resultsDir '_Meshed ' num2str(NGP) '.png']);
close all;  fprintf ('Done\nStarted Jacobian Calculation ..  ');
[~,dNdXGPS,~,~] = myMesh.Jacobian();    % Calculate the Jacobian
fprintf('Done for %d elements\nStart Assembly ...',length(M4.Elements));

%% The gradient:
DUDX = zeros(3,myMesh.NGP,myMesh.NELEMENTS);
for j=1:myMesh.NGP
    DUDX(1,j,:) = M4.E11(j,:);
    DUDX(2,j,:) = M4.E22(j,:);
    if NDIM == 3
        DUDX(3,j,:) = M4.E33(j,:);
        DUDX(4,j,:) = M4.E12(j,:);
        DUDX(5,j,:) = M4.E13(j,:);
        DUDX(6,j,:) = M4.E23(j,:);
    else
        DUDX(3,j,:) = M4.E12(j,:);
    end
end
% Create the data object:
myData = DataUm(NDIM,NNODE,NGP,'Strain',myMesh.NODES,...
                myMesh.ELEMENTS,myMesh.IX,dNdXGPS,DUDX);
[M4.K,M4.F]  = myData.Assembely();  % The stiffness matrix and the force vector:
save([resultsDir '_Integrated_Uxy.mat'],'M4','alldata','myMesh','NDIM',...
                                        'NNODE','NGP','dNdXGPS','DUDX','resultsDir');
% tic;[K,F]=Assembely(NDIM,NNODE,NGP,myMesh.NODES,myMesh.ELEMENTS,myMesh.IX,dNdXGPS,DUDX);toc
%% Determine the displacement:
fprintf(' Done\nStart Solving for Displacement ... \n');
mySolver = Dis_Solver(NDIM,NNODE,NGP,myMesh.NODES,myMesh.ELEMENTS,...
    myMesh.XNODES,myMesh.IX,M4.K,M4.F,M4.E11,M4.E22,M4.E12);
[M4.Ux,M4.Uy,M4.X,M4.Y,M4.Uz,M4.Z]   = mySolver.Plot_U();
save([resultsDir '_Integrated_Uxy.mat'],'M4','-append');
if M4.ScaleYN == 'Y' && NDIM == 2
    addScale([2 3 5],[M4.X(:) M4.Y(:)]);
end
% saveas(gcf,[resultsDir '_Displacement_GP_' num2str(NGP) '.tif'],'tiffn');
% saveas(gcf,[resultsDir '_Displacement_GP_' num2str(NGP) '.fig']);  close all
close
%% asign values to Nodes
for iV=1:size(M4.X,2)
    for iO=1:size(M4.X,1) %nodes
        if NDIM == 3
            for iJ=1:size(M4.X,3) %nodes
                A = ismember(M4.Nodes(:,2:4),[M4.X(iO,iV,iJ) M4.Y(iO,iV,iJ) M4.Z(iO,iV,iJ)],'row');
                NodeNumber = M4.Nodes(A,1);
                %                 OutData(NodeNumber,:) = [M4.X(iO,iV,iJ)  M4.Y(iO,iV,iJ)  M4.Z(iO,iV,iJ) ...
                %                                          M4.Ux(iO,iV,iJ) M4.Uy(iO,iV,iJ) M4.Uz(iO,iV,iJ)];
                M4.StepBCs((NodeNumber-1)*3+1,:) = [NodeNumber 1 1  M4.Ux(iO,iV,iJ)];
                M4.StepBCs((NodeNumber-1)*3+2,:) = [NodeNumber 2 2  M4.Uy(iO,iV,iJ)];
                M4.StepBCs((NodeNumber-1)*3+3,:) = [NodeNumber 3 3  M4.Uz(iO,iV,iJ)];
            end
            OutData = [M4.X(:) M4.Y(:) M4.Z(:) M4.Ux(:) M4.Uy(:) M4.Uz(:)];
        elseif NDIM == 2
            %             A = ismember(M4.Nodes(:,2:3),[M4.X(iO,iV) M4.Y(iO,iV)],'row');
            %             NodeNumber = M4.Nodes(A,1);
            %             OutData(NodeNumber,:) = [M4.X(iO,iV)  M4.Y(iO,iV) ...
            %                                      M4.Ux(iO,iV) M4.Uy(iO,iV)];
            OutData = [M4.X(:) M4.Y(:) M4.Ux(:) M4.Uy(:)];
        end
    end
end
save([resultsDir '_Integrated_Uxy.mat'],'OutData','M4','-append');

%% revert back to a matrix form, plot and remove rotation
if NDIM == 3
    tmp = sortrows(OutData,[3,1,2]);
    [~,dataum ] = reshapeData(tmp);
    X1 = dataum.X1;    Y1 = dataum.Y1;         Z1 = dataum.Z1;
    Ux = dataum.Ux;    Uy = dataum.Uy;         Uz = dataum.Uz;
    
    % remove rotation
    inpFile = fopen([fileparts(resultsDir) '\3D Data Uxy.dat'],'wt');
    fprintf(inpFile, 'TITLE = "3D Data Uxy"\n');
    fprintf(inpFile, 'VARIABLES = "x", "y", "z", "Vx", "Vy", "Vz", "isValid"\n');
    l = [length(unique(X1(:))),length(unique(Y1(:))),length(unique(Z1(:))),];
    fprintf(inpFile, 'ZONE T="Frame 0", I=%d, J=%d, K=%d, F=POINT\n',l(1),l(2),l(3));
    fprintf(inpFile,'%f %f %f %f %f %f %f\n',[X1(:)';Y1(:)';Z1(:)'-min(unique(Z1(:)));...
        Ux(:)';Uy(:)';Uz(:)';~isnan(X1(:))']);
    fclose(inpFile);
    %{
    try
        [M4.RadEulerAng,M4.rotCentre,tmp] = shoemake_3D_v04_07_02_Abdo...
            (fileparts(resultsDir),'3D Data Uxy.dat','um');
        
        % arranged corrected data
        tmp = sortrows(tmp,[3,1,2]);
        [~,dataum ] = reshapeData(tmp);
        X1 = dataum.X1;    Y1 = dataum.Y1;     Z1 = dataum.Z1;
        Ux = dataum.Ux;    Uy = dataum.Uy;     Uz = dataum.Uz;
    end
    if length(unique(Z1(:)))~=2 % a memberane
        Plot3D(sqrt(Ux.^2+Uy.^2+Uz.^2),Y1,X1,Z1,'\mum','U_{mag}');
        axis tight; axis image
        saveas(gcf,[resultsDir '_Umag_GP_' num2str(NGP) '_Raw.tif'],'tiffn');
        saveas(gcf,[resultsDir '_Umag_GP_' num2str(NGP) '_Raw.fig']);  close
        
        plotAllDis(Y1,X1,Z1,Ux,Uy,Uz,'\mum')
        saveas(gcf,[resultsDir '_Raw_DVC.fig']);
        saveas(gcf,[resultsDir '_Raw_DVC.tif']); close
        [~,ii] = reshapeStrainData( [M4.Xall(:) M4.Yall(:) M4.Zall(:) M4.E11(:) ...
                 M4.E22(:) M4.E33(:) M4.E12(:)  M4.E13(:)  M4.E23(:)]); %strain
        plotAllDis(ii.X1,ii.Y1,ii.Z1,ii.Exx,ii.Eyy,ii,'abs');
        saveas(gcf,[resultsDir '_Strain.fig']);
        saveas(gcf,[resultsDir '_Strain.tif']); close
    end
    %}
elseif NDIM == 2
    tmp = sortrows(OutData,[1,2]);
    [~,dataum ] = reshapeData(tmp);
    X1 = dataum.X1;     Y1 = dataum.Y1;	
    Ux = dataum.Ux;     Uy = dataum.Uy;
    
    % remove rotation
    [M4.theta,Ux,Uy,M4.rotCentre] = RotRemoving('true',X1, Y1, Ux, Uy);
    saveas(gcf,[resultsDir '_Corr_DIC.tif'],'tiffn');
    saveas(gcf,[resultsDir '_Corr_DIC.tif']); close
    
    subplot(1,1,1); contourf(X1,Y1,sqrt(Ux.^2+Uy.^2),'LineStyle','none');
    addScale([1 1 1],[M4.X(:) M4.Y(:)]);
    axis image;    axis xy;    colormap jet; box off;
    c = colorbar;        c.Label.String = 'U_{Mag}';
    set(gcf,'position',[10 50 669 489]);
    saveas(gcf,[resultsDir '_Umag_GP_' num2str(NGP) '.tif'],'tiffn');
    saveas(gcf,[resultsDir '_Umag_GP_' num2str(NGP) '.fig']);  close
end

save([resultsDir '_Integrated_Uxy.mat'],'Ux','Uy','Uz','X1','Y1','Z1','M4','-append');
fprintf(' Done\n');
end