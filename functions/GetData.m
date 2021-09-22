%% a function to get the data back after solving 
function [x,y,Ux,Uy] = GetData(X,Y,U1,U2,DIM)  
if DIM == 2 
    %Guess the dimensions of the image
    mesh.winodow = [length(unique(X,'first')),length(unique(Y,'first'))];
    mesh.winFE   = [1 1; mesh.winodow(1) mesh.winodow(2)]; % take 2 elements fromc coners

    %% Mesh for nodes & elements
    %Number of elements in the FE mesh
    szelnFE = (mesh.winFE(2,2)-mesh.winFE(1,2))*(mesh.winFE(2,1)-mesh.winFE(1,1));
    count   = 0;
    %For each FE element, find constitutive DIC nodes
    elnFE = zeros(szelnFE,4);
    for c =  mesh.winFE(1,2):1:mesh.winFE(2,2)-1
        for r =  mesh.winFE(1,1):1:mesh.winFE(2,1)-1
        count = count + 1;
        elnFE(count,:) = ...
            [c+(r)*mesh.winodow(2), ...
            c+1+(r)*mesh.winodow(2), ...
            c+1+(r-1)*mesh.winodow(2), ...
            c+(r-1)*mesh.winodow(2)];
        end
    end

    %% Assemble element nodes for FE    
X  = X';   XX  = zeros(max(max(elnFE)),4);	XX(XX==0)   = NaN;	XX(1:length(X),:)   = X; 
Y  = Y';   YY  = zeros(max(max(elnFE)),4);	YY(YY==0)   = NaN;	YY(1:length(Y),:)   = Y;
U1 = U1';  UU1 = zeros(max(max(elnFE)),4); 	UU1(UU1==0) = NaN;	UU1(1:length(U1),:) = U1; 
U2 = U2';  UU2 = zeros(max(max(elnFE)),4);	UU2(UU2==0) = NaN;	UU2(1:length(U2),:) = U2; 
    for i=1:4
        x(elnFE(:,i))  = XX(elnFE(:,i),i);       % nodes corrdinates in X
        y(elnFE(:,i))  = YY(elnFE(:,i),i);       % nodes coordinates in Y
        Ux(elnFE(:,i)) = UU1(elnFE(:,i),i);      % nodes coordinates in Data 1
        Uy(elnFE(:,i)) = UU2(elnFE(:,i),i);      % nodes coordinates in Data 2
    end
    %%
x(isnan(x)) = [];	y(isnan(y)) = [];   Ux(isnan(Ux)) = [];      Uy(isnan(Uy)) = []; 
try
    alldata = [x(:) y(:) Ux(:) Uy(:)];
    alldata = reshapeData(alldata);
    x = alldata.X;  y = alldata.Y;  Ux = alldata.Ux;    Uy = alldata.Uy;
catch err
    warning(err.message);
    % do nothing, report data as it's
end

% close all;      imagesc(alldata.X(1,:),alldata.Y(:,1),alldata.Uy);    
% c=colorbar;       set(gca,'Ydir','normal');           axis image;
% title('Answer in the command line');
% xlabel('X [m]','FontSize',20,'FontName','Times New Roman');          
% ylabel('Y [m]','FontSize',20,'FontName','Times New Roman');
% c.Label.String = 'U_Y';
% set(gcf,'position',[30 50 1300 950]); 
end