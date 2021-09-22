function [NGP,NNode,Xa,Ya,E11,E22,E12,Nodes,Elements] = ...
                FE_Mesh_Generator(alldata,NDim,ShapeFunOrder,a)
% alldata =[X Y E11 E22 E12] 2D Centroid
% alldata = [X Y Z E11 E22 E33 E12 E13 E23] 3D Centroid
% a is the crack length in X

%% Domain Dimensions
% The size:
X = unique(alldata(:,1));   Lex = X(2)-X(1);     X0 = min(X)-0.5*Lex;   
Lx     = max(unique(alldata(:,1)))-min(unique(alldata(:,1)))+Lex;
Ndivx  = length(unique(alldata(:,1)));

Y = unique(alldata(:,2));  Ley = Y(2)-Y(1);     Y0 = min(Y)-0.5*Ley;  
Ly     = max(unique(alldata(:,2)))-min(unique(alldata(:,2)))+Ley;
Ndivy  = length(unique(alldata(:,2)));

%% Linear or Quadratic Integration and Mesh parameters:
if  strcmp(ShapeFunOrder,'Linear') % NGP 1 reduced, 4 full
    NGP = 4;            NNode = 4;  
    if NGP == 1
        dX = [0];                   dY = [0];               
	elseif NGP == 4
        dX = (Lex/4)*[-1 1 -1 1];   dY = (Ley/4)*[-1 -1 1 1];  
    end
    % 2D , add z for 3d
    Lex = 2.0*Lex;      Ley = 2.0*Ley;      N = 1;
    if mod(Ndivx,2) == 0
        Ndivx  = Ndivx/2;       
    else
        Ndivx  = (Ndivx-1)/2;      
    end
    if mod(Ndivy,2) == 0
        Ndivy  = Ndivy/2;
    else
        Ndivy  = (Ndivy-1)/2;      
    end
elseif  strcmp(ShapeFunOrder,'Quadratic') % NGP 4 reduced, 9 full, 1 
    NGP = 9;                    NNode = 8;          N = 2;
    Lex = 3.0*Lex;              Ley = 3.0*Ley;
    Ndivx = floor(Ndivx/3);
    Ndivy = floor(Ndivy/3); 
	if NGP == 1
        dX = [0];                       dY = [0];               
	elseif NGP == 4
        dX = (Lex/4)*[-1  1 -1 1];      dY = (Ley/4)*[-1 -1  1 1];     
	elseif NGP == 9
        dX = (Lex/4)*[-1  0  1 -1 0 1 -1 0 1];
        dY = (Ley/4)*[-1 -1 -1  0 0 0  1 1 1];
	end
end
Tol = min([Lex Ley])/NNode;
NElements = Ndivx*Ndivy;

%% 3D
if NDim==3
    Z = unique(alldata(:,3));   
    Lez = Z(2)-Z(1);    
    Z0 = min(Z)-0.5*Lez; 
    Lz = max(unique(alldata(:,3)))-min(unique(alldata(:,3)))+Lez;
    Ndivz = length(unique(alldata(:,3)));
    % NGP 8 full, 1 reduced
    NGP = 8;           NNode = 8;           
    Lez = 2.0*Lez;
    if mod(Ndivx,2) == 0
        Ndivz  = Ndivz/2;
    else
        Ndivz  = (Ndivz-1)/2;      
    end
    if NGP == 1
        dX = [0];           dY = [0];           dZ = [0];
	elseif NGP == 8
        dX = (Lex/4)*[-1  1 -1  1 -1  1 -1 1];
        dY = (Ley/4)*[-1 -1  1  1 -1 -1  1 1];
        dZ = (Lez/4)*[-1 -1 -1 -1  1  1  1 1]; 
	end
    Tol = min([Lex Ley Lez])/NNode;
    NElements = Ndivx*Ndivy*Ndivz;
end

%% Discretisation:
% Nodes:
NNodes = (N*Ndivx+1)*(N*Ndivy+1);
if NDim==3
    NNodes = NNodes*(N*Ndivz+1);
end
Nodes  = zeros(NNodes,3);
for i=1:(N*Ndivx+1)   
    for j=1:(N*Ndivy+1) 
        if NDim == 2
            Nodes(i+(j-1)*(N*Ndivy+1),1) = i+(j-1)*(N*Ndivy+1);
            Nodes(i+(j-1)*(N*Ndivy+1),2) = X0+(i-1)*(Lex/N);
            Nodes(i+(j-1)*(N*Ndivy+1),3) = Y0+(j-1)*(Ley/N);
        elseif NDim == 3
            for o=1:(N*Ndivz+1) 
                Nodes(i+(j-1)*(N*Ndivy+1)+(o-1)*(N*Ndivz+1),1) = i+(j-1)*(N*Ndivy+1)+(o-1)*(N*Ndivz+1);
                Nodes(i+(j-1)*(N*Ndivy+1)+(o-1)*(N*Ndivz+1),2) = X0+(i-1)*(Lex/N);
                Nodes(i+(j-1)*(N*Ndivy+1)+(o-1)*(N*Ndivz+1),3) = Y0+(j-1)*(Ley/N);
                Nodes(i+(j-1)*(N*Ndivy+1)+(o-1)*(N*Ndivz+1),4) = Z0+(j-1)*(Lez/N);
            end
        end
    end
end

%% Elements:
Elements  = zeros(NElements,NNode+1);
NodesNill = [];
k         = 1;  countr=0;
if NDim == 2
  for i=1:Ndivx   
    for j=1:Ndivy   
        if strcmp(ShapeFunOrder,'Linear')
            countr = countr+1;
            Elements(i+(j-1)*Ndivy,1) = countr;
            Elements(i+(j-1)*Ndivy,2) = Nodes(i+(j-1)*(Ndivy+1),1); 
            Elements(i+(j-1)*Ndivy,3) = Nodes(i+1+(j-1)*(Ndivy+1),1);
            Elements(i+(j-1)*Ndivy,4) = Nodes(i+1+(j)*(Ndivy+1),1); 
            Elements(i+(j-1)*Ndivy,5) = Nodes(i+(j)*(Ndivy+1),1);  
            % Find element:
            % Output strain data:
            for l = 1:NGP
                Xc = X0+(i-0.5)*Lex+dX(l);      
                [~,Indx] = min(abs(X-Xc));      Xc = X(Indx);
                Yc = Y0+(j-0.5)*Ley+dY(l);
                [~,Indx] = min(abs(Y-Yc));      Yc = Y(Indx);
                
                A = ismember(alldata(:,1:2),[Xc Yc],'row');
                
                Xa(l,i+(j-1)*Ndivx)  = alldata(A,1);
                Ya(l,i+(j-1)*Ndivx)  = alldata(A,2);
                E11(l,i+(j-1)*Ndivx) = alldata(A,3);
                E22(l,i+(j-1)*Ndivx) = alldata(A,4);
                E12(l,i+(j-1)*Ndivx) = alldata(A,5);
            end                

            
        elseif strcmp(ShapeFunOrder,'Quadratic')
            countr = countr+1;
            Elements(i+(j-1)*Ndivx,1) = countr;
            Elements(i+(j-1)*Ndivx,2) = Nodes(N*i-1+(N*j-1-1)*(N*Ndivy+1),1);
            Elements(i+(j-1)*Ndivx,3) = Nodes(N*i-1+(N*j-1-1)*(N*Ndivy+1)+2,1);
            Elements(i+(j-1)*Ndivx,4) = Nodes(N*i-1+(N*j+1-1)*(N*Ndivy+1)+2,1);  
            Elements(i+(j-1)*Ndivx,5) = Nodes(N*i-1+(N*j+1-1)*(N*Ndivy+1),1);  
            Elements(i+(j-1)*Ndivx,6) = Nodes(N*i-1+(N*j-1-1)*(N*Ndivy+1)+1,1); 
            Elements(i+(j-1)*Ndivx,7) = Nodes(N*i-1+(N*j-1)*(N*Ndivy+1)+2,1);  
            Elements(i+(j-1)*Ndivx,8) = Nodes(N*i-1+(N*j+1-1)*(N*Ndivy+1)+1,1);  
            Elements(i+(j-1)*Ndivx,9) = Nodes(N*i-1+(N*j-1)*(N*Ndivy+1),1);  
            NodesNill(k)              = Nodes(N*i-1+(N*j-1)*(N*Ndivy+1)+1,1);
            k                         = k+1;
            
            % Output strain data:
            for l = 1:NGP
                Xc = X0+(i-0.5)*Lex+dX(l);     
                [~,Indx] = min(abs(X-Xc));      Xc = X(Indx);
                Yc = Y0+(j-0.5)*Ley+dY(l) ;
                [~,Indx] = min(abs(Y-Yc));      Yc = Y(Indx);
                
                A = ismember(alldata(:,1:2),[Xc Yc],'row');
                
                Xa(l,i+(j-1)*Ndivx)  = alldata(A,1);
                Ya(l,i+(j-1)*Ndivx)  = alldata(A,2);
                E11(l,i+(j-1)*Ndivx) = alldata(A,3);
                E22(l,i+(j-1)*Ndivx) = alldata(A,4);
                E12(l,i+(j-1)*Ndivx) = alldata(A,5);
            end 
        end
    end    
  end
elseif NDim ==3
    %% 3D
  for i=1:Ndivx   
    for j=1:Ndivy 
        for o=1:Ndivz
            countr = countr+1;
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,1) = countr;
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,2) = Nodes(i+(j-1)*(Ndivy+1)+(o-1)*(Ndivy+1),1); 
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,3) = Nodes(i+1+(j-1)*(Ndivy+1)+(o-1)*(Ndivy+1),1);
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,4) = Nodes(i+1+j*(Ndivy+1)+(o-1)*(Ndivy+1),1); 
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,5) = Nodes(i+j*(Ndivy+1)+(o-1)*(Ndivy+1),1); 
            
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,6) = Nodes(i+j*(Ndivy+1)+1+(o-1)*(Ndivy+1),1); 
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,7) = Nodes(i+j*(Ndivy+1)+1+o*(Ndivy+1),1); 
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,8) = Nodes(i+j*(Ndivy+1)+1+o*(Ndivy+1),1); 
            Elements(i+(j-1)*Ndivy+(o-1)*Ndivz,9) = Nodes(i+j*(Ndivy+1)+o*(Ndivy+1),1); 
            % Find element:
            % Output strain data:
            for l = 1:NGP
                Xc = X0+(i-0.5)*Lex+dX(l);      
                [~,Indx] = min(abs(X-Xc));      Xc = X(Indx);
                Yc = Y0+(j-0.5)*Ley+dY(l);
                [~,Indx] = min(abs(Y-Yc));      Yc = Y(Indx);
                Zc = Z0+(o-0.5)*Lez+dZ(l);
                [~,Indx] = min(abs(Z-Zc));      Zc = Z(Indx);
                
                A = ismember(alldata(:,1:2),[Xc Yc],'row');
                
                Xa(l,i+(j-1)*Ndivx)  = alldata(A,1);
                Ya(l,i+(j-1)*Ndivx)  = alldata(A,2);
                E11(l,i+(j-1)*Ndivx) = alldata(A,3);
                E22(l,i+(j-1)*Ndivx) = alldata(A,4);
                E12(l,i+(j-1)*Ndivx) = alldata(A,5);
            end
        end   
    end
  end
end

[NodesUn,IX,~] = setxor(Nodes(:,1),NodesNill);

%%  Define the crack:
if exist('a','var')
     a = 0.5*Lx;

% determine the crack face nodes and the cracktip node: 
NodesCr  = [];
j        = 1;
for i=1:NNodes
   if  (Nodes(i,2)>=-Tol) &&  (Nodes(i,2)<a-Tol)
       if (Nodes(i,3)>=0.5*Ly-Tol) &&  (Nodes(i,3)<=0.5*Ly+Tol)
           NodesCr(j,1) = Nodes(i,1);
           j         = j+1;       
       end
   elseif  (Nodes(i,2)>=a-Tol) &&  (Nodes(i,2)<=a+Tol)
       if (Nodes(i,3)>=0.5*Ly-Tol) &&  (Nodes(i,3)<=0.5*Ly+Tol)
           NodesTip = Nodes(i,1);
           j         = j+1;       
       end
   end
end

% The Elements centres:
Xc  = zeros(NElements,1);
Yc  = zeros(NElements,1);
for i=1:NElements
    for j=1:NNode
        Xc(i) = Xc(i)+Nodes(Elements(i,j),2);
        Yc(i) = Yc(i)+Nodes(Elements(i,j),3);
    end
    Xc(i) = Xc(i)/NNode;
    Yc(i) = Yc(i)/NNode;
end

% Determine the crack adjacent elements:
ElementsCrUp = [];
ElementsCrDo = [];
j            = 1;
k            = 1;
for i=1:NElements  
    if  (Yc(i)>=0.5*Ly-0.5*Ley-Tol) && (Yc(i)<=0.5*Ly-0.5*Ley+Tol)  
        if  (Xc(i)>=(j-0.5)*Lex-Tol) && (Xc(i)<=(j-0.5)*Lex+Tol) && (Xc(i)<=a-Tol) 
            ElementsCrDo(j,1) = i;    
            j                 = j+1;
        end
    elseif  (Yc(i)>=0.5*Ly+0.5*Ley-Tol) && (Yc(i)<=0.5*Ly+0.5*Ley+Tol) && (Xc(i)<=0.5*Ly-Tol) 
        if  (Xc(i)>=(k-0.5)*Lex-Tol) && (Xc(i)<=(k-0.5)*Lex+Tol) 
            ElementsCrUp(k,1) = i;    
            k                 = k+1;  
        end
    end
end

% Create the seam:
NodesCom = zeros(length(ElementsCrDo),N+1);
IXDo     = zeros(length(ElementsCrDo),N+1);
IXUp     = zeros(length(ElementsCrDo),N+1);
for i=1:length(ElementsCrDo) 
    [NodesCom(i,:),IXDo(i,:),IXUp(i,:)] = intersect(Elements(ElementsCrDo(i,1),:),Elements(ElementsCrUp(i,1),:)); 
end

% Duplicate the nodes:
NodesCrII = zeros(length(NodesCr),3);
for i=1:length(NodesCr)
    NodesCrII(i,1) = NNodes+i;   
    NodesCrII(i,2) = Nodes(NodesCr(i),2);   
    NodesCrII(i,3) = Nodes(NodesCr(i),3);         
end

% Add the duplicated nodes:
Nodes = [Nodes;NodesCrII];

% [NodesCom(i,:)
IDX  = zeros(length(ElementsCrUp),NNode);
for i=1:length(ElementsCrUp) 
    [~,IDX(i,:)] = ismember(Elements(ElementsCrUp(i,1),:),NodesCr);
    for j=1:NNode
        if IDX(i,j) ~= 0
            [~,k] = ismember(Elements(ElementsCrUp(i,1),j),NodesCr);
            Elements(ElementsCrUp(i,1),j) = NodesCrII(k,1);
        end
    end
end
end
%}

%% Plot Mesh:
% Create the coordinate matrix:        
X1 = zeros(NNode/N,NElements);
X2 = zeros(NNode/N,NElements);
for i=1:NElements
    for j=1:NNode/N
        X1(j,i) = Nodes(Elements(i,j+1),2);
        X2(j,i) = Nodes(Elements(i,j+1),3);
    end
end
figure();
fill(X1,X2,'w','EdgeColor','black','LineWidth',1.0);
hold on
plot(Nodes(NodesUn,2),Nodes(NodesUn,3),'.','MarkerFaceColor','black','MarkerEdgeColor','black')
if exist('a','var')
    plot(Nodes(NodesCr,2),Nodes(NodesCr,3),'.','MarkerFaceColor','red','MarkerEdgeColor','red')
    plot(Nodes(NodesTip,2),Nodes(NodesTip,3),'.','MarkerFaceColor','blue','MarkerEdgeColor','blue')
end
axis equal
axis tight
axis off
xlabel('$X_{1}$','Interpreter','LaTeX','FontSize',12,'FontName','Times New Roman');
ylabel('$X_{2}$','Interpreter','LaTeX','FontSize',12,'FontName','Times New Roman');
hold off

end


















