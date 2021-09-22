clear all
close all
clc
format long

%% Domain Dimensions:
% The geometry:
NDim      = 2;

% The size:
Lx        = 10;
Ly        = 10;
Tol       = 1.0E-3;

% Mesh parameters:
EleFlag   = '2nd';
if strcmp(EleFlag,'1st') == 1
    NNode = 4;
elseif strcmp(EleFlag,'2nd') == 1
    NNode = 8;
end

Ndivx     = 20;
Ndivy     = 20;
Lex       = Lx/Ndivx;
Ley       = Ly/Ndivy;

%% Discretisation:
% Nodes:
if strcmp(EleFlag,'1st') == 1
    N = 1;
elseif strcmp(EleFlag,'2nd') == 1
    N = 2;
end

NNodes = (N*Ndivx+1)*(N*Ndivy+1);

Nodes  = zeros(NNodes,3);
for i=1:(N*Ndivx+1)   
    for j=1:(N*Ndivy+1) 
        Nodes(i+(j-1)*(N*Ndivx+1),1) = i+(j-1)*(N*Ndivx+1);
        Nodes(i+(j-1)*(N*Ndivx+1),2) = (i-1)*(Lex/N);
        Nodes(i+(j-1)*(N*Ndivx+1),3) = (j-1)*(Ley/N);
    end
end

% Elements:
NElements = Ndivx*Ndivy;
Elements  = zeros(NElements,NNode);

NodesNill = [];
k         = 1;
for i=1:Ndivx   
    for j=1:Ndivy   
        if strcmp(EleFlag,'1st') == 1
            Elements(i+(j-1)*Ndivx,1) = Nodes(i+(j-1)*(Ndivx+1),1);
            Elements(i+(j-1)*Ndivx,2) = Nodes(i+(j-1)*(Ndivx+1)+1,1);
            Elements(i+(j-1)*Ndivx,3) = Nodes(i+(j)*(Ndivx+1)+1,1);  
            Elements(i+(j-1)*Ndivx,4) = Nodes(i+(j)*(Ndivx+1),1);    
        elseif strcmp(EleFlag,'2nd') == 1  
            Elements(i+(j-1)*Ndivx,1) = Nodes(N*i-1+(N*j-1-1)*(N*Ndivx+1),1);
            Elements(i+(j-1)*Ndivx,2) = Nodes(N*i-1+(N*j-1-1)*(N*Ndivx+1)+2,1);
            Elements(i+(j-1)*Ndivx,3) = Nodes(N*i-1+(N*j+1-1)*(N*Ndivx+1)+2,1);  
            Elements(i+(j-1)*Ndivx,4) = Nodes(N*i-1+(N*j+1-1)*(N*Ndivx+1),1);  
            Elements(i+(j-1)*Ndivx,5) = Nodes(N*i-1+(N*j-1-1)*(N*Ndivx+1)+1,1); 
            Elements(i+(j-1)*Ndivx,6) = Nodes(N*i-1+(N*j-1)*(N*Ndivx+1)+2,1);  
            Elements(i+(j-1)*Ndivx,7) = Nodes(N*i-1+(N*j+1-1)*(N*Ndivx+1)+1,1);  
            Elements(i+(j-1)*Ndivx,8) = Nodes(N*i-1+(N*j-1)*(N*Ndivx+1),1);  
            NodesNill(k)              = Nodes(N*i-1+(N*j-1)*(N*Ndivx+1)+1,1);
            k                         = k+1;
        end
    end    
end

[NodesUn,IX,~] = setxor(Nodes(:,1),NodesNill); 

% Define the crack:
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

%% Plot Mesh:
% Create the coordinate matrix:        
X1 = zeros(NNode/N,NElements);
X2 = zeros(NNode/N,NElements);
for i=1:NElements
    for j=1:NNode/N
        X1(j,i) = Nodes(Elements(i,j),2);
        X2(j,i) = Nodes(Elements(i,j),3);
    end
end
figure();
fill(X1,X2,'w','EdgeColor','black','LineWidth',1.0);
hold on
plot(Nodes(NodesUn,2),Nodes(NodesUn,3),'o','MarkerFaceColor','black','MarkerEdgeColor','black')
hold on
plot(Nodes(NodesCr,2),Nodes(NodesCr,3),'o','MarkerFaceColor','red','MarkerEdgeColor','red')
hold on
plot(Nodes(NodesTip,2),Nodes(NodesTip,3),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue')
hold on
axis equal
axis tight
axis off
xlabel('$X_{1}$','Interpreter','LaTeX','FontSize',12,'FontName','Times New Roman');
ylabel('$X_{2}$','Interpreter','LaTeX','FontSize',12,'FontName','Times New Roman');
hold off

























