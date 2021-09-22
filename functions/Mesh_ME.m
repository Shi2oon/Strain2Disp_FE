classdef Mesh_ME
    properties
        NDIM
        NNODE
        NGP
        ShapeFunOrder
        NODES
        ELEMENTS
        NNODES
        NELEMENTS
        XiN
        XiGP
        IX
        XNODES
    end
    
    methods
        % Constructor:
        function M = Mesh_ME(NDIM,NNODE,NGP,ShapeFunOrder,NODES,ELEMENTS)
            M.NDIM          = NDIM;
            M.NNODE         = NNODE;
            M.NGP           = NGP;
            M.ShapeFunOrder = ShapeFunOrder;
            M.NODES         = NODES;
            M.ELEMENTS      = ELEMENTS;
        end
        
        %% The number of nodes in the mesh:
        function NNODES = get.NNODES(M)
            NNODES = length(M.NODES(:,1));
        end
        
        %% The number of elements in the mesh:
        function NELEMENTS = get.NELEMENTS(M)
            NELEMENTS = length(M.ELEMENTS(:,1));
        end
        
        %% The local coordinates of the nodes:
        function XiN = get.XiN(M)
            if M.NDIM == 2 			% 2D:
                if M.NNODE == 4
                    XiN = [-1,-1; 1,-1 ;1,1 ;-1,1];
                elseif M.NNODE == 8
                    XiN = [-1,-1; 1,-1 ;1,1;-1,1;0,-1;1,0;0,1;-1,0];
                end
                
            elseif M.NDIM == 3  	% 3D:
                XiN = [-1,-1,-1; 1,-1,-1;1,1,-1;-1,1,-1;1,1,1;1,-1,1;1,1,1;-1,1,1];
            end
        end
        
        %% The local coordinates of the Gauss points:
        function XiGP = get.XiGP(M)
            if M.NDIM == 2        	% 2D:
                if M.NNODE == 4
                    if M.NGP == 1
                        XiGP = [0.0,0.0];
                    else%if M.NGP == 4
                        XiGP = 0.5*[-1,-1;	1,-1;   -1,1 ; 	1,1];
                    end
                elseif M.NNODE == 8
                    if M.NGP == 1
                        XiGP = [0.0,0.0];
                    elseif M.NGP > 1 && M.NGP <=4
                        XiGP = 0.5*[ -1,-1;	1,-1;   -1,1; 	1,1];
                    elseif M.NGP > 4
                        XiGP = 0.5*[-1,-1;	0,-1;	1,-1;   -1,0;
                            0,0;	1,0;    -1,1;	0,1;	1,1];
                    end
                end
                
            elseif M.NDIM == 3  	% 3D:
                if M.NGP == 1
                    XiGP = [0.0,0.0,0.0];
                else%if M.NGP == 8
                    XiGP = 0.5*[-1,-1,-1;	1,-1,-1;    -1,1,-1;	1,1,-1;
                        -1,-1,1;   	1,-1,1;     -1,1,1;    	1,1,1];
                end
            end
        end
        
        %% The assembely matrix:
        function IX = get.IX(M)
            [~,IX] = ismember(M.ELEMENTS(:,2:end),M.NODES);
            IX     = IX';
        end
        
        %% The physical coordinates of the nodes:
        function XNODES = get.XNODES(M)
            for j=1:M.NNODE
                for l=1:M.NDIM
                    XNODES(:,M.NDIM*(j-1)+l) = M.NODES(M.IX(j,:)',l+1);
                end
            end
            XNODES = XNODES';
        end
        %% The Jacobian and shape functions derivatives at the nodes and Gauss points:
        function [dNdXNODES,dNdXGPS,dXdXiNODES,dXdXiGPS] = Jacobian(M)
            % 1. Determine the shape function derivative matrix in the local coordinate:
            % 1.1. In the nodes:
            dNdXiNODESV = zeros(M.NDIM,M.NNODE,M.NNODE);
            dNdXiNODES  = zeros(M.NDIM*M.NDIM,M.NNODE*M.NDIM,M.NNODE);
            for i=1:M.NNODE
                if M.NDIM == 2
                    [~,~,dNdXiNODES(:,:,i),dNdXiNODESV(:,:,i)] = ...
                        M.ShapeFunctions2D(M.NNODE,M.XiN(i,1),M.XiN(i,2));
                elseif M.NDIM == 3
                    [~,~,dNdXiNODES(:,:,i),dNdXiNODESV(:,:,i)] = ...
                        M.ShapeFunctions3D(M.NNODE,M.XiN(i,1),M.XiN(i,2),M.XiN(i,3));
                end
            end
            
            % 1.2. In the Gauss points
            dNdXiGPSV   = zeros(M.NDIM,M.NNODE,M.NGP);
            dNdXiGPS    = zeros(M.NDIM*M.NDIM,M.NNODE*M.NDIM,M.NGP);
            for i=1:M.NGP
                if M.NDIM == 2
                    [~,~,dNdXiGPS(:,:,i),dNdXiGPSV(:,:,i)] = ...
                        M.ShapeFunctions2D(M.NNODE,M.XiGP(i,1),M.XiGP(i,2));
                elseif M.NDIM == 3
                    [~,~,dNdXiGPS(:,:,i),dNdXiGPSV(:,:,i)] = ...
                        M.ShapeFunctions3D(M.NNODE,M.XiGP(i,1),M.XiGP(i,2),M.XiGP(i,3));
                end
            end
            
            % 2. Determine the Jacobian:
            % 2.1. Nodes
            dXdXiNODES  = zeros(M.NDIM,M.NDIM,M.NNODE,M.NELEMENTS);
            dNdXiXNODES = zeros(M.NDIM*M.NDIM,M.NNODE,M.NELEMENTS);
            for j=1:M.NNODE
                dNdXiXNODES(:,j,:)  = dNdXiNODES(:,:,j)*M.XNODES;
            end
            if M.NDIM == 2
                dXdXiNODES(1,1,:,:) = dNdXiXNODES(1,:,:);
                dXdXiNODES(1,2,:,:) = dNdXiXNODES(3,:,:);
                dXdXiNODES(2,1,:,:) = dNdXiXNODES(2,:,:);
                dXdXiNODES(2,2,:,:) = dNdXiXNODES(4,:,:);
            elseif M.NDIM == 3
                dXdXiNODES(1,1,:,:) = dNdXiXNODES(1,:,:);
                dXdXiNODES(1,2,:,:) = dNdXiXNODES(4,:,:);
                dXdXiNODES(1,3,:,:) = dNdXiXNODES(7,:,:);
                dXdXiNODES(2,1,:,:) = dNdXiXNODES(2,:,:);
                dXdXiNODES(2,2,:,:) = dNdXiXNODES(5,:,:);
                dXdXiNODES(2,3,:,:) = dNdXiXNODES(8,:,:);
                dXdXiNODES(3,1,:,:) = dNdXiXNODES(3,:,:);
                dXdXiNODES(3,2,:,:) = dNdXiXNODES(6,:,:);
                dXdXiNODES(3,3,:,:) = dNdXiXNODES(9,:,:);
            end
            
            % 2.2. Gaussian
            dXdXiGPS    = zeros(M.NDIM,M.NDIM,M.NGP,M.NELEMENTS);
            dNdXiXGPS   = zeros(M.NDIM*M.NDIM,M.NGP,M.NELEMENTS);
            for j=1:M.NGP
                dNdXiXGPS(:,j,:)  = dNdXiGPS(:,:,j)*M.XNODES;
            end
            if M.NDIM == 2
                dXdXiGPS(1,1,:,:) = dNdXiXGPS(1,:,:);
                dXdXiGPS(1,2,:,:) = dNdXiXGPS(3,:,:);
                dXdXiGPS(2,1,:,:) = dNdXiXGPS(2,:,:);
                dXdXiGPS(2,2,:,:) = dNdXiXGPS(4,:,:);
            elseif M.NDIM == 3
                dXdXiGPS(1,1,:,:) = dNdXiXGPS(1,:,:);
                dXdXiGPS(1,2,:,:) = dNdXiXGPS(4,:,:);
                dXdXiGPS(1,3,:,:) = dNdXiXGPS(7,:,:);
                dXdXiGPS(2,1,:,:) = dNdXiXGPS(2,:,:);
                dXdXiGPS(2,2,:,:) = dNdXiXGPS(5,:,:);
                dXdXiGPS(2,3,:,:) = dNdXiXGPS(8,:,:);
                dXdXiGPS(3,1,:,:) = dNdXiXGPS(3,:,:);
                dXdXiGPS(3,2,:,:) = dNdXiXGPS(6,:,:);
                dXdXiGPS(3,3,:,:) = dNdXiXGPS(9,:,:);
            end
            
            % dNdX for Nodes and Gaussian points
            dNdXNODES   = zeros(M.NDIM,M.NNODE,M.NNODE,M.NELEMENTS);
            dNdXGPS     = zeros(M.NDIM,M.NNODE,M.NGP,M.NELEMENTS);
            for i=1:M.NELEMENTS
                for j=1:M.NNODE
                    for k=1:M.NNODE
                        dNdXNODES(:,k,j,i) = dXdXiNODES(:,:,j,i)\dNdXiNODESV(:,k,j);
                    end
                end
                for l=1:M.NGP
                    for m=1:M.NNODE
                        dNdXGPS(:,m,l,i) = dXdXiGPS(:,:,l,i)\dNdXiGPSV(:,m,l);
                    end
                end
            end
        end
        
        %% Plot Mesh
        function [X1,X2,X3] = Plot_Mesh(M)
            txtit = {[' A ' num2str(M.NDIM) 'D Mesh']; ['w/ ' num2str(M.NNODE)...
                ' nodes/element, ' num2str(M.NGP) ' Gaussian points/element'];...
                ['and a ' M.ShapeFunOrder ' Shape function Order']};
            if M.NDIM == 2		% 2D Mesh
                % Multiplier:
                if M.NNODE == 4;		FIX = [1 2 3 4];
                elseif M.NNODE == 8;	FIX = [1 5 2 6 3 7 4 8];
                end
                % Create the coordinate matrix:
                for j=1:M.NNODE
                    X1(j,1:M.NELEMENTS) = M.XNODES(M.NDIM*(FIX(j)-1)+1,1:M.NELEMENTS);
                    X2(j,1:M.NELEMENTS) = M.XNODES(M.NDIM*(FIX(j)-1)+2,1:M.NELEMENTS);
                end
                X3 = [];
                figure();
                fill(X1,X2,'w');						hold on
                plot(M.NODES(:,2),M.NODES(:,3),'.','MarkerFaceColor','r',...
                    'MarkerEdgeColor','r');		hold off;
                xlabel('X[\mum]');	ylabel('Y[\mum]');	title(txtit);	axis image;
                set(gcf,'position',[10 50 1900 950]);
                
            elseif M.NDIM == 3
                % The facets:
                FIX = [1 2 3 4;5 8 7 6;1 5 6 2;2 6 7 3;3 7 8 4;4 8 5 1];
                figure();
                for j=1:6
                    for k=1:4
                        X1(j,k,1:M.NELEMENTS) = M.XNODES(M.NDIM*(FIX(j,k)-1)+1,1:M.NELEMENTS);
                        X2(j,k,1:M.NELEMENTS) = M.XNODES(M.NDIM*(FIX(j,k)-1)+2,1:M.NELEMENTS);
                        X3(j,k,1:M.NELEMENTS) = M.XNODES(M.NDIM*(FIX(j,k)-1)+3,1:M.NELEMENTS);
                    end
                    fill3(squeeze(X1(j,:,:)),squeeze(X2(j,:,:)),squeeze(X3(j,:,:)),'w'); hold on
                end
                plot3(M.NODES(:,2),M.NODES(:,3),M.NODES(:,4),'.','MarkerFaceColor','r',...
                    'MarkerEdgeColor','r');		hold off;
                hold off;			axis image;			title(txtit);
                xlabel('X[\mum]');	ylabel('Y[\mum]');	zlabel('z[\mum]');
                set(gcf,'position',[10 50 1900 950]);
            end
        end
    end
    %% The shape functions:
    methods(Static)
        function [N,NV,dNdXi,dNdXiV] = ShapeFunctions2D(NNODE,Xi,Eta)
            if NNODE == 4
                N  = [0.25*(1-Xi)*(1-Eta)    0  0.25*(1+Xi)*(1-Eta) 0  0.25*(1+Xi)*(1+Eta) 0  0.25*(1-Xi)*(1+Eta) 0
                    0  0.25*(1-Xi)*(1-Eta) 0  0.25*(1+Xi)*(1-Eta) 0  0.25*(1+Xi)*(1+Eta) 0  0.25*(1-Xi)*(1+Eta)];
                NV = [0.25*(1-Xi)*(1-Eta) 0.25*(1+Xi)*(1-Eta) 0.25*(1+Xi)*(1+Eta) 0.25*(1-Xi)*(1+Eta)];
                dNdXi   = [-0.25*(1-Eta)    0  0.25*(1-Eta)  0   0.25*(1+Eta)  0   -0.25*(1+Eta)   0
                    -0.25*(1-Xi)     0 -0.25*(1+Xi)   0   0.25*(1+Xi)   0    0.25*(1-Xi)    0
                    0 -0.25*(1-Eta) 0  0.25*(1-Eta)  0   0.25*(1+Eta)  0   -0.25*(1+Eta)
                    0 -0.25*(1-Xi)  0 -0.25*(1+Xi)   0   0.25*(1+Xi)   0    0.25*(1-Xi)];
                dNdXiV  = [-0.25*(1-Eta)  0.25*(1-Eta)  0.25*(1+Eta)  -0.25*(1+Eta)
                    -0.25*(1-Xi)  -0.25*(1+Xi)   0.25*(1+Xi)    0.25*(1-Xi)];
            elseif NNODE == 8
                N = [0.25*(1-Xi)*(1-Eta)*(-Xi-Eta-1)    0  0.25*(1+Xi)*(1-Eta)*(Xi-Eta-1) 	0  ...
                    0.25*(1+Xi)*(1+Eta)*(Xi+Eta-1) 	0  0.25*(1-Xi)*(1+Eta)*(-Xi+Eta-1) 	0  ...
                    0.5*(1-Xi^2)*(1-Eta) 				0  0.5*(1+Xi)*(1-Eta^2) 			0  ...
                    0.5*(1-Xi^2)*(1+Eta) 				0  0.5*(1-Xi)*(1-Eta^2) 			0; ...
                    0  0.25*(1-Xi)*(1-Eta)*(-Xi-Eta-1) 0  0.25*(1+Xi)*(1-Eta)*(Xi-Eta-1) 	0  ...
                    0.25*(1+Xi)*(1+Eta)*(Xi+Eta-1) 	0  0.25*(1-Xi)*(1+Eta)*(-Xi+Eta-1) 	0  ...
                    0.5*(1-Xi^2)*(1-Eta) 				0  0.5*(1+Xi)*(1-Eta^2) 			0  ...
                    0.5*(1-Xi^2)*(1+Eta) 				0  0.5*(1-Xi)*(1-Eta^2)];
                
                NV = [0.25*(1-Xi)*(1-Eta)*(-Xi-Eta-1) 	0.25*(1+Xi)*(1-Eta)*(Xi-Eta-1)  ...
                    0.25*(1+Xi)*(1+Eta)*(Xi+Eta-1) 	0.25*(1-Xi)*(1+Eta)*(-Xi+Eta-1) ...
                    0.5*(1-Xi^2)*(1-Eta) 				0.5*(1+Xi)*(1-Eta^2) 			...
                    0.5*(1-Xi^2)*(1+Eta) 				0.5*(1-Xi)*(1-Eta^2)];
                
                dNdXi = [-0.25*(1-Eta)*(-Xi-Eta-1)-0.25*(1-Xi)*(1-Eta)  0	0.25*(1-Eta)*(Xi-Eta-1)+0.25*(1+Xi)*(1-Eta)  0 	...
                    0.25*(1+Eta)*(Xi+Eta-1)+0.25*(1+Xi)*(1+Eta)  	0  -0.25*(1+Eta)*(-Xi+Eta-1)-0.25*(1-Xi)*(1+Eta) 0 	...
                    -Xi*(1-Eta)  									0   0.5*(1-Eta^2)  								 0 	...
                    -Xi*(1+Eta)  									0  -0.5*(1-Eta^2)  								 0; ...
                    -0.25*(1-Xi)*(-Xi-Eta-1)-0.25*(1-Xi)*(1-Eta) 	0  -0.25*(1+Xi)-0.25*(1+Xi)*(1-Eta) 			 0 	...
                    0.25*(1+Xi)*(Xi+Eta-1)+0.25*(1+Xi)*(1+Eta) 	0   0.25*(1-Xi)*(-Xi+Eta-1)+0.25*(1-Xi)*(1+Eta)  0 	...
                    -0.5*(1-Xi^2) 									0  -Eta*(1+Xi) 									 0 	...
                    0.5*(1-Xi^2) 									0  -Eta*(1-Xi) 									 0; ...
                    0 -0.25*(1-Eta)*(-Xi-Eta-1)-0.25*(1-Xi)*(1-Eta) 0  0.25*(1-Eta)*(Xi-Eta-1)+0.25*(1+Xi)*(1-Eta) 0 	...
                    0.25*(1+Eta)*(Xi+Eta-1)+0.25*(1+Xi)*(1+Eta)  	0  -0.25*(1+Eta)*(-Xi+Eta-1)-0.25*(1-Xi)*(1+Eta) 0 	...
                    -Xi*(1-Eta)  									0   0.5*(1-Eta^2)  								 0 	...
                    -Xi*(1+Eta)  									0  -0.5*(1-Eta^2);	...
                    0 -0.25*(1-Xi)*(-Xi-Eta-1)-0.25*(1-Xi)*(1-Eta) 0 -0.25*(1+Xi)-0.25*(1+Xi)*(1-Eta) 			 0 	...
                    0.25*(1+Xi)*(Xi+Eta-1)+0.25*(1+Xi)*(1+Eta) 	0   0.25*(1-Xi)*(-Xi+Eta-1)+0.25*(1-Xi)*(1+Eta)  0 	...
                    -0.5*(1-Xi^2) 									0  -Eta*(1+Xi) 									 0  ...
                    0.5*(1-Xi^2) 									0  -Eta*(1-Xi)];
                
                dNdXiV = [-0.25*(1-Eta)*(-Xi-Eta-1)-0.25*(1-Xi)*(1-Eta)  0.25*(1-Eta)*(Xi-Eta-1)+0.25*(1+Xi)*(1-Eta) 	...
                    0.25*(1+Eta)*(Xi+Eta-1)+0.25*(1+Xi)*(1+Eta)  -0.25*(1+Eta)*(-Xi+Eta-1)-0.25*(1-Xi)*(1+Eta) 	...
                    -Xi*(1-Eta)			0.5*(1-Eta^2) 			-Xi*(1+Eta)  		-0.5*(1-Eta^2); 			...
                    -0.25*(1-Xi)*(-Xi-Eta-1)-0.25*(1-Xi)*(1-Eta)  -0.25*(1+Xi)-0.25*(1+Xi)*(1-Eta) 				...
                    0.25*(1+Xi)*(Xi+Eta-1)+0.25*(1+Xi)*(1+Eta) 	 0.25*(1-Xi)*(-Xi+Eta-1)+0.25*(1-Xi)*(1+Eta) 	...
                    -0.5*(1-Xi^2) 		-Eta*(1+Xi) 			 0.5*(1-Xi^2) -Eta*(1-Xi)];
                
                
            end
        end
        % The 3D shape functions:
        function [N,NV,dNdXi,dNdXiV] = ShapeFunctions3D(NNODE,Xi,Eta,Zeta)
            if NNODE == 8
                N = [0.125*(1-Xi)*(1-Eta)*(1-Zeta)    0  0  0.125*(1+Xi)*(1-Eta)*(1-Zeta) 0  0  ...
                    0.125*(1+Xi)*(1+Eta)*(1-Zeta)    0  0  0.125*(1-Xi)*(1+Eta)*(1-Zeta) 0  0  ...
                    0.125*(1-Xi)*(1-Eta)*(1+Zeta) 	  0  0  0.125*(1+Xi)*(1-Eta)*(1+Zeta) 0  0 	...
                    0.125*(1+Xi)*(1+Eta)*(1+Zeta)    0  0  0.125*(1-Xi)*(1+Eta)*(1+Zeta) 0  0;	...
                    0  0.125*(1-Xi)*(1-Eta)*(1-Zeta) 0  0  0.125*(1+Xi)*(1-Eta)*(1-Zeta) 0  0  ...
                    0.125*(1+Xi)*(1+Eta)*(1-Zeta) 	  0  0  0.125*(1-Xi)*(1+Eta)*(1-Zeta) 0  0  ...
                    0.125*(1-Xi)*(1-Eta)*(1+Zeta) 	  0  0  0.125*(1+Xi)*(1-Eta)*(1+Zeta) 0  0  ...
                    0.125*(1+Xi)*(1+Eta)*(1+Zeta)    0  0  0.125*(1-Xi)*(1+Eta)*(1+Zeta) 0; 0	...
                    0  0.125*(1-Xi)*(1-Eta)*(1-Zeta) 0  0  0.125*(1+Xi)*(1-Eta)*(1-Zeta) 0  0  ...
                    0.125*(1+Xi)*(1+Eta)*(1-Zeta) 	  0  0  0.125*(1-Xi)*(1+Eta)*(1-Zeta) 0  0  ...
                    0.125*(1-Xi)*(1-Eta)*(1+Zeta)    0  0  0.125*(1+Xi)*(1-Eta)*(1+Zeta) 0  0  ...
                    0.125*(1+Xi)*(1+Eta)*(1+Zeta) 	  0  0  0.125*(1-Xi)*(1+Eta)*(1+Zeta)];
                
                NV = [0.125*(1-Xi)*(1-Eta)*(1-Zeta) 0.125*(1+Xi)*(1-Eta)*(1-Zeta) 0.125*(1+Xi)*(1+Eta)*(1-Zeta) ...
                    0.125*(1-Xi)*(1+Eta)*(1-Zeta) 0.125*(1-Xi)*(1-Eta)*(1+Zeta) 0.125*(1+Xi)*(1-Eta)*(1+Zeta) ...
                    0.125*(1+Xi)*(1+Eta)*(1+Zeta) 0.125*(1-Xi)*(1+Eta)*(1+Zeta)];
                
                dNdXi = [-0.125*(1-Eta)*(1-Zeta)  	0  0 	 0.125*(1-Eta)*(1-Zeta)	0  0 	0.125*(1+Eta)*(1-Zeta) 	0  0	...
                    -0.125*(1+Eta)*(1-Zeta)	0  0   	-0.125*(1-Eta)*(1+Zeta) 0  0    0.125*(1-Eta)*(1+Zeta)	0  0 	...
                    0.125*(1+Eta)*(1+Zeta)   	0  0    -0.125*(1+Eta)*(1+Zeta)	0  0;  -0.125*(1-Xi)*(1-Zeta)  	0  0	...
                    -0.125*(1+Xi)*(1-Zeta)  	0  0   	 0.125*(1+Xi)*(1-Zeta)  0  0    0.125*(1-Xi)*(1-Zeta)   0  0    ...
                    -0.125*(1-Xi)*(1+Zeta)  	0  0    -0.125*(1+Xi)*(1+Zeta)  0  0    0.125*(1+Xi)*(1+Zeta)  	0  0	...
                    0.125*(1-Xi)*(1+Zeta)  	0  0;	-0.125*(1-Xi)*(1-Eta)   0  0   -0.125*(1+Xi)*(1-Eta) 	0  0    ...
                    -0.125*(1+Xi)*(1+Eta) 		0  0    -0.125*(1-Xi)*(1+Eta) 	0  0    0.125*(1-Xi)*(1-Eta) 	0  0	...
                    0.125*(1+Xi)*(1-Eta) 		0  0     0.125*(1+Xi)*(1+Eta)   0  0    0.125*(1-Xi)*(1+Eta) 	0  0; 	...
                    0 -0.125*(1-Eta)*(1-Zeta)	0  0     0.125*(1-Eta)*(1-Zeta) 0  0    0.125*(1+Eta)*(1-Zeta)  0  0	...
                    -0.125*(1+Eta)*(1-Zeta) 	0  0    -0.125*(1-Eta)*(1+Zeta) 0  0    0.125*(1-Eta)*(1+Zeta) 	0  0    ...
                    0.125*(1+Eta)*(1+Zeta)   	0  0    -0.125*(1+Eta)*(1+Zeta) 0; 0   -0.125*(1-Xi)*(1-Zeta)  	0  0	...
                    -0.125*(1+Xi)*(1-Zeta)  	0  0     0.125*(1+Xi)*(1-Zeta) 	0  0    0.125*(1-Xi)*(1-Zeta)  	0  0    ...
                    -0.125*(1-Xi)*(1+Zeta)  	0  0    -0.125*(1+Xi)*(1+Zeta)  0  0    0.125*(1+Xi)*(1+Zeta)  	0  0  	...
                    0.125*(1-Xi)*(1+Zeta)  	0; 0  	-0.125*(1-Xi)*(1-Eta) 	0  0   -0.125*(1+Xi)*(1-Eta) 	0  0    ...
                    -0.125*(1+Xi)*(1+Eta) 		0  0    -0.125*(1-Xi)*(1+Eta) 	0  0    0.125*(1-Xi)*(1-Eta) 	0  0  	...
                    0.125*(1+Xi)*(1-Eta) 		0  0     0.125*(1+Xi)*(1+Eta) 	0  0    0.125*(1-Xi)*(1+Eta) 	0; 0 	...
                    0 -0.125*(1-Eta)*(1-Zeta)	0  0     0.125*(1-Eta)*(1-Zeta) 0  0    0.125*(1+Eta)*(1-Zeta)  0  0 	...
                    -0.125*(1+Eta)*(1-Zeta) 	0  0    -0.125*(1-Eta)*(1+Zeta) 0  0    0.125*(1-Eta)*(1+Zeta) 	0  0    ...
                    0.125*(1+Eta)*(1+Zeta)   	0  0    -0.125*(1+Eta)*(1+Zeta);0  0   -0.125*(1-Xi)*(1-Zeta)  	0  0	...
                    -0.125*(1+Xi)*(1-Zeta)  	0  0     0.125*(1+Xi)*(1-Zeta) 	0  0    0.125*(1-Xi)*(1-Zeta)  	0  0   	...
                    -0.125*(1-Xi)*(1+Zeta)  	0  0    -0.125*(1+Xi)*(1+Zeta)  0  0    0.125*(1+Xi)*(1+Zeta)  	0  0	...
                    0.125*(1-Xi)*(1+Zeta); 	0  0    -0.125*(1-Xi)*(1-Eta) 	0  0   -0.125*(1+Xi)*(1-Eta)   	0  0   	...
                    -0.125*(1+Xi)*(1+Eta) 		0  0    -0.125*(1-Xi)*(1+Eta)   0  0    0.125*(1-Xi)*(1-Eta) 	0  0	...
                    0.125*(1+Xi)*(1-Eta)   	0  0     0.125*(1+Xi)*(1+Eta) 	0  0    0.125*(1-Xi)*(1+Eta)];
                
                dNdXiV = [-0.125*(1-Eta)*(1-Zeta)   0.125*(1-Eta)*(1-Zeta)   0.125*(1+Eta)*(1-Zeta)   -0.125*(1+Eta)*(1-Zeta)   ...
                    -0.125*(1-Eta)*(1+Zeta)   0.125*(1-Eta)*(1+Zeta)   0.125*(1+Eta)*(1+Zeta)   -0.125*(1+Eta)*(1+Zeta);	...
                    -0.125*(1-Xi)*(1-Zeta)   -0.125*(1+Xi)*(1-Zeta)    0.125*(1+Xi)*(1-Zeta)     0.125*(1-Xi)*(1-Zeta)  	...
                    -0.125*(1-Xi)*(1+Zeta)   -0.125*(1+Xi)*(1+Zeta)    0.125*(1+Xi)*(1+Zeta)     0.125*(1-Xi)*(1+Zeta);	...
                    -0.125*(1-Xi)*(1-Eta)    -0.125*(1+Xi)*(1-Eta)    -0.125*(1+Xi)*(1+Eta)     -0.125*(1-Xi)*(1+Eta) 	...
                    0.125*(1-Xi)*(1-Eta)     0.125*(1+Xi)*(1-Eta)     0.125*(1+Xi)*(1+Eta)      0.125*(1-Xi)*(1+Eta)];
            end
        end
    end
end