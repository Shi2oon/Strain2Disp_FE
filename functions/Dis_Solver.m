classdef Dis_Solver
   properties
       NDIM
       NNODE
       NGP
       NODES
       ELEMENTS
       XNODES
       IX
       K
       F
       NNODES
       NELEMENTS
       U
       E11
       E22
       E12
   end 
   
   methods
       %% Constructor:
       function S = Dis_Solver(NDIM,NNODE,NGP,NODES,ELEMENTS,XNODES,IX,K,F,E11,E22,E12)
         S.NDIM          = NDIM;
         S.NNODE         = NNODE;
         S.NGP           = NGP;
         S.NODES         = NODES;
         S.ELEMENTS      = ELEMENTS;
         S.XNODES        = XNODES;
         S.IX            = IX;
         S.K             = K;
         S.F             = F;
         S.E11           = E11; 
         S.E22           = E22;  
         S.E12           = E12;
       end
       
       % The number of elements in the mesh:
       function NELEMENTS = get.NELEMENTS(S)
           NELEMENTS = length(S.ELEMENTS(:,1)); 
       end
       
       % The nodal displacement:
  function U = get.U(S)
       
  KTK  = S.K'*S.K;
  try % Solve System of Linear Equations Using Pseudo-inverse.Note that One big
      % problem with pseudo-inverse; it’s a discontinuous mapping of the data  
      % when the matrix is not full rank. In other words, the pseudo-inverse of 
      %  a rank  deficient matrix is sensitive to noisy data. also relay
      % on SVD: Singular value decomposition (very slow)
      % Ref: Matrix Computation 4th edition section 5.5.5.
      fprintf('Solving using (first) Pseudo-inverse .. ');
      KKK  = pinv(full(KTK))*S.K'; 
      if isempty(KKK);         error('Solved matrix is empty');                 end
   catch err
      warning(err.message)
%       try  % ANOTHER TRY USING less expensive pseudoinverse
           % use Minimum norm least-squares solution to linear 
           % equation which used uses the COD (Complete Orthogonal
           % Decomposition) rather than SVD
           fprintf('failed\nnow using (2nd) lsqminnorm algorithm to solve .. ');
           KKK  = lsqminnorm(KTK,S.K'); 
           if isempty(KKK);         error('Solved matrix is empty');            end
              %{
       catch err
          warning(err.message)
          try % custom made pseudoinverse suitable for sparse matrices
              fprintf('failed\nnow using (3rd) pseudoinverse algorithm .. ')
              KKK  = pseudoinverse(KTK)*S.K';
              if isempty(KKK);         error('Solved matrix is empty');         end
          catch err
             warning(err.message)              
             try % more accurate and MATLAB standard and most efficient
                 % but the worse for deficient matrices 
                 fprintf('failed\nnow using (4th) MATLAB Backslash algorithm (not recommended)')
                 KKK  = KTK\S.K'; 
                 if isempty(KKK);         error('Solved matrix is empty');      end
             catch err
                warning(err.message)
                    % for non uniform matrix, similar oF backslah but check for 
                    % suitable conditioned method
                    fprintf('failed\nnow using (5th) MATLAB mldivide algorithm (suitable conditioned method)')
                    KKK  = mldivide(KTK,S.K');
                    if isempty(KKK);         error('Solved matrix is empty');	  end
             end
          end
       end
           %}
  end
  U  = KKK*S.F;   
  fprintf ('Done\n');
  end     
       
%% plot and separate displacement
       function [U1,U2,X1,X2,U3,X3] = Plot_U(S)
           % 2D Mesh
        if S.NDIM == 2
			% Multiplier:
            if S.NNODE == 4;					FIX = [1 2 3 4];       
            elseif S.NNODE == 8;				FIX = [1 5 2 6 3 7 4 8];             
            end
            % Create the coordinate matrix:        
            X1 = zeros(S.NNODE,S.NELEMENTS);	X2 = zeros(S.NNODE,S.NELEMENTS);
            U1 = zeros(S.NNODE,S.NELEMENTS);	U2 = zeros(S.NNODE,S.NELEMENTS);
            U_all  = S.U;
            for j=1:S.NNODE
				X1(j,1:S.NELEMENTS) = S.XNODES(S.NDIM*(FIX(j)-1)+1,1:S.NELEMENTS);
				X2(j,1:S.NELEMENTS) = S.XNODES(S.NDIM*(FIX(j)-1)+2,1:S.NELEMENTS);
				U1(j,1:S.NELEMENTS) = U_all(S.NDIM*(S.IX(FIX(j),1:S.NELEMENTS)-1)+1);
				U2(j,1:S.NELEMENTS) = U_all(S.NDIM*(S.IX(FIX(j),1:S.NELEMENTS)-1)+2);
            end
            %% plotting
			U3 = []; X3 = [];
			figure();
               s1=subplot(2,3,1); 	fill(X1,X2,S.E11(1:S.NNODE,:),'LineStyle','none');  brighten(0.5);
               s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; colormap jet;axis off;
               c = colorbar;        c.Label.String = '\epsilon_{11}';      caxis([-5e-3 5e-3]);	
               
               s1=subplot(2,3,2); 	fill(X1,X2,S.E12(1:S.NNODE,:),'LineStyle','none');  brighten(0.5);    
               s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; colormap jet;axis off;
               c = colorbar;        c.Label.String = '\epsilon_{12}';      caxis([-5e-3 5e-3]);	  
               
               s1=subplot(2,3,3); 	fill(X1,X2,S.E22(1:S.NNODE,:),'LineStyle','none');  brighten(0.5);    
               s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; colormap jet;axis off;
               c = colorbar;        c.Label.String = '\epsilon_{22}';      caxis([-5e-3 5e-3]);	
               
               s1=subplot(2,3,4);  	fill(X1,X2,U1,'LineStyle','none');     brighten(0.5);    
               s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; colormap jet;axis off;
               c = colorbar;        c.Label.String = ['u_{' num2str(1) '} [\mum]'];
               
               s1=subplot(2,3,6); 	fill(X1,X2,U2,'LineStyle','none');     brighten(0.5);    
               s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; colormap jet;axis off;
               c = colorbar;        c.Label.String = ['u_{' num2str(2) '} [\mum]'];
               
               s1=subplot(2,3,5);  	fill(X1,X2,'w');        		hold on
			   plot(S.NODES(:,2),S.NODES(:,3),'.','MarkerFaceColor','r',...
					'MarkerEdgeColor','r');		hold off;
			   set(gca,'Ydir','normal'); 
               s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; colormap jet;%axis off;
               xlabel('X[\mum]');  	ylabel('Y[\mum]');   
			   title('2D Mesh');    set(gca, 'YAxisLocation', 'right')
               set(gcf,'position',[10 50 1900 950]); 
			   
        elseif S.NDIM == 3
            % The facets:
            FIX = [1 2 3 4; 5 8 7 6; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
            % Create the coordinate matrix:        
            X1 = zeros(6,4,S.NELEMENTS);	X2 = zeros(6,4,S.NELEMENTS);	X3 = zeros(6,4,S.NELEMENTS);  
            U1 = zeros(6,4,S.NELEMENTS); 	U2 = zeros(6,4,S.NELEMENTS);    U3 = zeros(6,4,S.NELEMENTS);
            U_all  = S.U;
            % plot
            figure();                                  
            for j=1:6
                for k=1:4 
                    X1(j,k,1:S.NELEMENTS ) = S.XNODES(S.NDIM*(FIX(j,k)-1)+1,1:S.NELEMENTS );
                    X2(j,k,1:S.NELEMENTS ) = S.XNODES(S.NDIM*(FIX(j,k)-1)+2,1:S.NELEMENTS );
                    X3(j,k,1:S.NELEMENTS ) = S.XNODES(S.NDIM*(FIX(j,k)-1)+3,1:S.NELEMENTS );
                    U1(j,k,1:S.NELEMENTS ) = U_all(S.NDIM*(S.IX(FIX(j,k),1:S.NELEMENTS )-1)+1);
					U2(j,k,1:S.NELEMENTS ) = U_all(S.NDIM*(S.IX(FIX(j,k),1:S.NELEMENTS )-1)+2);
					U3(j,k,1:S.NELEMENTS ) = U_all(S.NDIM*(S.IX(FIX(j,k),1:S.NELEMENTS )-1)+3);
                end  
                subplot(1,4,1);	 fill3(squeeze(X1(j,:,:)),squeeze(X2(j,:,:)),...
                                    squeeze(X3(j,:,:)),'w');                   hold on
                subplot(1,4,2);  fill3(squeeze(X1(j,:,:)),squeeze(X2(j,:,:)),...
                    squeeze(X3(j,:,:)),squeeze(U1(j,:,:)),'LineStyle','none'); hold on
                subplot(1,4,3);  fill3(squeeze(X1(j,:,:)),squeeze(X2(j,:,:)),...
                    squeeze(X3(j,:,:)),squeeze(U2(j,:,:)),'LineStyle','none'); hold on
                subplot(1,4,4);  fill3(squeeze(X1(j,:,:)),squeeze(X2(j,:,:)),...
                    squeeze(X3(j,:,:)),squeeze(U3(j,:,:)),'LineStyle','none'); hold on
            end
            
			s1 = subplot(1,4,1);	
			plot3(S.NODES(:,2),S.NODES(:,3),S.NODES(:,4),'.r');	hold off;
			set(gca,'Ydir','normal');	s1.XDir='reverse';   	s1.YDir='reverse'; 
			axis image;axis xy; 		xlabel('Y[\mum]'); 		ylabel('X[\mum]');   
			zlabel('Z[\mum]');			title('3D Mesh');    	
			
			s2 = subplot(1,4,2);	brighten(0.5);          colormap jet;axis off;
			s2.XDir='reverse';      s2.YDir='reverse';      axis image;axis xy; 
            c  =colorbar;           cU(1,:) = c.Limits;     colorbar off; 
			title(['u_{' num2str(1) '}']);
			
			s3 = subplot(1,4,3);	brighten(0.5);      hold off; axis off;
			s3.XDir='reverse';  s3.YDir='reverse';      axis image;axis xy; 
            c  =colorbar;           cU(2,:) = c.Limits; colorbar off; 
			title(['u_{' num2str(2) '}']);
			
			s4 = subplot(1,4,4);	brighten(0.5);      hold off; axis off;
			s4.XDir='reverse';  s4.YDir='reverse';      axis image;axis xy;
            c  =colorbar;           cU(3,:) = c.Limits; colorbar off;
			title(['u_{' num2str(3) '}']);
            
			cbax  = axes('visible', 'off');             cU(abs(cU)==1)=0;
            caxis(cbax,[min(cU(:)) max(cU(:))]);
            h = colorbar(cbax, 'location', 'southoutside','position', [0.3513 0.0993 0.3 0.03] );
            h.Label.String = 'U [\mum]'; 
            set([s2 s3 s4],'clim',caxis);       set(gcf,'position',[1 41 1920 963]);  
		end
       end
   end 
end