classdef DataUm_v2
    properties
        NDIM
        NNODE
        NGP
        DATATYPE
        NODES
        ELEMENTS
        IX
        dNdX
        NNODES
        NELEMENTS
        DUDX
        NDATA
        KELE
        FELE
    end
    
    methods
        % Constructor:
        function D = DataUm_v2(NDIM,NNODE,NGP,DATATYPE,NODES,ELEMENTS,IX,dNdX,DUDX)
            D.NDIM          = NDIM;
            D.NNODE         = NNODE;
            D.NGP           = NGP;
            D.DATATYPE      = DATATYPE;
            D.NODES         = NODES;
            D.ELEMENTS      = ELEMENTS;
            D.IX            = IX;
            D.dNdX          = dNdX;
            D.DUDX          = DUDX;
        end
        
        % The number of gradients in the data:
        function NDATA = get.NDATA(D)
            NDATA = length(D.DUDX(:,1,1));
        end
        
        % The number of nodes in the mesh:
        function NNODES = get.NNODES(M)
            NNODES = length(M.NODES(:,1));
        end
        
        % The number of elements in the mesh:
        function NELEMENTS = get.NELEMENTS(M)
            NELEMENTS = length(M.ELEMENTS(:,1));
        end
        
        % The element stiffness matrix at the Gauss points:
        function KELE = get.KELE(D)
            DNELEMENTS = D.NELEMENTS;           DNNODE = D.NNODE;
            DNDATA     = D.NDATA;               DNGP   = D.NGP;
            KELE = zeros(DNDATA*DNGP,D.NDIM*DNNODE,DNELEMENTS);
            for j=1:DNGP
                for k=1:DNNODE
                    if D.NDIM == 2
                        KELE(D.NDATA*(j-1)+1,D.NDIM*(k-1)+1,1:D.NELEMENTS) ...
                            = D.dNdX(1,k,j,1:D.NELEMENTS);
                        KELE(D.NDATA*(j-1)+2,D.NDIM*(k-1)+2,1:D.NELEMENTS) ...
                            = D.dNdX(2,k,j,1:D.NELEMENTS);
                        
                        KELE(D.NDATA*(j-1)+3,D.NDIM*(k-1)+1,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(2,k,j,1:D.NELEMENTS);
                        KELE(D.NDATA*(j-1)+3,D.NDIM*(k-1)+2,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(1,k,j,1:D.NELEMENTS);
                    elseif D.NDIM == 3
                        KELE(D.NDATA*(j-1)+1,D.NDIM*(k-1)+1,1:D.NELEMENTS) ...
                            = D.dNdX(1,k,j,1:D.NELEMENTS);
                        KELE(D.NDATA*(j-1)+2,D.NDIM*(k-1)+2,1:D.NELEMENTS) ...
                            = D.dNdX(2,k,j,1:D.NELEMENTS);
                        KELE(D.NDATA*(j-1)+3,D.NDIM*(k-1)+3,1:D.NELEMENTS) ...
                            = D.dNdX(3,k,j,1:D.NELEMENTS);
                        
                        KELE(D.NDATA*(j-1)+4,D.NDIM*(k-1)+1,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(2,k,j,1:D.NELEMENTS);
                        KELE(D.NDATA*(j-1)+4,D.NDIM*(k-1)+2,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(1,k,j,1:D.NELEMENTS);
                        
                        KELE(D.NDATA*(j-1)+5,D.NDIM*(k-1)+1,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(3,k,j,1:D.NELEMENTS);
                        KELE(D.NDATA*(j-1)+5,D.NDIM*(k-1)+3,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(1,k,j,1:D.NELEMENTS);
                        
                        KELE(D.NDATA*(j-1)+6,D.NDIM*(k-1)+2,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(3,k,j,1:D.NELEMENTS);
                        KELE(D.NDATA*(j-1)+6,D.NDIM*(k-1)+3,1:D.NELEMENTS) ...
                            = (1/2)*D.dNdX(2,k,j,1:D.NELEMENTS);
                    end
                end
            end
        end
        
        % The element force at the Gauss points:
        function FELE = get.FELE(D)
            FELE = zeros(D.NDATA*D.NGP,D.NELEMENTS);
            for j=1:D.NGP
                if D.NDIM == 2
                    FELE(D.NDATA*(j-1)+1,1:D.NELEMENTS) = D.DUDX(1,j,1:D.NELEMENTS);
                    FELE(D.NDATA*(j-1)+2,1:D.NELEMENTS) = D.DUDX(2,j,1:D.NELEMENTS);
                    FELE(D.NDATA*(j-1)+3,1:D.NELEMENTS) = D.DUDX(3,j,1:D.NELEMENTS);
                elseif D.NDIM == 3
                    FELE(D.NDATA*(j-1)+1,1:D.NELEMENTS) = D.DUDX(1,j,1:D.NELEMENTS);
                    FELE(D.NDATA*(j-1)+2,1:D.NELEMENTS) = D.DUDX(2,j,1:D.NELEMENTS);
                    FELE(D.NDATA*(j-1)+3,1:D.NELEMENTS) = D.DUDX(3,j,1:D.NELEMENTS);
                    FELE(D.NDATA*(j-1)+4,1:D.NELEMENTS) = D.DUDX(4,j,1:D.NELEMENTS);
                    FELE(D.NDATA*(j-1)+5,1:D.NELEMENTS) = D.DUDX(5,j,1:D.NELEMENTS);
                    FELE(D.NDATA*(j-1)+6,1:D.NELEMENTS) = D.DUDX(6,j,1:D.NELEMENTS);
                end
            end
        end
        
        % The assembely of the stiffness matrix and force vector:
        function [K,F] = Assembely(D)
            K = sparse(D.NDATA*D.NGP*D.NELEMENTS,D.NDIM*D.NNODES);
            F = sparse(D.NDATA*D.NGP*D.NELEMENTS,1);
            DKELE      = D.KELE;                DFELE  = D.FELE;
            DNELEMENTS = D.NELEMENTS;           DNNODE = D.NNODE;
            DNDATA     = D.NDATA;               DNGP   = D.NGP;
            DNDIM      = D.NDIM;                DIX    = D.IX;
            for k=1:DNGP
                for l=1:DNDATA
                    X = DNDATA*DNGP*([1:DNELEMENTS]-1)+DNDATA*(k-1)+l;
                    x = DNDATA*(k-1)+l;
                    z = 1:DNELEMENTS;
                    for i = 1:DNDIM
                        Y = DNDIM*(DIX([1:DNNODE],1:DNELEMENTS)-1)+i;
                        y = DNDIM*([1:DNNODE]-1)+i;
                        XX = repmat(X,[DNNODE,1]);
                        K(sub2ind(size(K), XX(:),Y(:))) = squeeze(DKELE(x,y,z));
                    end
                    F(X,1) = DFELE(x,z);
                end
            end
        end
    end
    
end