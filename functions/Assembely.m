function [K,F] = Assembely(NDIM,NNODE,NGP,NODES,ELEMENTS,IX,dNdX,DUDX)
NDATA = length(DUDX(:,1,1));
NNODES = length(NODES(:,1));
NELEMENTS = length(ELEMENTS(:,1));
KELE = zeros(NDATA*NGP,NDIM*NNODE,NELEMENTS);
for j=1:NGP
    for k=1:NNODE
        if NDIM == 2
            KELE(NDATA*(j-1)+1,NDIM*(k-1)+1,1:NELEMENTS) ...
                = dNdX(1,k,j,1:NELEMENTS);
            KELE(NDATA*(j-1)+2,NDIM*(k-1)+2,1:NELEMENTS) ...
                = dNdX(2,k,j,1:NELEMENTS);
            
            KELE(NDATA*(j-1)+3,NDIM*(k-1)+1,1:NELEMENTS) ...
                = (1/2)*dNdX(2,k,j,1:NELEMENTS);
            KELE(NDATA*(j-1)+3,NDIM*(k-1)+2,1:NELEMENTS) ...
                = (1/2)*dNdX(1,k,j,1:NELEMENTS);
        elseif NDIM == 3
            KELE(NDATA*(j-1)+1,NDIM*(k-1)+1,1:NELEMENTS) ...
                = dNdX(1,k,j,1:NELEMENTS);
            KELE(NDATA*(j-1)+2,NDIM*(k-1)+2,1:NELEMENTS) ...
                = dNdX(2,k,j,1:NELEMENTS);
            KELE(NDATA*(j-1)+3,NDIM*(k-1)+3,1:NELEMENTS) ...
                = dNdX(3,k,j,1:NELEMENTS);
            
            KELE(NDATA*(j-1)+4,NDIM*(k-1)+1,1:NELEMENTS) ...
                = (1/2)*dNdX(2,k,j,1:NELEMENTS);
            KELE(NDATA*(j-1)+4,NDIM*(k-1)+2,1:NELEMENTS) ...
                = (1/2)*dNdX(1,k,j,1:NELEMENTS);
            
            KELE(NDATA*(j-1)+5,NDIM*(k-1)+1,1:NELEMENTS) ...
                = (1/2)*dNdX(3,k,j,1:NELEMENTS);
            KELE(NDATA*(j-1)+5,NDIM*(k-1)+3,1:NELEMENTS) ...
                = (1/2)*dNdX(1,k,j,1:NELEMENTS);
            
            KELE(NDATA*(j-1)+6,NDIM*(k-1)+2,1:NELEMENTS) ...
                = (1/2)*dNdX(3,k,j,1:NELEMENTS);
            KELE(NDATA*(j-1)+6,NDIM*(k-1)+3,1:NELEMENTS) ...
                = (1/2)*dNdX(2,k,j,1:NELEMENTS);
        end
    end
end

FELE = zeros(NDATA*NGP,NELEMENTS);
for j=1:NGP
    if NDIM == 2
        FELE(NDATA*(j-1)+1,1:NELEMENTS) = DUDX(1,j,1:NELEMENTS);
        FELE(NDATA*(j-1)+2,1:NELEMENTS) = DUDX(2,j,1:NELEMENTS);
        FELE(NDATA*(j-1)+3,1:NELEMENTS) = DUDX(3,j,1:NELEMENTS);
    elseif NDIM == 3
        FELE(NDATA*(j-1)+1,1:NELEMENTS) = DUDX(1,j,1:NELEMENTS);
        FELE(NDATA*(j-1)+2,1:NELEMENTS) = DUDX(2,j,1:NELEMENTS);
        FELE(NDATA*(j-1)+3,1:NELEMENTS) = DUDX(3,j,1:NELEMENTS);
        FELE(NDATA*(j-1)+4,1:NELEMENTS) = DUDX(4,j,1:NELEMENTS);
        FELE(NDATA*(j-1)+5,1:NELEMENTS) = DUDX(5,j,1:NELEMENTS);
        FELE(NDATA*(j-1)+6,1:NELEMENTS) = DUDX(6,j,1:NELEMENTS);
    end
end

% The assembely of the stiffness matrix and force vector:
K = sparse(NDATA*NGP*NELEMENTS,NDIM*NNODES);
F = sparse(NDATA*NGP*NELEMENTS,1);
for k=1:NGP
    for l=1:NDATA
        X = NDATA*NGP*([1:NELEMENTS]-1)+NDATA*(k-1)+l;
        x = NDATA*(k-1)+l;
        z = 1:NELEMENTS;
        for i = 1:NDIM
            Y = NDIM*(IX([1:NNODE],1:NELEMENTS)-1)+i;
            y = NDIM*([1:NNODE]-1)+i;
            XX = repmat(X,[NNODE,1]);
            K(sub2ind(size(K), XX(:),Y(:))) = squeeze(KELE(x,y,z));
        end
        F(X,1) = FELE(x,z);
    end
end
end