function [X,Y,Ux,Uy] = xgLOBAL(mesh,U1,U2,elX,elY)
UX = U1(:);			UY = U2(:);
 mesh.Data=NaN(size( mesh.Data));
% used_nodes(1,:) = unique(Elements(:));
% used_nodes(2,:) = 1:length(used_nodes);
% for ik=1:length(U1) %ELEMENTS WITH 4
%     elUx(Elements(:,2:5)==used_nodes(1,ik)) = UX(ik);
%     elUy(Elements(:,2:5)==used_nodes(1,ik)) = UY(ik);
% end

for i=1:length(mesh.xy)
    pos = find(elX(:)==mesh.xy(1,i) & elY(:)==mesh.xy(2,i));
    if(pos)
        mesh.Data(1,i) = UX(pos(1));
        mesh.Data(2,i) = UY(pos(1));
    end
end

% find non zero displacment nodes
logic = ((sqrt(mesh.Data(1,:).^2+mesh.Data(2,:).^2))~=0);

% global U
x = reshape(logic.*mesh.xy(1,:), mesh.winodow(2),mesh.winodow(1));
y = reshape(logic.*mesh.xy(2,:),mesh.winodow(2),mesh.winodow(1));
X = x(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));
Y = y(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));

% global d
ux = reshape(logic.*mesh.Data(1,:),mesh.winodow(2),mesh.winodow(1));
uy = reshape(logic.*mesh.Data(2,:), mesh.winodow(2),mesh.winodow(1));
Ux = ux(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));
Uy = uy(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));