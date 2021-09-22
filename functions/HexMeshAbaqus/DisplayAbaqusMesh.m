function HexMeshfig = DisplayAbaqusMesh(vert_ObjPtCloud, Faces_ObjPtCloud)

% Invoke patch to display mesh
clear figure
figure(1)
HexMeshfig = patch('Faces', Faces_ObjPtCloud, 'Vertices', vert_ObjPtCloud, ...
    'EdgeColor', 'green', 'FaceColor', 'blue', 'FaceAlpha', .3, 'LineWidth', 1);
% view(3)
view([135, 15])
axis equal; axis image; axis xy
set(gcf,'position',[30 50 1300 950]);   
xlabel('X [Pixels]');                
ylabel('Y [Pixels]'); 
zlabel('Z [Pixels]');
view([45 45 15])
title('Hexahedron Mesh with 8 Nodes');

end
