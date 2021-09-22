function Matlab2Abaqus(M4Nodes,M4Elements,resultsDir)
Nodes = M4Nodes(:,2:3);                                   
for i=1:size(M4Elements,1) 
    Elements{i}=M4Elements(i,2:end);                      
end

Elements_Sets{1}.Name='Set1';          
Elements_Sets{1}.Elements_Type='CPS4';
Elements_Sets{1}.Elements=1:length(M4Elements);

Filename=[resultsDir '\abaqus_model.inp'];
Bysction(Nodes,Elements,Elements_Sets,Filename)
end




function Bysction(Nodes,Elements,Elements_Sets,Filename)
fileID = fopen(Filename, 'w');
fprintf(fileID,'*NODE, NSET=NODE\n');   %Generate Nodes in Input File
[NNode, ND] = size(Nodes);
if ND == 2  %2D                               
    for i=1:1:NNode   
        fprintf(fileID,[num2str(i) ', ' num2str(Nodes(i,1)) ', '...
            num2str(Nodes(i,2)) '\n']); 
    end
elseif ND==3  %3D 
    for i=1:1:NNode
    fprintf(fileID,[num2str(i) ', ' num2str(Nodes(i,1)) ', ' ...
        num2str(Nodes(i,2)) ', ' num2str(Nodes(i,3)) '\n']); 
    end      
end
fprintf(fileID,'\n');

%% Generate Elements in Input File
for i=1:length(Elements_Sets)
    fprintf(fileID,strcat('*ELEMENT, ELSET=',Elements_Sets{i}.Name,...
        ', TYPE=',Elements_Sets{i}.Elements_Type,'\n'));
    for j=1:length(Elements_Sets{i}.Elements) %Loop for the elements in the elements set
        IE  = Elements_Sets{i}.Elements(j);        %Elements indices in elements sets
        NNN = [num2str(IE) ', '];
        for k=1:length(Elements{IE})    
            NNN = [NNN num2str(Elements{IE}(k)) ', '];
        end
    NNN = NNN(1:end-2);
    fprintf(fileID,[NNN '\n']);
    end
fprintf(fileID,'\n');
end
fprintf(fileID,'\n');

%%
fclose(fileID);
end