function Inp4Grain(M4Nodes,M4Elements,resultsDir)
Nodes = M4Nodes(:,2:3);                                   
for i=1:size(M4Elements,1) 
    Elements{i}=M4Elements(i,2:end);                      
end

Elements_Sets{1}.Name='Set1';          
Elements_Sets{1}.Elements_Type='CPS4';
Elements_Sets{1}.Elements=1:length(M4Elements);

Filename=[resultsDir '\abaqus_model.inp'];
Matlab2Abaqus(Nodes,Elements,Elements_Sets,Filename)