%% this funciton re-arrange maps line elememnt by data
function [DataOut] = ByNGP(DataIn,NGP)
%% getting in Number of gussian
round(length(DataIn)/4)
    for j=1:round(length(DataIn)/4)
        for k=1:NGP
            DataOut(k,j) = DataIn(k+(j-1)*NGP);  
        end
    end    


end