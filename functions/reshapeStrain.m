function [va] = reshapeStrain (DATAout)
    va.X   = reshape(DATAout(:,1),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
    va.Y   = reshape(DATAout(:,2),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
    va.E11 = reshape(DATAout(:,3),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
    va.E22 = reshape(DATAout(:,4),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
    va.E12 = reshape(DATAout(:,5),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
    
end