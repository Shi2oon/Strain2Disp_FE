function [Uxout,Uyout,Xout,Yout] = UByNGP(Uxin,Uyin,Xin,Yin,NNODE)
    for j=1:size(Uxin,2)
        for k=1:size(Uxin,1)
            Uxout(k+(j-1)*NNODE) = Uxin(k,j);  
            Uyout(k+(j-1)*NNODE) = Uyin(k,j); 
            Xout1(k+(j-1)*NNODE) = Xin(k,j);
            Yout1(k+(j-1)*NNODE) = Yin(k,j);
        end
    end    
Xout  = linspace(min(min(Xout1)),max(max(Xout1)),length(unique(Xout1))*NNODE/2-NNODE/2);
Yout  = linspace(min(min(Yout1)),max(max(Yout1)),length(unique(Yout1))*NNODE/2-NNODE/2);
Uyout = reshape(Uyout,length(Yout),length(Xout));
Uxout = reshape(Uxout,length(Yout),length(Xout));
end