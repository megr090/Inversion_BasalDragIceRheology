function [UserVar,C,m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)


m=3 ;
            sx = 10e3/4;
            sy = 10e3/4;
            x=MUA.coordinates(:,1) ;
            y=MUA.coordinates(:,2);
xC=x;
Lx=max(x)-min(x);
if CtrlVar.CisElementBased
    C=zeros(MUA.Nele,1)+1/20^m;
else
    C = 1;
    C=C.*(1+40*exp(-(x.^2/sx^2+y.^2/sy^2))+100*exp(-((x-7e4).^2/sx^2+(y-7e3).^2/sy^2))+5*exp(-((x+7e4).^2/sx^2+(y-7e3).^2/sy^2))+10*exp(-((x-7e4).^2/sx^2+(y+7e3).^2/sy^2))+5*exp(-((x+4e4).^2/sx^2+(y+7e3).^2/sy^2))+100*exp(-((x-2e4).^2/sx^2+(y-8e3).^2/sy^2))+20*exp(-((x+7e4).^2/sx^2+(y-3e3).^2/sy^2))+10*exp(-((x-6e4).^2/sx^2+(y+7e3).^2/sy^2))+70*exp(-((x+5e4).^2/sx^2+(y+5e3).^2/sy^2))+20*exp(-((x+3e4).^2/sx^2+(y-3e3).^2/sy^2))+10*exp(-((x-2e4).^2/sx^2+(y+7e3).^2/sy^2))+70*exp(-((x+1e4).^2/sx^2+(y+9e3).^2/sy^2)));
end

if CtrlVar.doDiagnostic
    
    switch lower(UserVar.RunType)
        
        case 'icestream'
            
            
            x=MUA.coordinates(:,1) ;
            y=MUA.coordinates(:,2);
            
            if CtrlVar.CisElementBased
                x=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
                y=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
            end
            
            sx=10e3/2 ; sy=10e3/2;
            
        case 'iceshelf'
            
    end
end





end
