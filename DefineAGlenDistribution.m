function  [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

n=3 ;
x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);
Lx = max(x)-min(x);
Ly=max(y)-min(y);
AGlen = AGlenVersusTemp(0);
AGlen=AGlen-0.5*cos(y*2*pi/Ly)*mean(AGlen); % define flow rate parameter distribution

if CtrlVar.AGlenisElementBased
    AGlen=AGlen+zeros(MUA.Nele,1);
else
    AGlen=AGlen+zeros(MUA.Nnodes,1);
end

if CtrlVar.doDiagnostic
    
    switch lower(UserVar.RunType)
        
        case 'iceshelf'
            
            
            x=MUA.coordinates(:,1) ;
            y=MUA.coordinates(:,2);
            
            if CtrlVar.AGlenisElementBased
                x=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
                y=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
            end
            
        case 'icestream'
            
    end
end
















end

