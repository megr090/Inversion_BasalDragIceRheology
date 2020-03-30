
function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

switch lower(UserVar.RunType)
    
    case 'icestream'
        
        
        x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
        alpha=0.;
        
        gamma=0.002;
        S=x*0 ;
        B0 = -500;
        B=B0+0*x ;
        b=B;
        s=270-gamma*x;
        
        
    case 'iceshelf'
        
        hmean=1000;
        b=zeros(MUA.Nnodes,1) ;
        S=zeros(MUA.Nnodes,1);
        B=S*0-1e10;
        s=hmean+b;
        
        alpha=0 ;
        
end


end

