
function [UserVar,s,b,S,B,alpha,h]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)
x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

load('Bedmap2GriddedInterpolantModifiedBathymetry.mat');
load('REMA_SurfaceElevationGriddedInterpolant_200m.mat')

s = FGel(x,y);
b = Fb(x,y);
S = zeros(length(s),1);
B = FB(x,y);
alpha = 0;
h = s - b;
end


