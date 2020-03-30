function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)



x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2);
Lx=max(x)-min(x); Ly=max(y)-min(y);

if CtrlVar.CisElementBased
    xC=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
    yC=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
else
    xC=x; yC=y;
end


if CtrlVar.AGlenisElementBased
    xA=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
    yA=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
else
    xA=x ; yA=y;
end
%

%% Define boundary conditions of adjoint problem
% Generally there is nothing that needs to be done here If BCsAdjoint is not
% modified, then Ua will define the BCs of the adjoint problem based on the BCs
% of the forward problem. 
%BCsAdjoint=BCs; % periodic BCs of forward model -> periodic BCs of adjoint model



%%  Covariance matrices of priors
% 
if CtrlVar.AGlenisElementBased
    CAGlen=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
else
    CAGlen=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
end

if strcmpi(CtrlVar.Inverse.Regularize.Field,'cov')
    Err=1e-2 ; Sigma=1e3 ; DistanceCutoff=10*Sigma;
    
    if CtrlVar.CisElementBased
        [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
    else
        [CC]=SparseCovarianceDistanceMatrix(xC,yC,Err,Sigma,DistanceCutoff);
    end
    
else
    if CtrlVar.CisElementBased
        CC=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
    else
        CC=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
    end
end

Priors.CovAGlen=CAGlen;
Priors.CovC=CC;

Priors.s=F.s;
Priors.b=F.b;
Priors.S=F.S;
Priors.B=F.B;


Priors.m = 3;
Priors.n = 3;

Priors.C = 0.01 + zeros(MUA.Nnodes,1);
Priors.AGlen = AGlenVersusTemp(-10)+zeros(MUA.Nnodes,1);

Priors.rho=F.rho;
Priors.rhow=F.rhow;

%% Define start values
InvStartValues.C = zeros(MUA.Nnodes,1)+0.01;    
InvStartValues.AGlen = AGlenVersusTemp(-10)+zeros(MUA.Nnodes,1);
InvStartValues.m=Priors.m;
InvStartValues.n=Priors.n;
%% Define measurements and measurement errors

fprintf(' Creating synthetic data for iC \n')

CtrlVar.doDiagnostic=1;

Priors.TrueC=[];
Priors.TrueAGlen=[];


load('L8-2015-GriddedInterpolants-1000m.mat');
x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
Meas.us = FGu(x,y);
Meas.vs = FGv(x,y);
Meas.ws = FGu(x,y)*0;


usError = 30;
vsError = 30;
wsError = 1;

Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,usError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,vsError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.wsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,wsError.^2,MUA.Nnodes,MUA.Nnodes);
 
% if add errors

Meas.us=Meas.us+UserVar.AddDataErrors*usError.*randn(MUA.Nnodes,1);
Meas.vs=Meas.vs+UserVar.AddDataErrors*vsError.*randn(MUA.Nnodes,1);
Meas.ws=Meas.ws+UserVar.AddDataErrors*wsError.*randn(MUA.Nnodes,1);

end
