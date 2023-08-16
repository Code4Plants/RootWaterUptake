function [ret] = SystemRoot(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Adapted from original problem definition for parameter dependent problem in BVPSuite Manual
global data h dh kr kx dkx 
LT = data.LT; rc = data.rc; r= data.r; PC = data.PC; Psis = data.Psis;  dL = data.dL;   
G = -data.G; ML = data.ML; Mb = data.Mb;  
krr = @(t)2*pi*kr(t)*r; 
%
switch request
    case 'n'
        ret = 2;
    case 'orders'
        ret = [ 2 0 ];
    case 'problem'
               ret =  [(krr(t) *(z(1,1) - h(z(2,1)))-kx(t)*z(1,3)-dkx(t)*z(1,2)) ;                                 
                        (z(2,1)- Mb + G*(krr(t) *(z(1,1) - h(z(2,1)))))   
                ];      
   case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        ret(1, 1, 1) = krr(t);
        ret(1, 1, 2) = -dkx(t); 
        ret(1, 1, 3) = -kx(t);
        ret(1, 2, 1) =  (krr(t))*(-dh(z(2,1))); 
                         
        ret(2, 1, 1) = G*krr(t);
        ret(2, 2, 1) = 1-G*krr(t)*dh(z(2,1));
                
    if sum(isnan(ret))+sum(isinf(ret))>0
        disp('asd')
    end
    case 'interval'
        ret = [0, LT];
    case 'linear'
        ret = 0;
    case 'parameters'
        ret = 0;
    case 'c'
        ret = [];
    case 'BV'
        ret=[
            za(1,2);
            zb(1,1)-PC;             
            ];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        ret(1,1,1,2) = 1;
        ret(2,2,1,1) = 1;
    case 'dP'
        ret = [];               
    case 'dP_BV'
        ret =   [];  
    case 'initProfile'
        points = floor(LT/dL)+1; 
        ret.initialMesh    = linspace(0,LT,points);
        ret.initialValues  = [(Psis + 0.4*(PC-Psis))*exp(ret.initialMesh.^2*log(PC/(Psis + 0.4*(PC-Psis)))/LT^2);
                             linspace(Mb*0.99,Mb*0.99-(Mb-ML)*0.3,points)];             
    case 'EVP'
        ret = 0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate=0;
    otherwise
        ret = 0;                
end
end








