function [ret] = settings(request)
global data h dh kr kx dkx 
% Solver settings for parameter dependent problem in Manual
switch request
    case 'mesh'
        ret = 0:1/floor(data.LT/data.dL):1;
    case 'collMethod'
        ret = 'gauss';
        %ret = 'lobatto';
        %ret = 'uniform';
        %ret = 'user';
    case 'collPoints'
        ret = 3;
    case 'meshAdaptation'
        ret = 0;
    case 'errorEstimate'
        ret = 0;
    case 'absTolSolver'
        ret = 1e-7;
    case 'relTolSolver'
        ret = 1e-7;
    case 'absTolMeshAdaptation'
        ret = 1e-7;
    case 'relTolMeshAdaptation'
        ret = 1e-7;
    case 'minInitialMesh'
        ret = 50;        
    case 'finemesh'
        ret = 0;
    case 'allowTRM'
        ret = 0;
    case 'maxFunEvalsTRM'
        ret = 90000;
    case 'maxIterationsTRM'
        ret = 9000;
    case 'lambdaMin'
        ret = 0.001;
    case 'maxAdaptations'
        ret = 8;
    case 'switchToFFNFactor'
        ret = 0.5;
    case 'updateJacFactor'
        ret = 0.5;
    case 'K'
        ret = 200;  
    case 'thetaMax'
        ret=0.01;
    case 'maxCorrSteps'
        ret=5;
    case 'maxSteplengthGrowth'
        ret=5;
    case 'angleMin'
        ret=0.75;
    case 'meshFactorMax'
        ret=1.5;
    case 'PredLengthFactor'
        ret=1;
    case 'CorrLengthGrowth'
        ret=8;
end

end