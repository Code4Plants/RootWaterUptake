% Driver for solution of a coupled soil-root water transport of a
% simple root stretch with variable hydraulic properties over length
%
% call to RootStretch.m  = Iterated matrix method of solution 
% call to bvpsuite2.0   =  DAE type reference solution using BVPSuite
% Install BVPsuite from  https://github.com/NumODEsTUW/bvpsuite2.0
%
% Jan Graefe, Richard Pauwels, Michael Bitterlich
% Water flow within and towards plant roots – a new concurrent solution
% In silico Plants (submitted)
%
clear all
close all
global data h dh kr kx dkx 
%
% Define and choose  soil parameters
% van Genuchten parameters from Alterra soil series
% %        alpha   Qr  Qs   n     lambda  Ks 
%
H(1,:)= [0.0144 0.02 0.46 1.534 -0.215  15.42];  %B3 sand
H(2,:)= [0.0195 0.01 0.59 1.109 -5.901  4.53];   %B11 clay
H(3,:)= [0.0099 0.01 0.43 1.288 -2.244  2.36];   %B8 loam
%
S       = 1;   %choose soil from above
%
data.LT = 30; % root length cm
RLD     = 1.0;  %root length density cm/cm^3  
data.rc = 1/(pi*RLD).^0.5; % equidistant parallel oriented and connected roots  
data.r  = 0.037;  % root radius,  cm
data.PC = -1.5*10^4;  % root collar potential  hPa
data.Psis    = -3000;   %  hPa   Soil bulk water potential 
data.dL = 1.0;     %  root discretisatin  
%
% Making tabular output for M(h) and h(M) functions (M = matrix flux
% potential) 
%
reten = []; pwp = 15000; 
Qr  = H(S,2); Qs = H(S,3); alpha = H(S,1); n = H(S,4); Ks = H(S,6); L = H(S,5);  m = 1-1/n; 
gam = @(h)(1./(1+(alpha*h).^n)).^m; 
K      = @(h)Ks*gam(h).^L.*(1-(1-gam(h).^(1/m)).^m).^2; 
Mflux  = @(h)integral(@(h1)K(h1),h,pwp,'ArrayValued',1) 
Theta  = @(h)Qr+(Qs-Qr)*gam(h);
for h = 10:10:pwp+50
      reten = [reten; [-h Mflux(h) Theta(h)]]  ;
end
% Define M(h) and h(M) based on tabular output
M = @(h)interp1(reten(:,1),reten(:,2),h,'pchip');
h = @(M)interp1(reten(:,2),reten(:,1),M,'pchip');  
%
data.Mb = M(data.Psis); 
%
% Piecewise linear defined hydraulic properties of roots 
% for exampel corn seminals from Meunier et al. 2018        
kxS   =   [1.76 1.91 5.45 14.85 15.85]*1E-4;  %cm^4 hPa-1 d-1      
kxLS  =   [7.0 10.0 21 40.0];  
krS   =   [0.47 0.47 0.43 0.28 0.28]*1E-4;  %cm hPa-1  d-1      
krLS  =   [19.0 20.0 30 100];  
kx    =   @(z)piecewk(z,kxS,kxLS,-1); 
dkx   =   @(z)piecewk(z,kxS,kxLS,1); 
kr    =   @(z)piecewk(z,krS,krLS,-1);

% Choose between two geometry functions
%steady rate, Schröder et al.  2008
%-----------------------------------------
rc = data.rc; r=data.r; 
a=0.53; p = rc/r; 
data.G = -((a^2*rc^2)/(2*r*(p^2 - 1))-(r*(p^2*log(p) + p^2*log(a) + 1/2))/(p^2 - 1))/(2*pi*r);
%
%steady state, Graefe et al. 2019 
%----------------------------------------
%data.G    =  0.5*(0.5*(r^2 + 0.53^2*rc^2)+(r^2+rc^2)*log(0.53*rc/r))/(pi*(rc^2-r^2));  
%
%--------------------------------------------------
% run simulation with approximate matrix solution
%--------------------------------------------------
[zy_A,j] = RootStretch(data, h, M, kr, kx, dkx)
flux_A =j; 

%--------------------------------------------------
% run simulation with reference solution 
%--------------------------------------------------
% Install the BVPsuite solver from  https://github.com/NumODEsTUW/bvpsuite2.0
% and add path for exact solution
addpath('E:\bvpsuite2.0\bvpsuite2.0');
%addpath('E:\bvpsuite2.0\Examples');
% Solve the DAE problem  for a corn seminal root with a length 
% of LT cm in an environment of 3 different soils at defined root length density RLD and root radius r
% Define addtionally required dh/dM function
ML = M(-14500);
data.ML = ML; 
dM0 = 0.01;
dM = max(dM0,dM0*(data.Mb-ML)); 
dh = @(M)(h(M+dM)-h(M))/dM; 
[z,y,s] = bvpsuite2('SystemRoot', 'MsettingsRoot');
%SystemRoot.m contains the system equations, Jacobi matrices and
%initializations
%MsettingsRoot.m is for general solver settings
Ps = zeros(length(z),1); 
 for i= 1:length(z)
   Ps(i) = h(y(2,i)); 
 end
zy_E = [z'  y(1,:)' Ps];
dx = (s.x1tau(end)-s.x1tau(end-1)); 
flux_E = -(s.valx1tau(1,end)-s.valx1tau(1,end-1))/dx*kx(data.LT); 
% Plot results
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot1 = plot(zy_A(:,1),zy_A(:,2),'b--'); hold on;  
plot2 = plot(zy_A(:,1),zy_A(:,3),'r--'); 
set(plot1,'DisplayName','Xylem-Matrix','Color',[0 0 1]);
set(plot2,'DisplayName','Root surface-Matrix','Color',[1 0 0]);
plot3 = plot(zy_E(:,1),zy_E(:,2),'b-'); 
plot4 = plot(zy_E(:,1),zy_E(:,3),'r-'); 
ylabel({'Water potential'});
xlabel({'Length from root tip'});
set(plot3,'DisplayName','Xylem-DAE','Color',[0 0 1]);
set(plot4,'DisplayName','Root surface-DAE','Color',[1 0 0]);
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.26649958624776 0.257402106706825 0.31578946808227 0.19914039568095]);

disp(['Flux at Collar (Matrix):',num2str(flux_A,4), ' L^3/T'])
disp(['Flux at Collar (DAE):',num2str(flux_E,4),  ' L^3/T'])


