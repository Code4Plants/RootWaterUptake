function [zy,j] = RootStretch(data, h, M, kr, kx, dkx)
%
% Iterated matrix solution for coupled soil-root water transport of a
% simple root stretch with variable hydraulic properties
%
% Jan Graefe,  Leibniz Institute of Vegetable and Ornamental Crops 
% Grossbeeren 
% 2023
%
% Jan Graefe, Richard Pauwels, Michael Bitterlich
% Water flow within and towards plant roots – a new concurrent solution
% In silico Plants (submitted)
% 
% EMail: graefe@igzev.de
%-------------------------
% Inputs
%data:  structure with input information
%h, M, kr, kx, dkx: function handles for h(M), M(h), kr(z), kx(z), d kx(z)/dz  
% M = Matrix flux potential 
%------------------------
% Oututs
%zy  Array with columns z , Xylem potential Px(z), Root surface potential Ps(z)
%j   water flux rate at collar (L^3/T) 
% -----------------------
% transfer data structure to local variables  
PC = data.PC; % Xylem water potential at root collar, i.e. Dirichlet BC
Mb = data.Mb; % Bulk matrix potential 
LT = data.LT; % Root length (multiple of dL)
r  = data.r;  % Root radius
dL  = data.dL;%  % preferred root segment length 
%
% algorithm parameters
%

mmx   = 15.0; %number of iterations to estimate root surface water potential Ps(z)  
beta  = 1.25;  % empirical distortion factor of effective kx    
optionsFZ = optimset('Display','off', 'TolX', 1.0e-6);   % Options structure for fzero
%
%
ck    = [0.0,0.25,0.5,0.75,1]; % Sampling points for hydraulic properties 
sen   = 1.0/length(ck); 
Pss   =  h(Mb); 
%
g = data.G ; 
%
jmax = LT/dL; 
for i=1:jmax
    L(i) = LT-(i-1)*dL; 
    k1(i)  =   sum(sen*2*pi*kr(L(i)+ck*dL)*r);
    Ps(i)  =   Pss; 
    kxc(i) =   kx(L(i)); 
    kxh(i)    =   sum(sen*kx(L(i)+ ck*dL));  
    dxh(i)    =   sum(sen*dkx(L(i) +ck*dL)); 
    tau(i) = (k1(i)/kxh(i))^0.5;
    delta(i) = exp(-tau(i)*dL)-exp(tau(i)*dL); 
    D(i)     = exp(-tau(i)*dL)+exp(tau(i)*dL); 
    v(i)     = tau(i)/delta(i); 
end
kxc(jmax+1)= kx(0.0);  
%
% Built matrix A   in   A*Px = f(b, Ps)
%//////////////////////////////
A = zeros(jmax,jmax); 
%Top row of A 
j = 2; 
A(1,1) = -(kxc(j)*D(j)*v(j) +kxc(j)*D(j-1)*v(j-1)); 
A(1,2) = 2*kxc(j)*v(j); 
%middle row of a 
for j=3:jmax 
    A(j-1,j-2) = 2*kxc(j)*v(j-1); 
    A(j-1,j-1) = -(kxc(j)*D(j)*v(j) +kxc(j)*D(j-1)*v(j-1)); 
    A(j-1,j)   = 2*kxc(j)*v(j);
end
j = jmax+1; 
%bottom row of A
A(jmax,j-2) = 2*kxc(j)*v(j-1);   
A(j-1,j-1) = -kxc(j)*D(j-1)*v(j-1); 

%Built right side f(b, Ps)
%//////////////////////////////
b = zeros(jmax,1); 
%Top row of b, PC = collar WP 
j = 2; 
b(1) = kxc(j)*v(j)*Ps(j)*(2-D(j))+ kxc(j)*v(j-1)*Ps(j-1)*(2-D(j-1))-2*kxc(j)*v(j-1)*PC; 
%middle rows of b 
for j=3:jmax 
    b(j-1) = kxc(j)*v(j)*Ps(j)*(2-D(j))   + kxc(j)*v(j-1)*Ps(j-1)*(2-D(j-1));
end
j = jmax+1; 

%bottom row of b
b(jmax) = kxc(j)*v(j-1)*Ps(j-1)* (2-D(j-1)); 

% Solve for xylem potentials
Pxx = A\b;  
%Iterate for Psoil and variable kx(z)
for  m=1:mmx
    % Update A
    for i=1:jmax
        if i>1
           Jbo=Jb; 
           Jb = -(Pxx(i-1)-Pxx(i))/dL;           
           Jb2 = (Jbo-Jb)/dL; 
            if abs(Jb2)>0.0
                kxe =   kxh(i) + dxh(i)*Jb/Jb2 * beta ;     
            else
                kxe =   kxh(i); 
            end
            tau(i) = (k1(i)/kxe)^0.5;
            delta(i) = exp(-tau(i)*dL)-exp(tau(i)*dL); 
            D(i)     = exp(-tau(i)*dL)+exp(tau(i)*dL); 
            v(i)     = tau(i)/delta(i);             
        else
            Jb = -(PC-Pxx(i))/dL;  
        end
    end    
    j = 2; 
    A(1,1) = -(kxc(j)*D(j)*v(j) +kxc(j)*D(j-1)*v(j-1)); 
    A(1,2) = 2*kxc(j)*v(j); 
    %middle row of a 
    for j=3:jmax 
        A(j-1,j-2) = 2*kxc(j)*v(j-1); 
        A(j-1,j-1) = -(kxc(j)*D(j)*v(j) +kxc(j)*D(j-1)*v(j-1)); 
        A(j-1,j)   = 2*kxc(j)*v(j);
    end
    j = jmax+1; 
    %bottom row of A
    A(jmax,j-2) = 2*kxc(j)*v(j-1);   
    A(j-1,j-1) = -kxc(j)*D(j-1)*v(j-1); 
       
    % Update b
    
     j=1;
     LL = [M(Pxx(ceil(j))-0.001) Mb]; 
     f2 = @(x)x-Mb+g*k1(j)*(h(x)-0.5*(Pxx(j)+PC));  
     Mr = fzero(f2,LL, optionsFZ);  
     Ps(j) = h(Mr);
  for j =2:jmax  
     LL = [M(Pxx(ceil(j))-0.001) Mb];  
     f2 = @(x)x-Mb+g*k1(j)*(h(x)-0.5*(Pxx(j)+Pxx(j-1)));  
     Mr = fzero(f2,LL, optionsFZ);  
     Ps(j) = h(Mr);
  end     
    j = 2; 
    b(1) = kxc(j)*v(j)*Ps(j)*(2-D(j))+ kxc(j)*v(j-1)*Ps(j-1)*(2-D(j-1))-2*kxc(j)*v(j-1)*PC; 
    %middle rows of b 
    for j=3:jmax 
        b(j-1) = kxc(j)*v(j)*Ps(j)*(2-D(j))   + kxc(j)*v(j-1)*Ps(j-1)*(2-D(j-1));
    end
    j = jmax+1; 
    b(jmax) = kxc(j)*v(j-1)*Ps(j-1)* (2-D(j-1)); 
    Pxx = A\b;  
end
z = (LT:-dL:0); 
LL = [M(PC-0.001) Mb];
k10 = 2*pi*kr(LT)*r ; 
f2 = @(x)x-Mb+g*k10*(h(x)-PC);  
Mr = fzero(f2,LL, optionsFZ);  
Ps0 = h(Mr);

zy = flip([z', [PC Pxx']', [Ps0 Ps]']);
%Eq 19
Ps_i = Pxx(1)-Ps(1); Ps_j = PC-Ps(1); 
j = -kxc(1)*v(1)*( (exp(-tau(1)*dL)*Ps_i-Ps_j)*exp(tau(1)*dL) ... 
                     -(-1*exp(tau(1)*dL)*Ps_i+Ps_j)*exp(-tau(1)*dL)) ; 


