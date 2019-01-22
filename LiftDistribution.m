function [gammaplot,dplot,lplot,wplot,aiplot,CDi,Di,CL,L,Elliptical] = LiftDistribution(a,b,c,N,Vinf,rho)

theta = [1:N]'/N*pi/2;
y=-cos(theta)*b;
n=1:2:2*N-1;
rhs=0;
q=(1/2)*rho*(Vinf^2)*(2*b*c);

% RHS equation from powerpoint
rhs=sin(theta*n).*((pi()*c*n)/(4*b) + repmat(sin(theta),1,N));

% LHS equation from powerpoint
lhs = (pi*c)/(4*b).*(a).*sin(theta);

% Reverse matrix division to solve for fourier coefficients
A = rhs\lhs;

% gamma calculation for semi-span
gamma=4*b*Vinf*sin(theta*n)*A;

% mirror concatnation of semi-span distribution
gammaplot = [gamma;flipud(gamma)];

% Calculation of aspect ration
AR=(2*b)/mean(c);
CL=AR*pi()*A(1)

% Fourier series implementation of downwash, seems to have issues 
% (is exactly one third the downwash from the equation in use).
% The downwash equation being used seems to allow the integrative CDi
% solvers to work.
% The integrative CL solver works regardless of which downwash equation I
% use because the gamma distribution is fine.
%w=-(sin(theta*n)*(A.*n'))./sin(theta)
w= (gamma-pi*Vinf*a.*c)./(pi.*c)


wplot = [w;flipud(w)]

%If 1, the imported semi-span chord/twist distribution is elliptical 
Elliptical = sum((gamma/(4*Vinf*A(1)*b)).^2 + (y/b).^2)/N;

%Induced angle of attack
ai=atan(-w./Vinf);
aiplot = [ai;flipud(ai)];


%Integrative CL finding method, to confirm fourier method
%CLInt = 2/(Vinf*2*b*c)*2*trapz(y,gamma) %0.4994 @ input twist for CL 0.5


%Integrative method, produces a CDi 0.0001 lower than ideal likely due to
%trapz error
%CDiInt1 = 2/((Vinf)*2*b*c)*2*trapz(y,gamma.*sin(ai))
%CdiInt2 = 2/(Vinf^2*(2*b*c))*2*trapz(-y,w.*gamma)

%method from prandtl.m, produces a CDi that is 0.005 higher, 
%probably because of rectangular integration.
%CDi2=2/(Vinf*c*2*b)*2*sum(gamma.*sin(ai)*0.5)

%ideal result based off fourier series, close enough to integrative methods
CDi=CL^2/pi/AR*(1+n(2:end)*(A(2:end).^2/A(1).^2));

%Integrative total lift equation
%L = rho*Vinf*2*trapz(y,gamma)

%Lift & Drag equations based off of fourier series derived CL
L = q*CL;
Di = q*CDi; 

%Actual lift & drag distribution
lplot = rho*Vinf*gammaplot;
dplot = -rho.*wplot.*gammaplot;

%plot(linspace(0,b,2*N),gammaplot)
