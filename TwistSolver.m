function twist = TwistSolver(b,c,N,CL)
theta = [1:N]'/N*pi/2;
%c=c*ones(N,1);
twist=ones(N,1);
y=-cos(theta).*b;
n=1:2:2*N-1;
AR=(2*b)/mean(c);

An = zeros(N,1);
An(1)=CL/(AR*pi());

% calculates left side of the equation
lhs = (pi*c)/(4*b).*sin(theta);

rhs=0;
% Old solver that I derived using fourier series from the book + theory,
% works but there is a scaling issue that I decided to ignore for the sake
% of time since the equation from the PLLT powerpoint works fine
% for i=1:2:(N-1)
%    rhs = rhs+sin(theta*n).*((pi*c*n)/(4*b) + sin(theta))*An
% end
% rhs = (1/L)*rhs

% rhs equation from powerpoint
rhs=sin(theta*n).*((pi()*c*n)/(4*b) + repmat(sin(theta),1,N))*An;

% element-wise division of the rhs by the lhs to product angle of attack
a=rhs./lhs;

% outputs twist in radians
twist=a;
end