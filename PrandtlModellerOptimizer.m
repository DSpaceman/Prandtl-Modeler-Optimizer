rho=0.002007956; %slug/ft^3
Vinf=70; %ft/s
b = 4.89; % semi-span
%c = linspace(5,3,5).';   % constant chord
c = 1.110429448
CL(1) = 0.5;
CL(2) = 1.0;
N = 5;
span = linspace(-b,b,2*N)';


twist(:,1) = TwistSolver(b,c,N,CL(1));
twist(:,2) = TwistSolver(b,c,N,CL(2));

[gammaplot(:,1),dplot(:,1),lplot(:,1),wplot(:,1),aiplot(:,1),CDi(1),Di(1),CLout(1),L(1),Elliptical(1)] = LiftDistribution(twist(:,1),b,c,N,Vinf,rho);
[gammaplot(:,2),dplot(:,2),lplot(:,2),wplot(:,2),aiplot(:,2),CDi(2),Di(2),CLout(2),L(2),Elliptical(2)] = LiftDistribution(twist(:,2),b,c,N,Vinf,rho);

hold on
clear fig

    subplot(2,3,1)
    plot(span,gammaplot,'-+')
    title('Circulation (gamma)')
    legend({'CL=0.5','CL=1.0'}) % Legend is the same for every plot

    subplot(2,3,2)
    plot(span,wplot,'-+')
    title('Downwash')
    
    subplot(2,3,3)
    plot(span,aiplot*180/pi,'-+')
    title('Induced AoA')
    
    subplot(2,3,4)
    plot(span,lplot,'-+')
    title('Lift Distribution')
    
    subplot(2,3,5)
    plot(span,dplot,'-+')
    title('Drag Distribution')
    
    subplot(2,3,6)
    plot(linspace(0,7.5,N)',flipud(twist*180/pi),'-+')
    title('Wing Twist')
    
hold off


dispCenterLine05 = ['Centerline AoA @ CL = 0.5: ', num2str(twist(N,1)*180/pi),' degrees'];
dispCenterLine10 = ['Centerline AoA @ CL = 1.0: ', num2str(twist(N,2)*180/pi),' degrees'];
dispTipAoA05 = ['Tip AoA @ CL = 0.5: ', num2str(twist(1,1)*180/pi),' degrees'];
dispTipAoA10 = ['Tip AoA @ CL = 1.0: ', num2str(twist(1,2)*180/pi),' degrees'];

dispCDi05 = ['CDi @ CL = 0.5: ', num2str(CDi(1))];
dispCDi10 = ['CDi @ CL = 1.0: ', num2str(CDi(2))];

dispDi05 = ['Di @ CL = 0.5: ', num2str(Di(1))];
dispDi10 = ['Di @ CL = 1.0: ', num2str(Di(2))];

dispL05 = ['L @ CL = 0.5: ', num2str(L(1))];
dispL10 = ['L @ CL = 1.0: ', num2str(L(2))];

dispEl05 = ['Elliptical @ CL = 0.5? ', num2str(Elliptical(1))];
dispEl10 = ['Elliptical @ CL = 1.0? ', num2str(Elliptical(2))];


disp(dispCenterLine05)
disp(dispCenterLine10)
disp(dispTipAoA05)
disp(dispTipAoA10)
disp(dispCDi05)
disp(dispCDi10)
disp(dispDi05)
disp(dispDi10)
disp(dispL05)
disp(dispL10)
disp(dispEl05)
disp(dispEl10)
disp(newline)
disp(['The wing twist optimal for the wing at CL=1.0 appears to have twice the twist at every element',newline, 'compared to the wing with a twist optimal at CL=0.5. Because of this, one would need to design',newline,'the aircraft to fly at a specific CL. The optimal twist distribution for a wing meant to fly',newline,'at one CL will not be the same for a wing meant to fly at another. What this means for manufacturing',newline,'is that in addition to the increase difficulty of manufacturing a variably twisted wing, it would only',newline,'be effective for specific flight conditions.'])