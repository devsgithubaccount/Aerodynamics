%% Diamond Airfoil Shock Profile
    % Theta-Beta-Mach Relation as seen in figure 9.9 in Anderson as well as
    % a model of lift and drag for a diamond airfoil in supersonic flow.
    % Finally, a brief analysis on the effects of angle of attack on the
    % diamond airfoil model.
        %Author: Dev Patel
        %Date: 12/4/2022
clear;clc;close all;

%% Problem 1: Theta-Beta-Mach Relation
%SET RANGES
n = 100;
M = [linspace(1,10,3*n/4) linspace(1.1,20,n/4 -1) 100000000] ;
theta = linspace(0,45.5,n);
gamma = 1.4;


for i = 1:n
    for j = 1:n
        Beta1(i,j) = ObliqueShockBeta(M(i),theta(j),gamma,'Weak');
        Beta2(i,j) = ObliqueShockBeta(M(i),theta(j),gamma,'Strong');
         if(abs(Beta2(i,j)) - abs(Beta1(i,j)) < 0.001)
             z = j;
           for k = z:n
           Beta2(i,k) = NaN;
           Beta1(i,k) = NaN;
           end
             break;
         end
    end
end

figure(1)
    hold on
    for i = 2:n
        plot(theta,real(Beta1(i,:)))
        plot(theta,real(Beta2(i,:))) 
        title('Oblique Shock Properties: Theta-Beta-Mach Relation')
        xlabel('Deflection Angle [deg]')
        ylabel('Shock wave Angle [deg]')
        axis([0 50 0 100]);
    end
hold off

%% Problem 2:Diamond Airfoil Function
[C_L,C_D] = DiamondAirfoil(2,5,10,5);
fprintf('The Sectional Lift and Wave Drag Coefficients for this Airfoil are: %d and %d, respectively',C_L,C_D)
%% Problem 3: Impact of Angle of Attack and Mach Number
M = [2,3,4,5];
n = 20;
e1 = 10;
e2 = 5;
alpha = linspace(-9.9,9.9,n);
cl = zeros(1,n);
cdw = zeros(1,n);
cl_lin = zeros(1,n);
cdw_lin = zeros(1,n);
    C1 = 1/(1+tand(e1)/tand(e2));
    C2 = 1/(1+tand(e2)/tand(e1));
    gl = C1*tand(e1)+C2*tand(e2);
    gu = gl;
for i = 1:length(M)
    for j = 1:n
    [cl(i,j), cdw(i,j)] = DiamondAirfoil(M(i),alpha(j),10,5);
    cl_lin(i,j) = (4*deg2rad(alpha(j))) / sqrt(M(i)^2 -1);
    cdw_lin(i,j) = 2*(2*deg2rad(alpha(j))^2+2*gl^2) / sqrt(M(i)^2-1);
    end
    figure (2)
    hold on
    plot(alpha,cl_lin(i,:),'--')
    plot(alpha,cl(i,:))
    legend('Linearized' ,'Predicted')
    title('Angle of Attack vs. Lift Coefficient')
    xlabel('Angle of Attack [deg]')
    ylabel('Lift Coefficient')
    hold off
    figure(3)
    hold on
    plot(alpha,cdw_lin(i,:),'--')
    plot(alpha,cdw(i,:))
    legend('Linearized','Predicted')
    title('Angle of Attack vs. Wave Drag Coefficient')
    xlabel('Angle of Attack [deg]')
    ylabel('Wave Drag Coefficient')
    hold off
end






