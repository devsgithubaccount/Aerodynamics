%%ASEN 3111 - Computational Assignment 1 - Main
%
% Author: Dev Patel
% Collaborators: Jack Barney, Hideyuki Nakanishi
% Date: 9/18/2022

%%Problem 1
clc
clear all
close all
%Function and Variable Definitions
theta = linspace(0,360, 360);
theta = deg2rad(theta);
C_p = @(angle) ((3/4)- 4.*(sin(angle)).^2 - 2.*sin(angle));
C_p_l = @(angle) -0.5.*((3/4)- 4.*(sin(angle)).^2 - 2.*sin(angle)).*sin(angle);
C_p_d = @(angle) -0.5.*((3/4)- 4.*(sin(angle)).^2 - 2.*sin(angle)).*cos(angle);
% Verification Case for Numerical Integration
c_l = integral(C_p_l,0,2*pi);
c_d = integral(C_p_d,0,2*pi);
%Trapazoidal Rule
trap_cL = trap(C_p_l,0,2*pi,100);
trap_cD = trap(C_p_d,0,2*pi,100);

simp_cL = simp(C_p_l, 0, 2*pi, 100);
simp_cD = simp(C_p_d, 0, 2*pi ,100);

x = 1:1:100;
for i = 1:100
    NcL_trap(i) = trap(C_p_l,0,2*pi,i);
    NcD_trap(i) = trap(C_p_d,0,2*pi,i);
    NcL_simp(i) = simp(C_p_l,0,2*pi,i);
    NcD_simp(i) = simp(C_p_d,0,2*pi,i);
end
figure(1)
hold on
plot(x,NcL_trap)
xlabel("Number of Panels, N")
ylabel("Sectional Lift Coefficient")
title("Sectional Lift Coefficient predicted by the Composite Trapazoidal Rule")
hold off

figure(2)
hold on
plot(x,NcD_trap)
xlabel("Number of Panels, N")
ylabel("Sectional Drag Coefficient")
title("Sectional Drag Coefficient predicted by the Composite Trapazoidal Rule")
hold off

figure(3)
hold on
plot(x,NcL_simp)
xlabel("Number of Panels, N")
ylabel("Sectional Lift Coefficient")
title("Sectional Lift Coefficient predicted by the Composite Simpson's Rule")
hold off

figure(4)
hold on
plot(x,NcD_simp)
xlabel("Number of Panels, N")
ylabel("Sectional Drag Coefficient")
title("Sectional Drag Coefficient predicted by the Composite Simpson' Rule")
hold off
%Error
trap_rule_panels = 0;
for i = 1:1000
    trap_cL_e = trap(C_p_l,0,2*pi,i);
    if(abs(c_l - trap_cL_e) < c_l / 100)
        trap_rule_panels = i;
        break
    end
end


simpsons_rule_panels = 0;
for i = 1:1000
    simp_cL_e = simp(C_p_l,0,2*pi,i);
    if(abs(c_l - simp_cL_e) < c_l / 100)
        simpsons_rule_panels = i;
        break
    end
end

%REFLECTION QUESTION 1: The trapazoidal rule, in this case, is able to
%reach a solution within one percent relative error with less panels than
%Simpson's rule. They would not perform the same way if given an array of
%pressure coefficient values these functions are made to use the
%trapazoidal/simpson's rule with function value at different points along
%the rotating cylinder. With an array of pressure values, the values would
%have to be matched with their respective postions along the cylinder, then
%a similar method can be used as the functions above.

%%Problem 2
%define variables

alpha = 15; %angle of attack[degrees]
alpha = deg2rad(alpha);
c = 2; % chord length [m]
v_inf = 30; %freestream airspreed [m/s]
rho = 1.225; %freestream Density [kg/m^3]
p = 101.3*10^3; %freestream Pressure [Pa]
R_e = 6000000; %Reynold's Number
load Cp.mat
n = 256;

x = linspace (0,c,n);
    Cp_u = fnval(Cp_upper,x/c);
    Cp_l = fnval(Cp_lower,x/c);

figure(5)
hold on 
plot(x,-Cp_u);
plot(x,-Cp_l);
xlabel("x/c")
ylabel("Coefficient of Pressure")
legend("Lower Surface","Upper Surface")
title("Coefficent of Pressure along the Upper and Lower Surfaces for a NACA0012 Airfoil")
hold off

Pu = (Cp_u .* rho .* v_inf^2 ./2) + p;
Pl = (Cp_l .* rho .* v_inf^2 ./2) + p;

xx = 12;
t = xx/100;
y = (t*c/0.2) .* ( 0.2969.*((x./c).^(1/2)) - 0.1260.*(x./c) - 0.3516.*((x./c).^2) + 0.2843.*((x./c).^3) - 0.1036.*((x./c).^4) );

[axial,normal,lift,drag] = LIFTANDDRAG(Pu,Pl,x,y,alpha);
[nL,nD] = ErrorCalc(v_inf,rho,p,c,alpha);
fprintf("Number of panels required to achieve a predicted sectional lift coefficient within 1 percent relative error:")
trap_rule_panels
simpsons_rule_panels
fprintf("Analytically Determined Lift and Drag coefficients:")
c_l
c_d

fprintf("Lift and Drag per Unit Span:")
lift
drag
fprintf("Integration points required for less than 1 percent relative error in Lift:");
nL
fprintf("Integration points required for less than 1 percent relative error in Drag:");
nD

%REFLECTIVE QUESTION 2: Based on my result for the minimum integration
%points required for less than one percent relative error, to
%experimentally calculate the lift and drag coefficients, one should
%measure the pressure along 195 points on each side of the airfoil. If the
%number of ports is limited, than more ports should be placed toward the
%leading edge of the airfoil, as the pressure gradient is larger toward the
%front of the airfoil. Also, this would help in getting a better
%approximation for half thickness, as the value varies increasingly less as
%x/c approaches 1. As a result, the lift and drag estimations would be more
%accurate since the placement of the ports would capture the pressure
%gradient more accurately in the data.