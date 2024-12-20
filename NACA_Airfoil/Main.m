%% ASEN 3111 - Computational Assignment 03 - Main
    % Use Vortex panel method to obtain cL for sections of an
    % airfoil and compare the results to those from thin airfoil theory for 
    % symmetric and cambered airfoils. Use Prandtl Lifting-Line Theory to 
    % calculate span efficiency, and coefficients of lift and drag. Finally,
    % understand how the number of terms affects the error in the solution.
        %Author: Dev Patel
        %Collaborators: Jack Barney, Ginger Olson, Hideyuki Nakanishi
        %Date: 11/6/2022

%NOTE: For a NACA ABCD Airfoil, A = 100m, B = 10p, while CD = 100t
clear;
close all;
clc;
%% PROBLEM 1
disp('BEGINNING PROBLEM 1')
%Constants and Things
ALPHA = 10; %[Degrees]
VINF = 10;
N = 100;
c = 4;
[x,y] = NACA_Airfoils(0,0,.12,c,N);

%FUNCTION CALL FOR PROBLEM 1
c_L = Vortex_Panel(x,y,VINF,ALPHA);
c_Ltrue = 1.2; %From Experimental Data

for i = 100:150
    [xt,yt] = NACA_Airfoils(0,0,.12,c,i);
    c_Ltest = Vortex_Panel(xt,yt,VINF,10);
    if (abs(c_Ltest-c_Ltrue) < (1.2/100))
        k = 2*i;
        % PRINT K TO THE COMMAND WINDOW
        break
    end
end

disp('The total number of panels required to predict a sectional lift coefficient with 1% relative error:');disp(k)


%% PROBLEM 2
disp('BEGINNING PROBLEM 2')
%NOTE: For a NACA ABCD Airfoil, A = 100m, B = 10p, while CD = 100t
m = 0;
p = 0;
t1 = .06;
t2 = .12;
t3 = .24;
%Range of angles of attack to evaluate
alpha = linspace(-5,15,100);
% Lift coefficent vs angle of attack for all three airfoils
for i = 1:100
    [x06,y06] = NACA_Airfoils(m,p,t1,c,131);
    cL06(i) = Vortex_Panel(x06,y06,VINF,alpha(i));
    [x12,y12] = NACA_Airfoils(m,p,t2,c,131);
    cL12(i) = Vortex_Panel(x12,y12,VINF,alpha(i));
    [x24,y24] = NACA_Airfoils(m,p,t3,c,131);
    cL24(i) = Vortex_Panel(x24,y24,VINF,alpha(i));
    cLthin(i) = (2*pi/100)*alpha(i);
end
%Use polyfit function to extract the slope assuming a linear relationship
%between cL and alpha
P1 = polyfit(alpha,cL06,1);
a06 = P1(1);
P2 = polyfit(alpha,cL12,1);
a12 = P2(1);
P3 = polyfit(alpha,cL24,1);
a24 = P3(1);

slopes = [a06,a12,a24];
print = 'The sectional lift slopes for the NACA 0006, 0012,and 0024 Airfoils are respectively:';
% PRINT to COMMAND WINDOW, a = 2pi and alpha_0 = 0 deg
disp(print)
disp(slopes)
% AS predicted by thin airfoil theory:
disp('Thin airfoil theory predicts a slope of .11 [/degree] and a zero lift angle of attack of 0 degrees, consistent with our findings')



figure(1)
hold on
    title('Sectional Lift Coefficient vs. Angle of Attack for symmetric NACA Airfoils')
    xlabel('Angle of Attack [degrees]')
    ylabel('Coefficient of Lift')
    plot(alpha,cL06,'g')
    plot(alpha,cL12,'r')
    plot(alpha,cL24,'b')
    plot(alpha,cLthin,'y')
    legend('NACA 0006','NACA 0012','NACA 0024','slope of 2\pi')
    grid on
hold off

%% PROBLEM 3
disp('BEGINNING PROBLEM 3')
 %NOTE: For a NACA ABCD Airfoil, A = 100m, B = 10p, while CD = 100t
for i = 1:100
    [x00,y00] = NACA_Airfoils(0,0,.12,c,131);
    cL00(i) = Vortex_Panel(x00,y00,VINF,alpha(i));

    [x2,y2] = NACA_Airfoils(.02,.4,.12,c,131);
    cL2(i) = Vortex_Panel(x2,y2,VINF,alpha(i));

    [x44,y44] = NACA_Airfoils(.04,.4,.12,c,131);
    cL44(i) = Vortex_Panel(x44,y44,VINF,alpha(i));
end
P1 = polyfit(alpha,cL00,1);
a00 = P1(1);
P2 = polyfit(alpha,cL2,1);
a2 = P2(1);
P44 = polyfit(alpha,cL44,1);
a24 = P3(1);
% PRINT to COMMAND WINDOW, a = 2pi and alpha_0 = 0 deg

slopes = [a00,a2,a24];
print = 'The sectional lift slopes for the NACA 0012, 2412,and 4412 Airfoils are respectively:';
% PRINT to COMMAND WINDOW, a = 2pi and alpha_0 = 0 deg
disp(print)
disp(slopes)
% AS predicted by thin airfoil theory:
disp('Thin airfoil theory predicts a slope of .11 [/degree]')
% AS predicted by thin airfoil theory: also find alpha_0



figure(2)
hold on
    title('Sectional Lift Coefficent vs. Angle of Attack for cambered NACA Airfoils')
    xlabel('Angle of Attack [degrees]')
    ylabel('Coefficient of Lift')
    plot(alpha,cL00,'g')
    plot(alpha,cL2,'r')
    plot(alpha,cL44,'b')
    plot(alpha,cLthin,'y')
    legend('NACA 0012','NACA 2412','NACA 4412','slope of 2\pi')
    grid on
hold off

%% Problem 5
disp('BEGINNING PROBLEM 5')
%%Calculating Lift and Drag
b = 33.33;
c_r = 5.33;
c_t = 3 + (8.5/12);
a12 = a12*180/pi;
a24 = a24*180/pi;

%Comparison case
[etrue,c_Ltrue,c_Ditrue] = PLLT(b,a12,a24,c_t,c_r,0,deg2rad(-2),0,deg2rad(1),500);

rho = 17.56*10^-4;
V = 82*1.68781;
ss = b*(c_r+c_t)/2;
qinf = 0.5*rho*V^2;
L = c_Ltrue*qinf*ss;
D = c_Ditrue*qinf*ss;


%Error calculations
for i = 1:1000
    [et,c_Lt,c_Dit] = PLLT(33.33,a12,a24,3.7083,5.33,0,deg2rad(-2),0,deg2rad(1),i);
    if( abs(c_Lt-c_Ltrue)< (.001*c_Ltrue) && abs(c_Dit-c_Ditrue)< (.001*c_Ditrue))
        k1 = i;
        break
    end
end
for i = 1:500
    [et,c_Lt,c_Dit] = PLLT(33.33,a12,a24,3.7083,5.33,0,deg2rad(-2),0,deg2rad(1),i);
    if( abs(c_Lt-c_Ltrue)< (.01*c_Ltrue) && abs(c_Dit-c_Ditrue)< (.01*c_Ditrue))
         k2 = i;
        break
    end
end
for i = 1:500
       [et,c_Lt,c_Dit] = PLLT(33.33,a12,a24,3.7083,5.33,0,deg2rad(-2),0,deg2rad(1),i);
    if( abs(c_Lt-c_Ltrue)< (.1*c_Ltrue) && abs(c_Dit-c_Ditrue)< (.1*c_Ditrue))
         k3 = i;
        break
    end
end
disp('The total number of odd terms required to obtain lift and drag solutions with 10% relative error:');disp(k1)
disp('The total number of odd terms required to obtain lift and drag solutions with 1% relative error:');disp(k2)
disp('The total number of odd terms required to obtain lift and drag solutions with 0.1% relative error:');disp(k3)


