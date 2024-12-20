%% PROBLEM 4
function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    %OUTPUTS:
        % e: Span Efficiency Factor
        % c_L: Lift Coefficient
        % c_Di: Induced Drag Coefficient
    %INPUTS:
        % a0_t: cross sectional lift slope at the tips [/radian]
        % a0_r: cross sectional lift slope at the root [/radian]
        % c_t: chord at the tip [ft]
        % c_r: chord at the root [ft]
        % aero_t: zero lift angle of attack at the tip [degrees]
        % aero_r: zero lift angle of attack at the root [degrees]
        % geo_t: geometric angle of attack at the tip [degrees]
        % geo_r: geometric angle of attack at the root [degrees]
        % N: number off odd terms to include in series expansion
    %This function solves for the span efficiency factor ,coefficient of
    %lift and induced drag using Prandtl Lifting-Line Theory. The function
    %uses N odd terms of the approximating series to approximate vortex
    %strength distribution along the span. N test locations are evaluated
    %from 0 to pi/2.
%define values across span
%Theta, c, alpha, and slope values across span, and PLL Eq in Matrix Form
for i = 1:N
    theta(i) = pi*i/(2*N);
end
    c = c_r + (c_t-c_r)*cos(theta);
    alpha = (geo_r + (geo_t-geo_r)*cos(theta));
    alpha_zero = (aero_r + (aero_t-aero_r)*cos(theta));
    a0 = a0_r + (a0_t-a0_r)*cos(theta);

for i = 1:N
    for j = 1:N
        K(i,j) = ((4*b/(a0(i)*c(i)))*sin((2*j-1)*theta(i)))+( (2*j-1)*sin((2*j-1)*theta(i))/sin(theta(i)));
    end
    F(i) = alpha(i)-alpha_zero(i);
end
%Solve for Fourier Coefficients
A = K\F';
% for i = 1:N/2
%     A(i*2) = 0;
% end
delta = zeros(1,N);
sum = 0;
%Solve for Delta
for j = 2:N
    delta(j) = (2*j-1)*(A(j)/A(1))^2;
    sum = sum + delta(j);
end

%Solve for CL e and CDi
S = b*(c_r+c_t)/2;
AR = b^2/S;
c_L = A(1)*pi*AR;
e = 1/(1+sum);
c_Di = c_L^2/(pi*e*AR);
end
