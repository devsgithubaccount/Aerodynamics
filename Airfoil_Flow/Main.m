%% ASEN 3111 - Computational Assignment 02 - Main
    % Approximate Vortex strength at points along the chord. Use this value to
    % find the stream and velocity potential functions using the principles of 
    % superposition, as well as functions for the radius and angle of each point
    %from the the center of each vortex. Add the contributions from the uniform
    %flow. Find final Velocity components and solve for Cp and Pressure
    %Plot Findings.
        %Author: Dev Patel
        %Collaborators: Hideyuki Nakanishi
        %Date: 10/9/2022
clear
close all
clc
%% Inputs/Constants
c = 2; % [m]
alpha = deg2rad(12);%[rad]
V_inf = 68; % [m/s]
p_inf = 101.3*10^3; % [Pa]
rho_inf = 1.225; % [kg/m^3]
N = 100; % Number of Vortices
%% Calling Functions
Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N,1)
errorCalc(c,alpha,V_inf,p_inf,rho_inf,100)

%% Variations in Parameters
c_ = linspace(1,5,5);
alpha_ = deg2rad(linspace(5,20,5));
V_inf_ = linspace(50,100,5);
%% Plotting Variations in Chord Length: Figures 3-7
for i = 1:5
    Plot_Airfoil_Flow(c_(i),alpha,V_inf,p_inf,rho_inf,N,i+2)
end
%% Plotting Variations in Angle of Attack: Figures 8-12
for i = 1:5
    Plot_Airfoil_Flow(c,alpha_(i),V_inf,p_inf,rho_inf,N,i+7)
end
%% Plotting Variations in Freestream Velocity: Figures 13-17
for i = 1:5
    Plot_Airfoil_Flow(c,alpha,V_inf_(i),p_inf,rho_inf,N,i+12)
end





        
