function [V,P] = Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N)
%PLOT AIRFOIL STREAMLINES
%%CONSTANTS & LIMITS
    q_inf = 0.5*rho_inf*V_inf^2;
    del_x = c/N;
    xmin = -c;
    xmax = 2*c;
    ymin = -2*c;
    ymax = 2*c;
    y_ = zeros(N,1);
    x_c = linspace(del_x/2,c - del_x/2,N);
%%Define Mesh Grid of points in the area of interest
    [x,y] = meshgrid(linspace(xmin,xmax,N),linspace(ymin,ymax,N));
%%Functions for vortex strength,radius/angle from vortex center
    strength = @(V_inf,alpha,x_c) 2*alpha*V_inf*sqrt((1-x_c)/(x_c));
    r = @(x,y,x1,y1) sqrt((x-x1).^2 + (y).^2);
    theta = @(x,y,xc,yc) atan2(y,(x-xc));
    Vx = 0;
    Vy = 0;
    %% Solve Velocity in the X,Y directions
    for i = 1:N
        Vx = Vx + (sin(theta(x,y,x_c(i),y_(i))).*strength(V_inf,alpha,x_c(i)./c) *del_x./(2*pi*r(x,y,x_c(i),y_(i))));
        Vy = Vy - (cos(theta(x,y,x_c(i),y_(i))).*strength(V_inf,alpha,x_c(i)./c) *del_x./(2*pi*r(x,y,x_c(i),y_(i))));
    end
    %%Total X and Y velocities using principles of superposition
    Vx = Vx + V_inf*cos(alpha);
    Vy = Vy + V_inf*sin(alpha);
    %%Solve for Pressure Distribution using Velocity and Pressure
    %%Coefficient
    V = sqrt(Vx.^2 + Vy.^2);
    cP = 1 - (V./V_inf).^2;
    P = cP.*q_inf + p_inf;
