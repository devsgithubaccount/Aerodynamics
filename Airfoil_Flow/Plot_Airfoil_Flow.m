%% PLOT AIRFOIL STREAMLINES
function Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N,b)%%PLOT AIRFOIL STREAMLINES
%% CONSTANTS & LIMITS
    q_inf = 0.5*rho_inf*V_inf^2;
    del_x = c/N;
    xmin = -c;
    xmax = 2*c;
    ymin = -2*c;
    ymax = 2*c;
    xoverc = linspace(0,c,N);
    y_ = zeros(N,1);
    x_c = linspace(del_x/2,c - del_x/2,N);
%% Define Mesh Grid of points in the area of interest
    [x,y] = meshgrid(linspace(xmin,xmax,N),linspace(ymin,ymax,N));
%% Functions for vortex strength,radius/angle from vortex center
    strength = @(V_inf,alpha,x_c) 2*alpha*V_inf*sqrt((1-x_c)/(x_c));
    r = @(x,y,x1,y1) sqrt((x-x1).^2 + (y).^2);
    theta = @(x,y,xc,yc) atan2(y,(x-xc));
%% Uniform stream and velocity potential throughout grid
    u_stream = V_inf*cos(alpha)*y-V_inf*sin(alpha)*x;
    u_potential = V_inf*cos(alpha)*x+V_inf*sin(alpha)*y;
%% Initialize values for vortex stream, potential, as well as velocity.
    v_stream = 0;
    v_potential = 0;
    Vx = 0;
    Vy = 0;
    %% Solve for Vortex Stream/Potential and Velocity in the X,Y directions
    for i = 1:N
        v_stream = v_stream + strength(V_inf,alpha,x_c(i) ./c) * del_x / (2*pi) * log(r(x,y,x_c(i),y_(i)));
        v_potential = v_potential + -1 * strength(V_inf,alpha,x_c(i) ./c) * del_x / (2*pi) .*mod(theta(x,y,x_c(i),x_c(i)*sin(alpha)),2*pi);
        Vx = Vx + (sin(theta(x,y,x_c(i),y_(i))).*strength(V_inf,alpha,x_c(i)./c) *del_x./(2*pi*r(x,y,x_c(i),y_(i))));
        Vy = Vy - (cos(theta(x,y,x_c(i),y_(i))).*strength(V_inf,alpha,x_c(i)./c) *del_x./(2*pi*r(x,y,x_c(i),y_(i))));
    end
    %% Total X and Y velocities using principles of superposition
    Vx = Vx + V_inf*cos(alpha);
    Vy = Vy + V_inf*sin(alpha);
    %% Solve for Pressure Distribution using Velocity and Pressure
    %% Coefficient
    V = sqrt(Vx.^2 + Vy.^2);
    cP = 1 - (V./V_inf).^2;
    P = cP.*q_inf + p_inf;
    %% Total Stream,Potential using Principles of Superposition
    stream = (u_stream + v_stream);
    potential = (v_potential + u_potential);

    %% Define levels to be passed into the contourf function

    slevmin = stream(1,N);
    slevmax = stream(N,N/2);
    pmin = min(min(potential));
    pmax = max(max(potential));
    plevels= linspace(pmin,pmax,50);
    lmin = min(min(P));
    lmax = max(max(P));
    level = linspace(lmin,lmax,60)';
    slevels = linspace(slevmin,slevmax,50)';
%%Plots for Streamlines, Equipotential Lines, and Pressure Distributions
%%passing in stream, velocity and pressure matrices into the counterf
%%function.
        figure (b)
        subplot(1,3,1)
    hold on
    contourf(x,y,stream,slevels)
    plot(xoverc,y_,'r','linewidth',1);
    axis equal
    axis([xmin xmax ymin ymax])
    title('Streamlines')
    colorbar;
    xlabel('x')
    ylabel('y')
    hold off
        subplot(1,3,2)
     hold on
     contourf(x,y,(potential),plevels)
     plot(xoverc,y_,'r','linewidth',1);
     axis equal
     axis([xmin xmax ymin ymax])
     title('Equipotential Lines')
     colorbar;
     xlabel('x')
     ylabel('y')
     hold off
         subplot(1,3,3)
      hold on 
      contourf(x,y,P,level)
      plot(xoverc,y_,'r','linewidth',5);
      axis equal
      axis([xmin xmax ymin ymax])
      title('Pressure Distribution along x/c')
      colorbar;
      xlabel('x')
      ylabel('y')
      hold off
end     