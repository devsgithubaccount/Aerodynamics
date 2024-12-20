function [x,y] = NACA_Airfoils(m,p,t,c,N)
    %OUTPUTS:
        % x: vector containing the x locations of the boundary points
        % y: vector containing the y locations of the boundary points
    %INPUTS:
        % m: maximum camber
        % p: location of max camber
        % t: thickness
        % c: chord length
        % N: number of employed panels
    % NOTE: points in x and y should be ordered around the surface in a
    % clockwise direction starting at the trailing edge. First and last
    % points should meet at the trailing edge
    x = linspace(0,c,N) /c;
    X = x*c;
    %PREALLOCATE VECTORS
    zero = zeros(1,N);
    y_t = zero;
    y_c = zero;
    x_U = zero;
    x_L = zero;
    y_U = zero;
    y_L = zero;
    xi = zero;

    for i = 1:N
        y_t(i) = (t/0.2)*c*( 0.2969*sqrt(x(i)) - 0.1260*x(i) - 0.3516*x(i)^2 + 0.2843*x(i)^3 - 0.1036*x(i)^4);
    % Calculate mean camber line, y_c, derivitive of camber line, and angle
    % the camber line makes with the horizontal.
        if (x(i) < p*c && x(i) > 0)
            y_c(i) = (m*X(i)/p^2) * (2*p -x(i));
            dy = 2*m/p*(1-X(i)/(p*c));
            xi(i) = atan(dy);
        else
            y_c(i) = (m*(c-X(i))/(1-p)^2) * (1+x(i)-2*p);
            dy = m/((1-p)^2)*(-2*X(i)/c+2*p);
            xi(i) = atan(dy);
        end
        x_U(i) = X(i) - y_t(i)*sin(xi(i));
        x_L(i) = X(i) + y_t(i)*sin(xi(i));
        y_U(i) = y_c(i) + y_t(i)*cos(xi(i));
        y_L(i) = y_c(i) - y_t(i)*cos(xi(i));
    end
    
%Put data in a form so the path is clockwise around the airfoil surface
    x_L = flip(x_L);
    y_L = flip(y_L);
    x = [x_L x_U(2:end)];
    y = [y_L y_U(2:end)];
end