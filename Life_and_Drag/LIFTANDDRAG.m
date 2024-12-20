function [A,N,L,D] = LIFTANDDRAG(Pu,Pl,x,y,alpha)
%Function for Lift and Drag
    %INPUTS:
        %Pressure along the upper surface
        %Pressure along lower surface
        %array of values x/c
        %Airfoil half thickness, y
        %Angle of Attack
    %OUTPUTS:
        %Axial Force
        %Normal Force
        %Lift per unit Span
        %Drag per unit span
    %Process:
    %1. Apply Trapazoidal rule along points x/c for Normal and Axial Force
    %creating "Panel"
    %2. Sum panels for Normal and Axial force then apply Lift and Drag Eqs
        N_u = zeros(length(x));
        A_u = zeros(length(x));
        N_l = zeros(length(x));
        A_l = zeros(length(x));

    for i = 1: length(x) -1
        del_x = abs(x(i+1)-x(i));
        del_y = y(i+1) - y(i);
        N_u(i) = -1.*(del_x) .* (Pu(i+1)+Pu(i))/2;
        A_u(i) = -1.*(del_y) .* (Pu(i+1)+Pu(i))/2;
        N_l(i) = (del_x) .* (Pl(i+1)+Pl(i))/2;
        A_l(i) = (del_y) .* (Pl(i+1)+Pl(i))/2;
    end

    N = sum(N_l,"all") + sum(N_u,"all");
    A = sum(A_l,"all") - sum(A_u,"all");
    L = N*cos(alpha) - A*sin(alpha);
    D = N*sin(alpha) + A*cos(alpha);

end