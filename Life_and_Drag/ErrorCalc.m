function [n_L, n_D] = ErrorCalc(v_inf,rho,p,c,alpha)
%Function for Minimizing Error in Lift and Drag values
    %INPUTS:
        %Freestream Velocity
        %Freestream Density
        %Freestream pressure
        %Chord Length
        %Angle of Attack
    %OUTPUTS:
        %Number of integration points for a less than one percent relative error in Lift per unit Span
        %Number of integration points for a less than one percent relative error in Drag per unit span
    %Process:
    %1. Load Data
    %2. Create test case using large number of equispaced points along the
    %airfoil (5000).
    %3. Use lift and drag function in a for loop to increase integration
    %points with each iteration. (WILL HAVE TO RECALCULATE PRESSURE ALONG
    %UPPER AND LOWER SURFACES AND HALF THICKNESS FOR EACH ITERATION).
    %Once the values become within 1 percent of the test case, break the
    %for loop and return the number of integration points.
load cp.mat
    n_L = 0;
    n_D = 0;
    xx = 12;
    t = xx/100;
    x2 = linspace(0,c,5000); %

        Cp_u = fnval(Cp_upper,x2/c);
        Cp_l = fnval(Cp_lower,x2/c);

        Pu = (Cp_u .* rho .* v_inf^2 ./2) + p;
        Pl = (Cp_l .* rho .* v_inf^2 ./2) + p;

        y = (t*c/0.2) .* ( 0.2969.*((x2./c).^(1/2)) - 0.1260.*(x2./c) - 0.3516.*((x2./c).^2) + 0.2843.*((x2./c).^3) - 0.1036.*((x2./c).^4) );

    [A,N,L,D] = LIFTANDDRAG(Pu,Pl,x2,y,alpha);

    for i = 1:10000
        k = linspace(0,c,i);

        y = (t*c/0.2) .* ( 0.2969.*((k./c).^(1/2)) - 0.1260.*(k./c) - 0.3516.*((k./c).^2) + 0.2843.*((k./c).^3) - 0.1036.*((k./c).^4) );

        Cp_u = fnval(Cp_upper,k/c);
        Cp_l = fnval(Cp_lower,k/c);

        Pu = (Cp_u .* rho .* v_inf^2 ./2) + p;
        Pl = (Cp_l .* rho .* v_inf^2 ./2) + p;

        [A1,N1,L1,D1] = LIFTANDDRAG(Pu,Pl,k,y,alpha);
        if (abs(L-L1) < 0.01*L)
           n_L = i;
           break
        end
    end

    for i = 1:10000
        k = linspace(0,c,i);
        y = (t*c/0.2) .* ( 0.2969.*((k./c).^(1/2)) - 0.1260.*(k./c) - 0.3516.*((k./c).^2) + 0.2843.*((k./c).^3) - 0.1036.*((k./c).^4) );

        Cp_u = fnval(Cp_upper,k/c);
        Cp_l = fnval(Cp_lower,k/c);

        Pu = (Cp_u .* rho .* v_inf^2 ./2) + p;
        Pl = (Cp_l .* rho .* v_inf^2 ./2) + p;

        [A1,N1,L1,D1] = LIFTANDDRAG(Pu,Pl,k,y,alpha);
         if (abs(D-D1) < 0.01*D)
            n_D = i;
         break
         end
    end

end