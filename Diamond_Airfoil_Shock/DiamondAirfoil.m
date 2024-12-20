%% Problem 2: Lift and Drag coefficients for Diamond-Wedge Airfoil
function [c_l,c_dw] = DiamondAirfoil(M,alpha,epsilon1,epsilon2)
%INPUTS:
    %Mach: M
    %Angle of Attack: Alpha [deg]
    %Leading Edge half-angle: epsilon1 [deg]
    %Trailing Edge half-angle: epsilon2 [deg]
%OUTPUTS:
    %Sectional Lift coefficient: c_l
    %Sectional Wave Drag coefficient: c_dw
%Process: use the normal shock, isentropic flow, and Prandtl-Meyer
%properties to solve for changes in the flow across regions of the airfoil.
%Use the pressure ratios and geometry of the airfoil to solve for
%coefficients.

%simplifications
    e1 = epsilon1;
    e2 = epsilon2;
    gamma = 1.4;
%Solve for oblique shock angles at LE:
    shock_deflection_upper = abs(e1 - alpha); %Deflection vs freestream
    shock_deflection_lower = abs(e1 + alpha); 
    %Conditions where this approach does not apply
    if (shock_deflection_upper > 45.5)
        disp('This Configuration will lead to a bow shock')
        c_l = NaN;
        c_dw = Nan;
    end
    if (shock_deflection_upper > 45.5)
        disp('This Configuration will lead to a bow shock')
        c_l = NaN;
        c_dw = Nan;
    end
    fan_deflection = e1 + e2;

    [Beta_upper] = ObliqueShockBeta(M,shock_deflection_upper,gamma,'Weak');
    [Beta_lower] = ObliqueShockBeta(M,shock_deflection_lower,gamma,'Weak');
   %Normal shock Relations
    M_n1_upper = M*sind(Beta_upper);
    M_n1_lower = M*sind(Beta_lower);

    [~, ~, P_12, ~, Mn2, ~, P_pitot_1_2] = flownormalshock([gamma gamma],[M_n1_upper M_n1_lower], 'mach');
   
    P21_upper = P_12(1);
    P21_lower = P_12(2);

    M2 = [0 0];
    M2(1) = Mn2(1)/ (sind(Beta_upper - shock_deflection_upper));
    M2(2) = Mn2(2)/ (sind(Beta_lower - shock_deflection_lower));
%Expansion Fans: Solve for New Prandtl Meyer angles on upper and lower surfaces
    [~, v, ~] = flowprandtlmeyer([gamma gamma],M2,'mach');
    v_new_upper = v(1) + fan_deflection;
    v_new_lower = v(2) + fan_deflection;

    [M3_u, ~ ,~] = flowprandtlmeyer(gamma,v_new_upper,'nu');
    [M3_l, ~, ~] = flowprandtlmeyer(gamma,v_new_lower,'nu'); %Mach Number after expansion
    M3 = [M3_u M3_l];

    [~,~,P32,~,~] = flowisentropic([gamma gamma],M3); %NOTE P32 is local static over stagnation, isentropic so P_o is the same as region 2
    
    %Pressure Ratios (Static upstream over downstream)
    P21 = P21_upper;
    P41 = P21_lower;
   
    %Pressure Ratios (Local Static over Stagnation)
    P303 = P32(1);
    P505 = P32(2);
    
    %Rayleigh Pitot Pressure Ratios(static upstream over total downstream)
    P021 = 1 / P_pitot_1_2(1);
    P041 = 1 / P_pitot_1_2(2);
   
    %Solve for Pressure Ratios with P1 as denominator; NOTE*ISENTROPIC
    %EXPANSION TOTAL PRESSURE IS THE SAME ACROSS THE TOP SURFACE AFTER
    %SHOCK
    P31 = P303*P021;
    P51 = P505*P041;
    
    %Find Chord lengths in terms of Epsilon 1&2
    C1 = 1/(1+tand(e1)/tand(e2));
    C2 = 1/(1+tand(e2)/tand(e1));

    V = tand(e1)*(C1);
    % SOLVE FOR CL AND CD
    Ca = (P21*V - P31*V + P41*V - P51*V) / (0.5*gamma*M^2);
    Cn = (-P21*C1 - P31*C2 + P41*C1 + P51*C2) / (0.5*gamma*M^2);
    
    c_l = Cn*cosd(alpha)-Ca*sind(alpha);
    c_dw = Cn*sind(alpha)+Ca*cosd(alpha);

end