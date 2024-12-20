function errorCalc(c,alpha,V_inf,p_inf,rho_inf,n)
%% Comparison Case
    [Vtrue,Ptrue] = Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,500);
    %% Maximum Pressure and Velocity Values
    Vmax = max(max(Vtrue));
    Pmax = max(max(Ptrue));

    errV = zeros(1,n);
    errP = zeros(1,n);
    %% Calculate Error for varying number of vortices
    for i = 1:n
         [V,P] = Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,i);
         errV(i) = abs(Vmax - max(max(V))) / Vmax;
         errP(i) = abs(Pmax - max(max(P))) / Pmax;
    end
    %% PLOTS
    figure(2)
    subplot(1,2,1)
    hold on
    plot(linspace(1,n,n),errV)
    title('Error in Velocity vs Number of Panels, N')
    xlabel('Number of Panels')
    ylabel('Error in Velocity')
    hold off
    subplot(1,2,2)
    hold on
    plot(linspace(1,n,n),errP)
    title('Error in Pressure vs Number of Panels, N')
    xlabel('Number of Panels')
    ylabel('Error in Pressure')
    hold off
end