function z_eff = TwoSPP_Weighted_Div_Func(lambda ,w_0, l_int, SPP_phases, dist_aft_SPP1,dist_aft_SPP2,theta_1, theta_2, measured_divergence)

% Theoretical Weighted Diveregence (tested against confirmed weighted diveregence and data)

% Only built for two SPPs currently but z_eff can be expanded using the
% developed theory.

% nominal, only need two random lengths to calc, change if wanted
sim_lengths = [1,14.5];

if measured_divergence == 0
    thetas = [];
    for i = 1:numel(SPP_phases)
        waists = [];
        for m = 1:numel(sim_lengths)
            k = 2 * pi / lambda; 
            z_r = (pi * w_0^2) / lambda;  
            w_of_z = sqrt(2*(z_r^2 + sim_lengths(m)^2) / k * z_r); 
            waists(m) = w_of_z;
        end
        thetas(i) = 2 * atan((waists(2) - waists(1)) / (sim_lengths(2) - sim_lengths(1)))
    end
    z_eff = ((dist_aft_SPP1 * acos(thetas(1))) + ((dist_aft_SPP2) * (SPP_phases(2) / l_int) * acos((thetas(1) + (thetas(1) + thetas(2))))))


else
% Tested and validated with measured divergence
    z_eff = ((dist_aft_SPP1 * acos(theta_1)) + ((dist_aft_SPP2) * (SPP_phases(2) / l_int) * acos((theta_1 + (theta_1 + theta_2))))); 

end