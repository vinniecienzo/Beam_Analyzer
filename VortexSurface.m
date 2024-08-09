function vortex = vortexsurface(x, y, len, w_0, p_int, l_int, lambda, cent_x, cent_y)
    % Define necessary constants and parameters
    r = sqrt((x - cent_x).^2 + (y - cent_y).^2) * 6.5 *10^-6 ; % verfied correct (px to m)
    k = 2 * pi / lambda;  % wave number
    z_r = (pi * w_0^2) / lambda;  % Rayleigh range
    z_bar = (z_r^2 + len^2) / len; % complex conjugate
    w_of_z = sqrt(2*(z_r^2 + len^2) / k * z_r); % waist of beam at l distance away; if seq SPP uses Z-eff
    psi_z = atan(len / z_r); % diveragnce angle
    constant = 1;


    % Calculate Laguerre polynomial L_l_int(2 * r^2 / w_of_z^2)
    L = zeros(size(r));
    for m = 0:l_int
        log_nchoosek = gammaln(l_int + 1) - gammaln(m + 1) - gammaln(l_int - m + 1);
        L = L + exp(log_nchoosek) * ((-1)^m / factorial(m)) * (2 * r.^2 / w_of_z^2).^m;
    end

    % Define the vortex function as described in the paper
    vortex = constant * ((-1)^p_int / w_of_z) * ((r * sqrt(2)) / w_of_z).^l_int .* L .*(2*(r/ w_of_z).^2* ...
        exp(-1i * (2 * p_int + l_int + 1) * psi_z) .* ...
        exp(-1i * (k * (r.^2) / (2 * z_bar))) .* exp(-((r.^2) / w_of_z^2)));
end
