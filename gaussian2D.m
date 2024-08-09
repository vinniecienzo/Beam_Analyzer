function z = gaussian2D(x0, y0, sigma_x, sigma_y, x, y)
    z = 1 * exp(-((x - x0).^2 / (2 * sigma_x^2) + (y - y0).^2 / (2 * sigma_y^2)));
end