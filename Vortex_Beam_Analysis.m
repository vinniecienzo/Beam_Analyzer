clc;
clear all;
close all;

%% READ ME

% Program uses directory to find a working folder and images that are
% stored there, turns the images into matrices, averages the matrices. The
% average matrix is used to find the absolute max intensity. From the max
% intensity the image matrix is normalized.

% As this is a beam with orbital angular momentum the centroid is the
% center of the circle, a diffraction filter is used to clearly define the
% circle to find the centroid and radius. A vortex surface from L. Allens 
% 1992 paper is fit over the same position to be compared for a fit for
% intensity and radius

% IMPORTANT: when generating a ideal vortex surface, if the total measured
% intensity is within (inside) the ideal surface it was considered acceptable even
% if the SSE gave a bad fit; especially if the radii fit test was really good.
% essentially the scaling of the complex scalar matrix of E_field predicts
% a larger FWHM, it's beleieved we have a beam with smaller FWHM and higher
% Q with an ideal radius.

% Following the discovery of weighted divergence of the OAM beam, a new feature
% function of weighted divergence was created; this function is only for
% two SPP, I confirmed this works but the full theory hasn't been validated
% due to experimental limitations. Further work on this code should include
% higher order terms in a function form for multiple SPPs and LG
% transformations to allow for a more precise numerical validation of the theory
%
% Therefore, I have included a draft function in the directory, but it hasnt been evaluated
% yet.

% any questions on this program email vmg6966@nyu.edu

%% Instructions 

% Drop PNGs in respective folder "Vortex Working Folder", adjust the input
% parameters below for respective beam, if using multiple SPPs use weighted
% divergence (you may need to do some  with code a bit), run program.

% if there seems to be an error likely it is in the experimental set-up
% check there first (normally the pin hole), this is the whole point of
% comparing the numerical simulation to the generated vortex beam.

% if not there, it's likely cent_finders where circles are taken for the
% nomarlized image and the simulation; try adjusting the diffraction filter strength, it bites of on small circles.

% If not there check the weighted divergence and associated measurements
% (this is how I discovered the theory, still not 100% veriefied)

%% Input Vortex Parameters (Known or Measured)
p_int = 0; % typically zero
l_int = 64; % total OAM phase after all transformations 
lambda = 834.7 * 10^-9; % wavelength of beam in meters
w_0 = 344 * 10^-6; % this is the beams waist at first LG transformation, in meters (measured 343/344 for summer 2024 experiment)
measured_len = 9.5; % total length from first transformation to beam profiler, in inches

%% Weighted Diveregnce Calculations for Linear Combintions of SPP
% measured diveregnce in this case; in the future should be able to use
% Allens w(z) function to generate a truly ideal vortex as a consequnce of
% SPPs

% if not measured use zero to calculate thoeretical divergence.
measured_divergence =  2 * atan((2.43-1.51)/(14.75*25.4)); % CAO 18 Jul 24; pinhole affects divergence appartently :(

theta_1 = 2 * atan((2.43-1.51)/(14.75*25.4));
theta_2 = theta_1;
dist_aft_SPP1 = 1.25; 
dist_aft_SPP2 = 8.25;
SPP_phases = [32,32]; % inital phase shift that scales amplitude

z_eff = TwoSPP_Weighted_Div_Func(lambda ,w_0, l_int, SPP_phases, dist_aft_SPP1, dist_aft_SPP2, theta_1, theta_2, measured_divergence);


% effective z component added to measured lengh to account for additional divergence as a ratio of OAM azimuthal component.
len = (measured_len + z_eff) * 2.54 * 10^-2; %m (inch convert cm) equation uses from first SPP total distance, should be fixed

%% Extra Features

% % Just for reference to check where you place beam profiler
% rayleigh_range = pi * w_0^2 / lambda;

diffraction_filt_out = 10; % in px, diffraction filter buffer outside vortex, adjust as needed for fit 
diffraction_filt_in = 75; % in px, diffraction filter buffer inside vortex, adjust as needed for fit 
diff_filt_strength = 0.4; %detemines strength of filter needed to remove diffraction to find circle

%% IMAGE PROCESSING
% FIGURE ATTRIBUTES
fig.text = 8;
fig.title = 24;
fig.label = 20;
fig.axis = 20;
fig.leg = 20;
fig.save = true; % save figure or not
fig.dir = 'img';

% CREATE IMG DIRECTORY TO STORE IMAGES
if ~exist(fig.dir , 'dir')
    mkdir(fig.dir)
end

% Images taken from data storage to run this code
exp_name = 'Vortex Working Folder';
% test_name = 'm_32_noL_834nm_test';

if ~exist([fig.dir '/' exp_name], 'dir')
    mkdir([fig.dir '/' exp_name])
end

% Make sure all files are in Working Folder

%files used for image analysis 
 files = dir([exp_name '/*.png']);

 fileCount = sum(~[files.isdir]);

 Exp_im_tmp = [];

% This loop is for main image processing  
for i = 1:fileCount

    % Reads the image as a n * m matrix 
    tmp.img = double(imread([exp_name '/' files(i).name]));
    
    % Converts RGB to Greyscale, if needed 
    %tmp.img = mean(tmp.img(:,:,1:3),3);

    % Applies filter, if needed
    % tmp.img = imgaussfilt(tmp.img, 0.5);

    % Creates a list of each images as a matrix, each px is a data point in
    % the matrix
    Exp_im_tmp(:,:,i) = tmp.img;
end

% averages the matrixes
Exp_Im = mean(Exp_im_tmp, 3);   

% Compute the maximum value in the background_image matrix
maxValue = max(Exp_Im(:));

% Normalize background_image
norm_im = Exp_Im./ maxValue; 

% populates the variable for later use
norm_im_no_diff = norm_im;

% Find the max value in the norm matrix, sanity check its one i guess
maxNormValue = max(norm_im(:));

%% CIRCLE DETECTIOION; USED TO CENTER SIMULATION
% Convert to grayscale if needed
if size(norm_im, 3) == 3
    cent_finder = rgb2gray(norm_im);

else cent_finder = norm_im;
end

% Gaussian filter & pixel dropper, if needed (adjust std-dev for smoothing)
cent_finder = imgaussfilt(cent_finder, 5);
cent_finder(cent_finder < diff_filt_strength) = 0;

% Binarize the image
bw = imbinarize(cent_finder);

% Fill holes and remove small objects
bw = imfill(bw, 'holes');
bw = bwareaopen(bw, 50);

% Detect the object
stats = regionprops(bw, 'Centroid', 'EquivDiameter');

% Calculate centroid and radius
if ~isempty(stats)
    centroid = stats(1).Centroid;
    radii = stats(1).EquivDiameter / 2;
else
    error('No object detected in the image.');
end

radii_mm = radii * 0.0065;

% Had trouble with region stats... double check your gaussian filter so it
% doesn't look for a little circles; if it gives an array of multiple
% centroids, broaden the filter. This figure section allows you to to
% verify theres no little cricles ruining your code, important for
% centering you simulated vortex

% figure;
% imagesc(norm_im);
% title('Normalized Averaged Beam');
% xlabel('Position (px)');
% ylabel('Position (px)');
% axis equal;
% axis tight;
% colorbar;
% 
%  text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f\n radii: %.2f', l_int, p_int, radii_mm), ...
%     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

% %Shows the circle that cent finder is using from norm im 
% hold on
% viscircles(centroid,radii)
% 
% hold off

%% Build Vortex Surface

% intensity profile along a horizontal line (x values for set, y value)
% Kept a profile with mm if wanted to employ it in an image, not crucial

cent_y = floor(centroid(1,2)); % Changes based on identity of Max Intensity
intensity_profile_vert = norm_im(cent_y, :);

% pixel to mm conversion
intensity_profile_row_mm = [];

for i = 1:length(intensity_profile_vert)
    intensity_profile_row_mm(i) = i * .0065;

end

%Samesteps as above for the Vertical Gauss Profile

cent_x = floor(centroid(1,1)); 
intensity_profile_horz = norm_im(:, cent_x);

intensity_profile_horz_mm = [];

for i = 1:length(intensity_profile_horz)
    intensity_profile_horz_mm(i) = i * .0065;

end

%establishes pixels grid for function
height = numel(intensity_profile_horz);
width = numel(intensity_profile_vert);
[x, y] = meshgrid(1:width, 1:height);

% Main Function which generates the surface
sim_vortex = VortexSurface(x, y, len, w_0, p_int, l_int, lambda, cent_x, cent_y);

sim_vortex = imgaussfilt(abs(sim_vortex),1);

% Gen sim for radius measurement 
max_sim_vortex = max(max(abs(sim_vortex))); % must do twice for row and col

norm_sim_vortex = abs(sim_vortex) / max_sim_vortex;

cent_finder2 = norm_sim_vortex;

% find the sims center and radii to compare
cent_finder2(cent_finder2 < diff_filt_strength) = 0;

bw2 = imbinarize(cent_finder2);

stats2 = regionprops(bw2, 'Centroid', 'EquivDiameter');
sim_centroid = stats2(1).Centroid;
sim_radii = stats2(1).EquivDiameter;
sim_radii_mm = sim_radii * 0.0065;

% % check cent finder / sims if needed
% figure;
% mesh(abs(norm_sim_vortex));
% title('Sim Vortex Beam');
% xlabel('Position (px)');
% ylabel('Position (px)');
% zlabel('Intensity')
% axis tight;
% axis equal;
% colorbar;
% 
% %Shows the circle that cent finder is using from norm im 
% hold on
% viscircles(centroid,radii)
% 
% hold off

%% Correcting Norm Im for Diffraction
% Fileters norm vortex image for diffraction noise

r_0 = sqrt(((x - cent_x)).^2 + ((y - cent_y)).^2);

norm_im_no_diff(r_0 > radii + diffraction_filt_out) = 0; % is a pixel buffer filtering out diffration
norm_im_no_diff(r_0 < radii - diffraction_filt_in) = 0; % inside the radius for noise 

%% Fit Calculations

% intensity profile along a horizontal line
line_y = round(centroid(1,2)); % Changes based on identity of Max Intensity
intensity_profile_row = norm_im(line_y, :);

% pixel to mm conversion
intensity_profile_row_mm = [];

for i = 1:length(intensity_profile_row)
    intensity_profile_row_mm(i) = i * .0065;

end

%Samesteps as above for the Vertical Gauss Profile

line_x = round(centroid(1,1));
intensity_profile_col = norm_im(:, line_x);

intensity_profile_col_mm = [];

for i = 1:length(intensity_profile_col)
    intensity_profile_col_mm(i) = i * .0065;

end

% Evaluate goodness of fit (2D)
sse_surface = sum((double(norm_im_no_diff(:)) - double(abs(norm_sim_vortex(:)))).^2); % Sum of squared errors for the surface
sse_tot = sum((double(norm_im_no_diff(:)) - mean(norm_im_no_diff(:))).^2); % Total sum of squares

% Goodness of fit
ideal_vortex_fit = 1 - (sse_surface / sse_tot);

radii_fit = radii / sim_radii; % radii given in pixels, converts to mm then fits

% Display the fit percentage
%disp(['Goodness of fit (2D Gaussian): ', num2str(fit_percent_2D_Gauss * 100), '%']);

%displays averaged beam image
%Side by side comparison of simlated beam and averaged exp beam
%% Figures 

% These are built to build awareness of your vortex and how ideal it is; I
% wouldn't use or save all these images.

figure;
subplot(1,2,1);
imagesc( intensity_profile_row_mm, intensity_profile_col_mm, norm_im);
title('Normalied Averaged Beam');
axis equal;
axis tight;
xlabel('Position (mm)');
ylabel('Position (mm)');
colormap('parula');
colorbar;

% %Shows the circle that cent finder is using from norm im 
% hold on
% viscircles(centroid,radii)
% 
% hold off

 text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f\n radii: %.2f,\n radii fit: %.2f', l_int, p_int, radii_mm, radii_fit), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

subplot(1,2,2);
imagesc(intensity_profile_row_mm, intensity_profile_col_mm, abs(norm_sim_vortex));
colormap('parula');
colorbar;
axis equal;
axis tight;
title('Simulated Vortex Beam Intensity');
xlabel('Position (mm)');
ylabel('Position (mm)');

 text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f\nideal radii: %.2f', l_int, p_int, sim_radii_mm), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

 % figure comparison of norm in no diffraction used for the 
 figure;
subplot(1,2,1);
imagesc(intensity_profile_row_mm, intensity_profile_col_mm, norm_im_no_diff);
title('Normalied Averaged Beam; Diffraction Filter');
axis equal;
axis tight;
xlabel('Position (mm)');
ylabel('Position (mm)');
colormap('parula');
colorbar;

 text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f\n radii: %.2f,\n radii fit: %.2f', l_int, p_int, radii_mm, radii_fit), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

subplot(1,2,2);
imagesc(intensity_profile_row_mm, intensity_profile_col_mm, abs(norm_sim_vortex));
colormap('parula');
colorbar;
title('Simulated Vortex Beam Intensity');
axis equal;
axis tight;
xlabel('Position (mm)');
ylabel('Position (mm)');

 text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f\nideal radii: %.2f', l_int, p_int, sim_radii_mm), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');


figure;
mesh(abs(norm_sim_vortex), 'FaceAlpha', 0.1, 'EdgeAlpha', 0.2, 'FaceColor', 'interp');
title('Ideal Vortex Intensity Surface per L. Allen');
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Intensity');
colorbar; 
colormap("parula"); %cbrewer not working?

text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f', l_int, p_int), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

figure;
mesh(abs(norm_im), 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'FaceColor', 'interp');
title('Measured Vortex Intensity Surface');
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Intensity');
colorbar; 
colormap("parula"); %cbrewer not working?

text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f', l_int, p_int), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

% Simulated Vortex Surface on the filtered data set; adjust transculnce 
figure;
[x, y] = meshgrid(1:width, 1:height);
mesh(norm_im_no_diff,'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', 'interp');

hold on;

% Plot the second surface with colormap 'parula'
mesh(abs(norm_sim_vortex), 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'FaceColor', 'interp'); 
title('Ideal Vortex Surface overlaid on top of a Normalized Averaged Vortex');
xlabel('X (px)');
ylabel('Y (px)');
axis tight;
zlabel('Intensity');
colormap('winter'); % fix cbrewer
colorbar;

text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f,\nPercentage Ideal Fit: %.2f', l_int, p_int, ideal_vortex_fit), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

hold off;

% Can save automatically, but I opt to save manually after checking
% displayed figs in MATLAB; setup as needed

% % Optional: Save the plots to files
% saveas(gcf, 'Intensity_Profile_Gaussian_Fit_Vertical.png');
% 
% % Save the plots to files
% if fig.save
%     saveas(gcf, [fig.dir '/' exp_name '/Intensity_Profile_Gaussian_Fit.png']);
% end
