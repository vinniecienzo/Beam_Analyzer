%% Read me

% Divergence checker is just a simplified code from Vortex and Gaussian
% analysis to measure a beams collimation or a SPP's measured divergence.
% Does not calculate divergence since must be measured in the experiment
% apparatus, but gives respective diameters for divergence.

% Instructions: take two measurements at two different Z lengths, save the
% figures below and use difference in diameters (or radii) and measured lengths in the below formula 
% Theta is the divergence angle in the rayleigh range.

% theta = 2 * atan((diameter_2 - diamteter_2)/(length_difference(z2 - z1));
% check your units

%% Change as needed
% Images taken from data storage to run divergence check
exp_name = 'Vortex Working Folder';

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

figure;
imagesc(norm_im);
title('Normalized Averaged Beam');
xlabel('Position (px)');
ylabel('Position (px)');
axis equal;
axis tight;
colorbar;

 text(0.05, 0.95, sprintf('L: %.2f\nP: %.2f\n radii: %.2f', l_int, p_int, radii_mm), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

%Shows the circle that cent finder is using from norm im 
hold on
viscircles(centroid,radii)

hold off