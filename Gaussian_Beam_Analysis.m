clc;            % Clear command window
clear all;      % Clear all variables from workspace
close all;      % Close all figure windows

% Program uses directory to find a working folder and images that are
% stored there, turns the images into matrices, averages the matrices. The
% average matrix is used to find the absolute max intensity. From the max
% intensity the image matrix is normalized and a X and Y contour is taken
% to find the intensity distribution. The intensity dist is then compared
% to a gaussian fit over the same position and compared for a fit. The
% contours are plotted side by side along with the avergared beam image.

% Instruction: Drop PNGs in "Gaussian Working Folder" in directory (create
% if not present), run program.

% any questions on this program email vmg6966@nyu.edu

%% INITIALIZE SEQUENCE
% FIGURE ATTRIBUTES
fig.text = 12;      % Font size for text
fig.title = 24;     % Font size for titles
fig.label = 20;     % Font size for axis labels
fig.axis = 20;      % Font size for axis ticks
fig.leg = 20;       % Font size for legends
fig.save = true;    % Flag to save figures
fig.dir = 'img';    % Directory to save images

% exp_name = "Diveragance Test Position 2: 27 in";

% CREATE IMG DIRECTORY TO STORE IMAGES
if ~exist(fig.dir , 'dir')
    mkdir(fig.dir)
end

% Images taken from data storage to run this code
folder_name = 'Gaussian Working Folder';
%test_name = 'm_32_noL_834nm_test';

if ~exist([fig.dir '/' folder_name], 'dir')
    mkdir([fig.dir '/' folder_name])
end

% Make sure all files are in Working Folder

%files used for image analysis 
 files = dir([folder_name '/*.png']);

 fileCount = sum(~[files.isdir]);

  % %Noise/test comparison directory if needed
 % %files for testing comparison, dust, etc.
 % noise_files =  dir([test_name '/*.png']);
 % 
 % noise_file_count = sum(~[noise_files.isdir]);
 % 
 % % disp(fileCount) % For troubke shooting purposes

%% Image Processing
 Exp_im_tmp = [];

% Loop through each image based on # no of files; ensure same pixels l * w

for i = 1:fileCount

    % Read image, convert to double !!!! ensure naming convention is solid
    tmp.img = double(imread([folder_name '/' files(i).name]));

    %RGB Conversion if needed
    %tmp.img = mean(tmp.img(:,:,1:3),3)

     % Applies soothing gaussian filter, if needed
    % tmp.img = imgaussfilt(tmp.img, 0.5);
    
    % Store processed image in BKG array
    Exp_im_tmp(:,:,i) = tmp.img;
end

%noise_tmp = [];

% this is for background image processing to reduce noise like dust etc.;
% This was in the original verision of the code and not sure how it was employed 
% but I'm keeping it here as I'm sure there is way to benefit from this.

% for i=1:noise_file_count
% 
%     % same processing as above for test images to reduce noise
%     tmp.img = double(imread(test_name '/' files(i).name]));
%     tmp.img = mean(tmp.img(:,:,1:3),3);
%     tmp.img = imgaussfilt(tmp.img, 0.5);
%     noise_tmp(:,:,i) = tmp.img;
% %     tmp = tmp./max(tmp(:));
% %     tmp(tmp<0.09) = 0;
% end

% Compute overall background image by averaging EXP images subtracting from
% noise, noise not accounted for but could be, somehow, someway?

Exp_image = mean(Exp_im_tmp, 3); %- mean(noise_temp, 3);

% Compute the maximum value in the background_image matrix

maxValue = max(Exp_image(:));

% Normalize background_image
norm_im = Exp_image / maxValue;

imwrite(norm_im, 'Working Folder\average_normalized_image.png');

% Find the maximum value in the normalized matrix; sanity check should be 1
maxNormValue = max(norm_im(:));

if size(norm_im, 3) == 3
    cent_finder = rgb2gray(norm_im);

else cent_finder = norm_im;
end

% Gaussian filter if needed (adjust std-dev for smoothing)
cent_finder = imgaussfilt(cent_finder, 10);

% Had trouble with region stats... double check your gaussian filter so it
% doesn't look for a little circles; if it gives an array of multiple
% centroids, broaden the filter. This figure section allows you to to
% verify theres no little cricles ruining your code

% Binarize the image
bw = imbinarize(cent_finder);

% Fill holes and remove small objects
bw = imfill(bw, 'holes');
bw = bwareaopen(bw, 50);

% Detect the object
stats = regionprops(bw, 'EquivDiameter', 'Circularity', 'Centroid');

% Calculate centroid and radius
if ~isempty(stats)
    diameter = stats(1).EquivDiameter * .00634615;
    circularity = stats(1).Circularity;
    radii = stats(1).EquivDiameter /2;
    centers = stats(1).Centroid;
else
    error('No object detected in the image.');
end

% figure;
% imagesc(norm_im);
% title('Normalized Averaged Beam');
% xlabel('Position (px)');
% ylabel('Position (px)');
% colorbar;


%% Building an and Goodness of Fit with SSE

[maxRow, maxCol] = find(norm_im == maxNormValue);

% Display the maximum normalized value and its indices, if needed
% disp(['Maximum normalized value: ' num2str(maxNormValue)]);
% disp(['Indices of maximum value: Row ' num2str(maxRow) ', Column ' num2str(maxCol)]);

% intensity profile along a horizontal line (x values for set, y value)
line_y = maxRow; % Changes based on identity of Max Intensity
intensity_profile_row = norm_im(line_y, :);

% pixel to mm conversion
intensity_profile_row_mm = [];

for i = 1:length(intensity_profile_row)
    intensity_profile_row_mm(i) = i * .00634615;

end

%Samesteps as above for the Vertical Gauss Profile

line_x = maxCol; 
intensity_profile_col = norm_im(:, line_x);

intensity_profile_col_mm = [];

for i = 1:length(intensity_profile_col)
    intensity_profile_col_mm(i) = i * .00634615;

end

height = numel(intensity_profile_col);
width = numel(intensity_profile_row);
[x, y] = meshgrid(1:width, 1:height);


% Reshape the grid and the normalized image into column vectors
x_data = x(:);
y_data = y(:);
xy_data = [x_data, y_data];

% Define the custom 2D Gaussian model
gauss2D_model = fittype(@(x0, y0, sigma_x, sigma_y, x, y) ...
    gaussian2D(x0, y0, sigma_x, sigma_y, x, y), ...
    'independent', {'x', 'y'}, 'dependent', 'z');

% Improved initial guesses for the parameters
initial_guess = [mean(x_data), mean(y_data), std(x_data)/2, std(y_data)/2];

% Set fitting options to improve convergence
fit_options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'StartPoint', initial_guess, ...
                         'Lower', [0, 0, 0, 0], ...
                         'Upper', [Inf, Inf, Inf, Inf], ...
                         'MaxIter', 1000, ...
                         'TolFun', 1e-6);

% Fit the model to the data
[fitresult, gof] = fit([x_data, y_data], norm_im(:), gauss2D_model, fit_options);

% Evaluate the fit for plotting
Ideal_Gauss_Surf = fitresult(x, y);

% Use to see the reuslts of the iterative fit
% disp(fitresult);

% Evaluate goodness of fit (2D)
sse_surface = sum((double(norm_im(:)) - double(Ideal_Gauss_Surf(:))).^2); % Sum of squared errors for the surface
sse_tot = sum((double(norm_im(:)) - mean(norm_im(:))).^2); % Total sum of squares

% Goodness of fit
fit_percent_2D_Gauss = 1 - (sse_surface / sse_tot);
avg_std_dev = sqrt(fitresult.sigma_x * fitresult.sigma_y) * .00634615;

% Display the fit percentage
%disp(['Goodness of fit (2D Gaussian): ', num2str(fit_percent_2D_Gauss * 100), '%']);

%% Plotting 

%displays averaged beam image
figure;
imagesc(norm_im);
title('Normalied Averaged Beam');
xlabel('Position (px)');
ylabel('Position (px)');
colormap('parula')
colorbar;
axis equal;
axis tight;
text(0.05, 0.95, sprintf('Diameter (mm): %.2f\nCircularity: %.2f', diameter, circularity), ...
     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

% Shows the found circle
hold on
viscircles(centers,radii)

hold off


%Display the entire normalized data
figure;
mesh(norm_im, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'FaceColor', 'interp');
title('Normalized Background Image');
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Intensity');
colorbar;  % Display a colorbar to show intensity scale
colormap("parula") %changes color of of actual data for visualization

% Plot the Gaussian surface
figure;
mesh(Ideal_Gauss_Surf, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'FaceColor', 'interp'); 
title('2D Ideal Gaussian Surface');
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Intensity');
colormap('winter');  % Change colormap for the second surface
colorbar;

% fitted Surface on the actual data set; adjust transculnce 
figure;
[x, y] = meshgrid(1:width, 1:height);
mesh(norm_im,'FaceAlpha', 0.15, 'EdgeAlpha', 0.15, 'FaceColor', 'interp');

hold on;

% Plot the second surface with colormap 'parula'
mesh(Ideal_Gauss_Surf, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'FaceColor', 'interp'); 
title('2D Ideal Gaussian Surface overlaid on top of a Normalized Averaged Beam');
xlabel('X (px)');
ylabel('Y (px)');
zlabel('Intensity');
colormap('winter');  % Change colormap for the second surface
colorbar;

text(0.05, 0.95, sprintf('2D Surface Avg Std Dev: %.2f\nFit: %.2f%%', avg_std_dev, fit_percent_2D_Gauss*100), ...
     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');
hold off;

%% Figure saving

% Can save automatically, but I opt to save manually after resizing in
% MATLAB

% Display parameters of the Gaussian fit (vertical)
%fprintf('Vertical Fit - Mean: %.2f, Standard deviation: %.2f\n', gaussian_fit_col.b1, gaussian_fit_col.c1);

% % Optional: Save the plots to files
% saveas(gcf, 'Intensity_Profile_Gaussian_Fit_Vertical.png');
% 
% % Save the plots to files
% if fig.save
%     saveas(gcf, [fig.dir '/' folder_name '/' exp_name '.png']);
% end
