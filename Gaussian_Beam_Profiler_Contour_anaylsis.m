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

% Improvements: Find a method to better fit the gaussian to an idealized
% beam via simulation using the same parameters; do this in a 3D manner
% using the RGB data given with the average instead of contours.

%% INITIALIZE SEQUENCE
% FIGURE ATTRIBUTES
fig.text = 12;      % Font size for text
fig.title = 24;     % Font size for titles
fig.label = 20;     % Font size for axis labels
fig.axis = 20;      % Font size for axis ticks
fig.leg = 20;       % Font size for legends
fig.save = true;    % Flag to save figures
fig.dir = 'img';    % Directory to save images

% CREATE IMG DIRECTORY TO STORE IMAGES
if ~exist(fig.dir , 'dir')
    mkdir(fig.dir)
end

% Images taken from data storage to run this code
exp_name = 'Working_folder';   
if ~exist([fig.dir '/' exp_name], 'dir')
    mkdir([fig.dir '/' exp_name])
end

% Make sure all files are in Working Folder in Dir

 files = dir([exp_name '/*.png']);

 fileCount = sum(~[files.isdir]);

 % disp(fileCount) % For troubke shooting purposes

% Loop through each image based on # no of files; ensure same pixels l * w

for i = 1:fileCount

    % Read image, convert to double !!!! ensure naming convention is solid
    tmp.img = double(imread([exp_name '/' files(i).name]));
    
    % Store processed image in BKG array
    BKG(:,:,i) = tmp.img;
end

% Compute overall background image by averaging BKG images
background_image = mean(BKG, 3);

% Compute the maximum value in the background_image matrix
maxValue = max(background_image(:));

% Normalize background_image
norm_background_im = background_image / maxValue;

% Find the maximum value in the normalized matrix
maxNormValue = max(norm_background_im(:));

[maxRow, maxCol] = find(norm_background_im == maxNormValue);

% Display the maximum normalized value and its indices, if needed
% disp(['Maximum normalized value: ' num2str(maxNormValue)]);
% disp(['Indices of maximum value: Row ' num2str(maxRow) ', Column ' num2str(maxCol)]);

% Display the entire normalized background image with color
% figure;
% imagesc(norm_background_im);
% title('Normalized Background Image; No RGB');
% colorbar;  % Display a colorbar to show intensity scale

% intensity profile along a horizontal line
line_y = maxRow; % Changes based on identity of Max Intensity
intensity_profile_row = norm_background_im(line_y, :);

% pixel to mm conversion
intensity_profile_row_mm = [];

for i = 1:length(intensity_profile_row)
    intensity_profile_row_mm(i) = i * .00634615;

end

% Fit Gaussian distribution to the intensity profile (horizontal)
x_row = 1:numel(intensity_profile_row);
gaussian_fit_row = fit(intensity_profile_row_mm', double(intensity_profile_row)', 'gauss1');

% Evaluate goodness of fit
gaussian_values_row = gaussian_fit_row(intensity_profile_row_mm);


sse_row = sum((double(intensity_profile_row) - gaussian_values_row').^2); % Sum of squared errors
sse_tot = sum((double(intensity_profile_row) - mean(intensity_profile_row).^2));

fit_percent_row = 1 - (sse_row / sse_tot);

%Samesteps as above for the Vertical Gauss Profile

line_x = maxCol; 
intensity_profile_col = norm_background_im(:, line_x);

intensity_profile_col_mm = [];

for i = 1:length(intensity_profile_col)
    intensity_profile_col_mm(i) = i * .00634615;

end

% (vertical)
gaussian_fit_col = fit(intensity_profile_col_mm', double(intensity_profile_col), 'gauss1');
gaussian_values_col = gaussian_fit_col(intensity_profile_col_mm);


sse_col = sum((double(intensity_profile_col') - gaussian_values_col').^2); % Sum of squared errors
sse_tot = sum((double(intensity_profile_col') - mean(intensity_profile_col)).^2);



fit_percent_col = 1 - (sse_col / sse_tot);

%displays averaged beam image
figure;
imagesc(norm_background_im);
title('Normalied Averaged Beam');
xlabel('Position (px)');
ylabel('Position (px)');
colorbar;


% Step 6: Plot the results for visualization (horizontal)
figure;

subplot(1,2,1);
plot(intensity_profile_row_mm, intensity_profile_row, 'b',intensity_profile_row_mm, gaussian_values_row, 'r');
title(['X-axis Beam Profile: Intensity Profile and Gaussian Fit (Peak Measured =' num2str(line_y) ' px)']);
legend('Intensity Profile', 'Gaussian Fit');
xlabel('Position (mm)');
ylabel('Intensity (au)');
legend("show");

text(0.05, 0.95, sprintf('Mean: %.2f\nStd Dev: %.2f\nFit: %.2f%%', gaussian_fit_row.b1, gaussian_fit_row.c1, fit_percent_row*100), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

% Display parameters of the Gaussian fit (horizontal)
%fprintf('Horizontal Fit - Mean: %.2f, Standard deviation: %.2f\n', gaussian_fit_row.b1, gaussian_fit_row.c1);

%plot the results for visualization (vertical)
subplot(1,2,2)
plot(intensity_profile_col_mm, intensity_profile_col, 'b', intensity_profile_col_mm, gaussian_values_col, 'r');
title(['Y-axis Beam Profile: Intensity Profile and Gaussian Fit (Peak Measured =' num2str(line_x) ' px)']);
legend('Intensity Profile', 'Gaussian Fit');
xlabel('Position');
ylabel('Intensity');
legend("show");

text(0.05, 0.95, sprintf('Mean: %.2f\nStd Dev: %.2f\nFit: %.2f%%', gaussian_fit_col.b1, gaussian_fit_col.c1, fit_percent_col*100), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

% Can save automatically, but I opt to save manually after resizing in
% MATLAB

% Display parameters of the Gaussian fit (vertical)
%fprintf('Vertical Fit - Mean: %.2f, Standard deviation: %.2f\n', gaussian_fit_col.b1, gaussian_fit_col.c1);

% % Optional: Save the plots to files
% saveas(gcf, 'Intensity_Profile_Gaussian_Fit_Vertical.png');
% 
% % Save the plots to files
% if fig.save
%     saveas(gcf, [fig.dir '/' exp_name '/Intensity_Profile_Gaussian_Fit.png']);
% end