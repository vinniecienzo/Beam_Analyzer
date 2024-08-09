clc;
clear all;
close all;

%opengl software;

%% INITIALIZE SEQUENCE
% FIGURE ATTRIBUTES
fig.text = 20;
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
test_name = 'm_32_noL_834nm_test';

if ~exist([fig.dir '/' exp_name], 'dir')
    mkdir([fig.dir '/' exp_name])
end

% Make sure all files are in Working Folder

%files used for image analysis 
 files = dir([exp_name '/*.png']);

 fileCount = sum(~[files.isdir]);

 %files for testing comparison, dust, etc.
 noise_files =  dir([test_name '/*.png']);

 noise_file_count = sum(~[noise_files.isdir]);

 % disp(fileCount) % For troubke shooting purposes

 Exp_im_tmp = [];

% This loop is for main image processing  
for i = 1:fileCount

    % Reads the image as a n * m matrix 
    tmp.img = double(imread([exp_name '/' files(i).name]));
    
    % Converts RGB to Greyscale, if needed 
    %tmp.img = mean(tmp.img(:,:,1:3),3);

    % Applies soothing gaussian filter, if needed
    % tmp.img = imgaussfilt(tmp.img, 0.5);

    % Creates a list of each images n * m matrix
    Exp_im_tmp(:,:,i) = tmp.img;
end

% averages them
Exp_Im = mean(Exp_im_tmp, 3);   

noise_tmp = [];

% this is for background image processing to reduce noise like dust etc.;
% I'm not sure how it was originally employed but I'm keeping it here for
% now
% for i=1noise_file_count
% 
%     % same processing as above for test images to reduce noise
%     tmp.img = double(imread(test_name '/' files(i).name]));
%     tmp.img = mean(tmp.img(:,:,1:3),3);
%     tmp.img = imgaussfilt(tmp.img, 0.5);
%     noise_tmp(:,:,i) = tmp.img;
% %     tmp = tmp./max(tmp(:));
% %     tmp(tmp<0.09) = 0;
% end

% Attempts to remove noise
IMAGE = Exp_Im; %- mean(noise_tmp,3);   Ensure how to employ noise


% Compute the maximum value in the background_image matrix
maxValue = max(IMAGE(:));

% Normalize background_image
norm_im = IMAGE./ maxValue;

% Find the maximum value in the normalized matrix
maxNormValue = max(norm_im(:));

%Section used for Gaussian distribution

% [maxRow, maxCol] = find(norm_im == maxNormValue);
% 
% %Display the maximum normalized value and its indices, if needed
% disp(['Maximum normalized value: ' num2str(maxNormValue)]);
% disp(['Indices of maximum value: Row ' num2str(maxRow) ', Column ' num2str(maxCol)]);

%% Circle Detection; only use to find centroid 
% Convert to grayscale if needed
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

% figure;
% imagesc(cent_finder);
% title('Normalized Averaged Beam');
% xlabel('Position (px)');
% ylabel('Position (px)');
% colorbar;

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
    radius = stats(1).EquivDiameter / 2;
else
    error('No object detected in the image.');
end

% intensity profile along a horizontal line
line_y = round(centroid(1,2)) % Changes based on identity of Max Intensity
intensity_profile_row = norm_im(line_y, :);

% pixel to mm conversion
intensity_profile_row_mm = [];

for i = 1:length(intensity_profile_row)
    intensity_profile_row_mm(i) = i * .00634615;

end

% This section will need to be reevaluted for fit parameters when we
% determine how determine OAM fit.

% Fit Gaussian distribution to the intensity profile (horizontal)
% x_row = 1:numel(intensity_profile_row);
% gaussian_fit_row = fit(intensity_profile_row_mm', double(intensity_profile_row)', 'gauss1');

% Evaluate goodness of fit
%gaussian_values_row = gaussian_fit_row(intensity_profile_row_mm);


% sse_row = sum((double(intensity_profile_row) - gaussian_values_row').^2); % Sum of squared errors
% sse_tot = sum((double(intensity_profile_row) - mean(intensity_profile_row).^2));
% 
% fit_percent_row = 1 - (sse_row / sse_tot);

%Samesteps as above for the Vertical Gauss Profile

line_x = round(centroid(1,1))
intensity_profile_col = norm_im(:, line_x);

intensity_profile_col_mm = [];

for i = 1:length(intensity_profile_col)
    intensity_profile_col_mm(i) = i * .00634615;

end

% Again goodness of fit will be evaluted later

% (vertical)
% gaussian_fit_col = fit(intensity_profile_col_mm', double(intensity_profile_col), 'gauss1');
% gaussian_values_col = gaussian_fit_col(intensity_profile_col_mm);
% 
% 
% sse_col = sum((double(intensity_profile_col') - gaussian_values_col').^2); % Sum of squared errors
% sse_tot = sum((double(intensity_profile_col') - mean(intensity_profile_col)).^2);
% 
% 
% 
% fit_percent_col = 1 - (sse_col / sse_tot);

%displays averaged beam image

figure;
imagesc(norm_im);
title('Normalized Averaged Beam with Contours');
xlabel('Position (px)');
ylabel('Position (px)');
colorbar;

hold on;
line([line_x,line_x], [1, size(norm_im, 1)], 'Color', 'r', 'LineWidth', 2);
line([1, size(norm_im, 2)], [line_y,line_y], 'Color', 'r', 'LineWidth', 2);
hold off;

figure;
imagesc(norm_im);
title('Normalized Averaged Beam');
xlabel('Position (px)');
ylabel('Position (px)');
colorbar;

% Plot the results for visualization (horizontal)
figure;

subplot(1,2,1);
plot(intensity_profile_row_mm, intensity_profile_row, 'b') %,intensity_profile_row_mm, gaussian_values_row, 'r'); no gauss fit ywt
title(['X-axis Beam Profile: Intensity Profile and Gaussian Fit (Peak Measured =' num2str(line_y) ' px)']);
legend('Intensity Profile') %, 'Gaussian Fit');
xlabel('Position (mm)');
ylabel('Intensity (au)');
legend("show");

% Again fit will be accounted for later 

% text(0.05, 0.95, sprintf('Mean: %.2f\nStd Dev: %.2f\nFit: %.2f%%', gaussian_fit_row.b1, gaussian_fit_row.c1, fit_percent_row*100), ...
%     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

% Display parameters of the Gaussian fit (horizontal)
%fprintf('Horizontal Fit - Mean: %.2f, Standard deviation: %.2f\n', gaussian_fit_row.b1, gaussian_fit_row.c1);

%plot the results for visualization (vertical)
subplot(1,2,2)
plot(intensity_profile_col_mm, intensity_profile_col, 'b') %, intensity_profile_col_mm, gaussian_values_col, 'r');
title(['Y-axis Beam Profile: Intensity Profile and Gaussian Fit (Peak Measured =' num2str(line_x) ' px)']);
legend('Intensity Profile') %, 'Gaussian Fit');
xlabel('Position (mm)');
ylabel('Intensity (au)');
legend("show");

% Again fit will be accounted for later 

% text(0.05, 0.95, sprintf('Mean: %.2f\nStd Dev: %.2f\nFit: %.2f%%', gaussian_fit_col.b1, gaussian_fit_col.c1, fit_percent_col*100), ...
%     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', fig.text, 'BackgroundColor', 'w');

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
