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
