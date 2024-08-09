function [hPlot, ax] = fig_create(name,size, pos, Text)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    hPlot = figure('NumberTitle' ,'off', 'Name', name,...
                   'units','centimeters','outerposition', size, 'Visible', 'on'); %
    clf(hPlot);
    set(hPlot, 'PaperUnits', 'centimeters');
    set(hPlot, 'Color', 'White');
    ax = axes('units','centimeters','XGrid','off','XMinorGrid','off','FontSize',Text.axis,'Box','on','Layer','top');
    ax.Position = pos;
    hold(ax, 'on');
    
end

