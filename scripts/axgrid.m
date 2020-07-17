function [ ax ] = axgrid(figh, figw, mb, mt, ml, mr, rows, cols, units)
%AXGRID Return axis positions that fit a given area as well as possible.
% Takes figure height, figure width, margin height bottom, margin height top, 
% margin width left, margin width right, rows, columns, and units to use 
% ('pixels', 'inches', 'centimeters', 'normalized').

height = (figh - rows*(mt+mb)) / rows; % height per subplot, based on total figure size
width = (figw - cols*(ml+mr)) / cols; % width per subplot, based on total figure size

ax = gobjects(cols, rows);
for i = 1:cols
    for j = 1:rows
        left = ml + (i-1)*width + (i-1)*(ml+mr);
        bottom = mb + (j-1)*height + (j-1)*(mt+mb);
        ax(i,j) = axes('units', units, 'Position', [left bottom width height]);
    end
end

% Rearrange so ax(1) is plot 1, ax(2) is plot 2, etc., from left to right
% and from up to down
ax_temp = ax;
k = 0;
for j = rows:-1:1
    for i = 1:cols
        k = k+1;
        ax(k) = ax_temp(i,j);
    end
end

end

