function shiftedcolor = plot_colorshift(color)

% PLOT_COLORSHIFT.m converts a 1, 2 or 3-dimensional indicator of vertex
% occupancy into a 3-dimensional RGB value. This function also allows the
% permutation of color channels or the use of false colors.
%
% VARIABLES
%
% ncolors        Permitted values: 1, 2, 3.
% color          ncolor-dimensional column vector, indicator variables
%                representing the coloring (i.e. occupancy) of a single vertex.
% shiftedcolor   3-dimensional row vector representing a color in RGB coding,
%                each element is in [0,1].
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 7/2/2017

ncolors=length(color);
switch ncolors
    case 1 % The color of plotting is blue.
        shiftedcolor=[0 0 color];
    case 2 % The colors of plotting are green and blue.
        shiftedcolor=[0 color'];
    case 3
        shiftedcolor=color'; % There is no shift.
        %shiftedcolor=[color(3) color(1) color(2)]; % False colors.
    otherwise % For more than 3 colors, we can try plotting any 3 color channels, with some blacking out effects when unplotted colors are present.
        shiftedcolor=color([1,2,4])';
end
end
