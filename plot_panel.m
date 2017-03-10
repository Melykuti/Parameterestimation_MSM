function plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, wells_on, edges_on)

% PLOT_PANEL.m is called by PLOT_FIGURE.m to plot a single subplot in a figure
% with possibly multiple panels. 
%
% INPUTS
%
% wells_on==0 displays tiny thin white outlines with black interior for vertices
% wells_on==1 displays balls of color with identical outline
% wells_on==2 displays thin white outlines with black interior for empty
% wells and balls of color with identical outline
%
% edges_on==0 displays no edges
% edges_on==1 displays edges underneath filled circles
% edges_on==2 displays edges on top of filled circles
%
% Permitted combinations: (wells_on, edges_on) == (0,1), (1,0), (2,0), (2,1), (2,2)
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 26-27/4/2016, 7/2/2017

period_xtick=5; % period of labelled rows and columns
period_ytick=period_xtick;

d=1.5; % distance between centres of circles, when the diameter of a circle is 1
d_sqrt32=d*sqrt(3)/2;
gc=[0.5 0.5 0.5]; % gray color

axis equal
axis([d_sqrt32-d d_sqrt32*nrows+1+d 0 d*(ncols+0.5)+1+d])

set(gca,'color','k')
% Labelling rows:
set(gca,'XTick',d_sqrt32*(1:period_xtick:nrows)+0.5, 'XTickLabel', 1:period_xtick:nrows)
% Labelling columns:
set(gca,'YTick',d*(1:period_ytick:(ncols+1))+0.5, 'YTickLabel', 1:period_ytick:(ncols+1))


view(90,90)
hold on
if wells_on==0 & edges_on==1
    plot_edges(nrows, ncols, originalrows, shiftedrows, d, d_sqrt32, experiment_edgelist, [1 1 1])
    plot_wells(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, d, d_sqrt32, [1 1 1], 0)
    %plot_edges(nrows, ncols, originalrows, shiftedrows, d, d_sqrt32, experiment_edgelist, gc)
end

if wells_on==1 & edges_on==0
    plot_wells(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, d, d_sqrt32, gc, 1)
end

if wells_on==2 & edges_on==0
    plot_wells(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, d, d_sqrt32, gc, 2)
end

if wells_on==2 & edges_on==1
    plot_edges(nrows, ncols, originalrows, shiftedrows, d, d_sqrt32, experiment_edgelist, [1 1 1])
    plot_wells(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, d, d_sqrt32, gc, 2)
end

if wells_on==2 & edges_on==2
    plot_wells(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, d, d_sqrt32, gc, 2)
    plot_edges(nrows, ncols, originalrows, shiftedrows, d, d_sqrt32, experiment_edgelist, [1 1 1])
end

hold off
