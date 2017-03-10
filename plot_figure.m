function plot_figure(n)

% PLOT_FIGURE.m creates the two figures of the paper and provides a general
% tool for plotting our datasets of triangular lattices with colored vertices
% and randomly open edges.
%
% INPUTS
%
% n==1 creates Plot 1 of paper
% n==2 creates Plot 2 of paper
% n==3 creates a modified Plot 2 with panels over each other, not side by side
%
% n==10 displays a single panel
% n==11 displays two panels side by side: seeding (i.e. before contamination)
%  and after contamination
% n==12 is similar to n==11 but for the case of missing variable edgelist. It
%  first creates edgelist from the matrix open_edges using EDGEMATRIX_TO_EDGELIST.m
%  at the bottom of this file.
%
% Save image with e.g.
% >> print('data_0.png','-dpng', '-r300')
% >> set(gcf, 'PaperSize', [9 6.5]); print('dat8x12_ver1.pdf','-dpdf','-bestfit')
% >> print('data_0.eps','-depsc')
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 26-28/4/2016, 30/1-7/2/2017

figure;
set(gcf, 'InvertHardCopy', 'off', 'Color', [1 1 1]); % to save black panel background but prevent gray figure background
switch n
    case 1
        mu=0; % We need to initialise mu as a variable before loading it because it is also a function name.
        load('dat8x12.mat')
        %load('data40x30_highcont.mat')
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(experiment_matrix, shape);
        subplot(2,2,1);plot_panel(experiment_seeds, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 0);
        subplot(2,2,2);plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 0, 1);
        subplot(2,2,3);plot_panel(experiment_seeds, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 1);
        subplot(2,2,4);plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 0);

    case 2
        mu=0; % We need to initialise mu as a variable before loading it because it is also a function name.
        load('dat40x30_lowcont.mat')
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(experiment_matrix, shape);
        subplot(1,2,1);plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 1, 0);

        load('dat40x30_highcont.mat')
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(experiment_matrix, shape);
        subplot(1,2,2);plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 1, 0);

    case 3
        mu=0; % We need to initialise mu as a variable before loading it because it is also a function name.
        load('dat40x30_lowcont.mat')
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(experiment_matrix, shape);
        subplot(2,1,1);plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 1, 0);

        load('dat40x30_highcont.mat')
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(experiment_matrix, shape);
        subplot(2,1,2);plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 1, 0);

    case 10
        %load('dat8x12.mat')
        load('Data/25x25_1color_m2.mat'); experiment_matrix=wells_after_contamination;
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(experiment_matrix, shape);
        plot_panel(experiment_matrix, nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 1)

    case 11
        load('test4.mat');
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(wells_after_contamination(:,:,:,1), shape);
        subplot(1,2,1);plot_panel(wells_before_contamination(:,:,:,1), nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 1)
        subplot(1,2,2);plot_panel(wells_after_contamination(:,:,:,1), nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 1)
        
    case 12 % when there are no edges saved in experiment_edgelist, only in open_edges
        load('25x25.mat');
        [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(wells_after_contamination(:,:,:,1), shape);
        experiment_edgelist=edgematrix_to_edgelist(open_edges, originalrows, shiftedrows);
        subplot(1,2,1);plot_panel(wells_before_contamination(:,:,:,1), nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 1)
        subplot(1,2,2);plot_panel(wells_after_contamination(:,:,:,1), nrows, ncols, ncolors, originalrows, shiftedrows, experiment_edgelist, 2, 1)
end
end


function edgelist=edgematrix_to_edgelist(edges, originalrows, shiftedrows)

% This originates from CONTAMINATION.m
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 12/8/2015, 7/2/2017

edgelist=zeros(4,0);

% For right arrows
[i, j]=find(edges(:,:,1)); % Find row and column indices of all nonzero elements of edges(:,:,1)
i=i(:); j=j(:);
% 1st coordinate: row indices of 1st vertices of edges
% 2nd coordinate: column indices of the same vertices
% 3rd coordinate: row indices of 2nd vertices of edges
% 4th coordinate: column indices of the same vertices
edgelist=[edgelist [i';j';i';j'+1]];

% Analogue for right-down arrows
[i, j]=find(edges(shiftedrows(:),:,2)); i=i(:); j=j(:);
edgelist=[edgelist [shiftedrows(i');j';shiftedrows(i')+1;j'+1]];
[i, j]=find(edges(originalrows(:),:,2)); i=i(:); j=j(:);
edgelist=[edgelist [originalrows(i');j';originalrows(i')+1;j']];

% Analogue for left-down arrows
[i, j]=find(edges(shiftedrows(:),:,3)); i=i(:); j=j(:);
edgelist=[edgelist [shiftedrows(i');j';shiftedrows(i')+1;j']];
[i, j]=find(edges(originalrows(:),:,3)); i=i(:); j=j(:);
edgelist=[edgelist [originalrows(i');j';originalrows(i')+1;j'-1]];

end
