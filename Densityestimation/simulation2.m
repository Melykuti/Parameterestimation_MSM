function [iniwells, edges, wells, edgelist] = simulation2(ncolors, nrows, ncols, shape, lambda, mu, simstep_max, wells, indexed_wells, totalwells, rpwells, max_edges, indexed_edges, totaledgelocations, rpedges, shiftedrows, originalrows)

% SIMULATION2.m creates the wells and the contamination edges for the
% simulations. The wells before the contamination, as well as the
% wells after contamination, which are an output of contamination.m, are
% visualized.
% This uses method 2, that is, random permutations.
%
% EXPLANATION OF VARIABLES
%
% ncolors           Number of colours (most functions that have not been updated require 3)
% nrows             Number of rows
% ncols             Number of columns
% shape             Shape of the grid
% originalrows      Rows in wells that are 'further left'.
% shiftedrows       Rows in wells that are shifted right by a half unit.
% simstep_max       Number of simulations
% simstep           Number of current simulation
% wells             simstep_max matrices of size (nrows,ncols,ncolours+1) -
%                   representing the grid we get from the data set - with
%                   entries of either 0 or 1.
%                   First coordinate: row
%                   Second coordinate: column
%                   Third coordinate: 1:Indicator whether well is part of
%                                       the area of interest. Ie if
%                                       wellsds.w(i,j,1)==1 that well is
%                                       going to be part of the final grid
%                                       and of the calculations. If
%                                       wellsds.w(i,j,1)==0, we ignore that
%                                       well for all our calculations.
%                   For ncolors==3,   2:Indicator if red (1=well is filled
%                                       with red seed, 0=well is not filled
%                                       with red seed)
%                                     3:Indicator of green
%                                     4:Indicator of blue
% totalwells        Number of wells that could possibly have a color.
% threshwell        Number of seeds that will be put into wells
% indexed_wells     wells(:,:,1,:) has indicator variables 0 and 1;
%                   indexed_wells will index the locations of 1s in
%                   increasing order, 0s will remain 0
% rpwells           Random permutation of totalwells
% iniwells          The wells before contamination (after seeding). Needed
%                   for visualization
% edges             simstep_max matrices of size (nrows, ncols, 3)
%                   matrix whose entries (i,j,k) are indicator
%                   variables of edges between well (i,j) and its
%                   neighbor to the (k=1) right, (k=2) right down,
%                   (k=3) left down.
% totaledgelocations Number of edges that could possibly be open.
% threshedge        Number of edges that will be put into edges.
% indexed_edges     max_edges has indicator variables 0 and 1;
%                   indexed_edges will index the locations of 1s in
%                   increasing order, 0s will remain 0.
% rpedges           Random permutation of totaledgelocations
% lambda            Seeding rate for the ncolors colors
% mu                Contamination rate (probability of any given
%                   undirected edge being open)
%
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 16/7/2015, 25-26/4/2016, 1/11/2016

wells=repmat(wells,[1,1,1,simstep_max]);
edges=zeros(nrows,ncols,3,simstep_max);

threshwell=round(lambda*totalwells);
% We use the random permutations rpwells so that the proportions of seeded wells are very accurately the respective parameters.
% rpwells(1:threshwell(i),i,simstep);  -- Gives index values; the places in indexed_wells where these values are found are those that get the seeds of colour i.
% The entries of wells with the four coordinates
% [j1 j2]=find(indexed_wells==rpwells(j,i,simstep)), i+1, simstep
% will get a seed of colour i.
for simstep=1:simstep_max
    for i=1:ncolors

        for j=1:threshwell(i)
            [j1, j2]=find(indexed_wells==rpwells(j,i,simstep));
            wells(j1,j2,i+1,simstep)=1;
        end

    end
    wells(:,:,1,simstep)=wells(:,:,1,1); % Adjust area of interest for all simulations
end
iniwells=wells; % seeding state of vertices saved before contamination


threshedge=round(mu*totaledgelocations);
% We use the random permutations rpedges so that the proportion of open edges is very accurately the respective parameter.
% rpedges(1:threshedge,simstep);  -- Gives index values; the places in indexed_edges where these values are found are those that will get the open edges.
% The entries of edges with these four coordinates:
% ind2sub([nrows,ncols,3], find(indexed_edges==rpedges(i,simstep))), simstep
% will get an edge. Here we cannot use find on its own as there are 3 coordinates, not only 2. The way we use it gives a linear index.

for simstep=1:simstep_max
    for i=1:threshedge
        [i1, i2, i3]=ind2sub([nrows,ncols,3], find(indexed_edges==rpedges(i,simstep)));
        edges(i1,i2,i3,simstep)=1;
    end
    
    edges3d=edges(:,:,:,simstep); % Transform edges into 3D to give correct input for contamination.m
    % Find all the contaminated wells
%    if simstep==1 % For visualization of the first simulation
        [edgelist, components1]=contamination(edges3d, nrows, ncols, shiftedrows, originalrows);
%    else
%        [~, components1]=contamination(edges3d, nrows, ncols, shiftedrows, originalrows);
%    end


    wells3d=wells(:,:,:,simstep); % Transform wells for each simulation into 3D
    for i=1:size(components1,2) % Loop for all connected components of contamination
        for j=1:ncolors % Loop for all colors RGB
            comp=components1{i}+j*nrows*ncols; % j*nrows*ncols because of dimension of wells => 'Jumps' into next matrix
            wells3d(comp)=max(wells3d(comp));   % Apply contamination for connected component i and color j
        end
    end
    wells(:,:,:,simstep)=wells3d; %Transform wells back into 4D
end


end
