function [iniwells, edges, wells, edgelist] = simulation1(ncolors, nrows, ncols, shape, lambda, mu, simstep_max, rvRGB, wells, rvedges, max_edges, shiftedrows, originalrows)
% SIMULATION1.m creates the wells and the contamination edges for the
% simulations. The wells before the contamination, as well as the
% wells after contamination, which are an output of contamination.m, are
% visualized.
% This uses method 1, that is, independent Bernoulli variables.
%
% EXPLANATION OF VARIABLES
%
% ncolors           Number of colours (most functions that have not been updated require 3)
% nrows             Number of rows
% ncols             Number of columns
% shape             Shape of the grid
% originalrows      Rows in wells that are 'further left'.
% shiftedrows       Rows in wells that are shifted right by a half unit.
% simstep_max       Total number of simulations
% simstep           Number of current simulation
% wells             simstep_max matrices of size (nrows,ncols,ncolors+1) -
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
% iniwells          The wells before contamination (after seeding). Needed
%                   for visualization
% edges             simstep_max matrices of size (nrows, ncols, 3)
%                   matrix whose entries (i,j,k) are indicator
%                   variables of edges between well (i,j) and its
%                   neighbor to the (k=1) right, (k=2) right down,
%                   (k=3) left down.
% lambda            Seeding rate for the ncolors colors
% mu                Contamination rate (probability of any given
%                   undirected edge being open)
%
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 16/7/2015, 25-26/4/2016, 1-3/11/2016

wells=repmat(wells,[1,1,1,simstep_max]);
edges=zeros(nrows,ncols,3,simstep_max);

for simstep=1:simstep_max
    for i=1:ncolors
        wells(:,:,1+i,simstep)=wells(:,:,1,simstep).*(rvRGB(:,:,i,simstep)<lambda(i));  
    end
end

iniwells=wells; % seeding state of vertices saved before contamination

for simstep=1:simstep_max
    for i=1:3
     edges(:,:,i,simstep)=max_edges(:,:,i).*(rvedges(:,:,i,simstep)<mu);
    end
end


for simstep=1:simstep_max

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
