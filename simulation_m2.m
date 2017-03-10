function [inivertices, vertices, edges, edgelist] = simulation_m2(nrows, ncols, ncolors, shape, originalrows, shiftedrows, lambda, mu, simstep_max, max_vertices, max_edges, totalvertices, totaledgelocations, rplabelled_colors, rplabelled_edges)

% SIMULATION_M2.m creates the seeding, the contamination edges and the resulting 
% synthetic dataset by simulation.
%
% This uses method 2, that is, random permutations.
%
% EXPLANATION OF VARIABLES
%
% inivertices       Size: (nrows, ncols, 1+ncolors, simstep_max), essentially
%                   simstep_max independent random realizations of the seeding,
%                   before contamination would have taken place. Primary use
%                   is in visualization. See also vertices.
% vertices          Size: (nrows, ncols, 1+ncolors, simstep_max), essentially
%                   simstep_max independent random realizations of the dataset.
%                   They store the synthetic dataset, similarly to how
%                   experiment_matrix stores the experimental data.
%                   vertices(:,:,1,simstep) are indicators of the area of interest. 
%                   (These are identical for simstep=1:simstep_max.)
%                   vertices(i,j,2:ncolors+1) represent the presence 
%                   of individual colors in the dataset in the given vertex (i,j).
% edges             Size (nrows,ncols,3,simstep_max) matrix whose entries
%                   (i,j,k,simstep) are, for every fixed simstep, indicator
%                   variables of open edges between well (i,j) and its neighbor
%                   to the (k==1) right, (k==2) right down, (k==3) left down.
% edgelist          List of the edges (connected neighbors) from the first of
%                   simstep_max simulations (simstep==1). It has size (4, no.
%                   of open edges). Each column stores an edge (v_1, v_2) as
%                   [rowindex(v_1); colindex(v_1); rowindex(v_2); colindex(v_2)].
%
% lambda            Row vector of seeding rates for the ncolors colors
% mu                Contamination rate (probability of any given
%                   undirected edge being open)
% simstep_max       Total number of simulations
% simstep           Index of current simulation in 1:simstep_max.
%
% totalvertices     The total number of vertices in the lattice that belong to
%                   the area of interest. Defined as the sum of all entries of
%                   experiment_matrix(:,:,1), which is identical to max_vertices.
% threshvertices    Row vector of length ncolors, the numbers of seeds for
%                   individual colors that will be put into vertices(:,:,:,simstep).
% totaledgelocations The total number of edges, whether open or closed,
%                   whose both ends belong to the area of interest. Defined
%                   as the sum of all entries of max_edges.
% threshedges       Number of edges that will be put into edges(:,:,:,simstep).
%                   These are distributed randomly across edges(:,:,:,simstep)
%                   for every fixed simstep.
% rplabelled_colors Matrix of size (nrows, ncols, ncolors, simstep_max) whose
%                   entries are 0 outside the area of interest (i.e. where
%                   max_vertices==0) and a random permutation of (1, 2, ...,
%                   totalvertices) for every fixed i and simstep in
%                   rplabelled_colors(:,:,i,simstep).
% rplabelled_edges  Matrix of size (nrows, ncols, 3, simstep_max) whose entries
%                   are 0 where there can be no edge (i.e. where max_edges==0)
%                   and a random permutation of (1, 2, ..., totaledgelocations)
%                   for every fixed simstep in rplabelled_edges(:,:,:,simstep).
%
% components1{i}    Matrix which contains in one row the elements of the
%                   i'th connected component (linear indexing). Its purpose
%                   is that each vertex within a connected component will
%                   ultimately have the same color.
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 16/7/2015, 2/2/2017


vertices=cat(3, max_vertices, zeros(nrows,ncols,ncolors));
vertices=repmat(vertices, [1,1,1,simstep_max]);
edges=zeros(nrows,ncols,3,simstep_max);

threshvertices=round(lambda*totalvertices);
% We use the random permutations rplabelled_colors so that the proportions
% of seeded vertices are very accurately the respective parameters.

% We add a seed of color i if rplabelled_colors is both non-zero (greater
% than 0) and it is at most threshvertices(i).
vertices(:,:,2:1+ncolors,:) = (0<rplabelled_colors) .* (rplabelled_colors<=repmat(reshape(threshvertices,[1,1,ncolors,1]),[nrows,ncols,1,simstep_max]));

inivertices=vertices; % seeding state of vertices saved before contamination


threshedges=round(mu*totaledgelocations);
% We use the random permutations rplabelled_edges so that the proportion of
% open edges is very accurately the respective parameter.

% We add an edge if rplabelled_edges is both non-zero (greater than 0) and
% it is at most threshedges.
edges = (0<rplabelled_edges) .* (rplabelled_edges<=threshedges);

for simstep=1:simstep_max

% Find all the contaminated vertices

    if simstep==1 % For later visualization of the first simulation, we store edgelist.
        [edgelist, components1]=contamination(edges(:,:,:,simstep), nrows, ncols, originalrows, shiftedrows);
    else
        [~, components1]=contamination(edges(:,:,:,simstep), nrows, ncols, originalrows, shiftedrows);
    end

    if length(components1{1})>0 % Only do the following if there is at least one edge, one nontrivial component. Otherwise there is a size mismatch between (0,0) and (0,2) defining comp.

    % Transform vertices for each simulation into 3D for linear indexing
    vertices3d=vertices(:,:,:,simstep);

    for i=1:size(components1,2) % Loop for all connected components of contamination

% for j=1:ncolors, j*nrows*ncols is justified by the dimensions of the
% matrix vertices; this way it 'jumps' into the (:,:,2:ncolors+1) matrices
% which represent all the colors
%{
% More transparent, original version applied contamination for connected
% component i and color j:
        for j=1:ncolors
            comp=components1{i}+j*nrows*ncols;
            vertices3d(comp)=max(vertices3d(comp)); 
        end
%}

        % We pair elements of components1{i} with 1:ncolors in all possible ways.
        comp=repmat(components1{i}',[1, ncolors])+repmat(1:ncolors, [length(components1{i}), 1])*nrows*ncols;

        % Apply contamination for connected component i and color j for j=1:ncolors in one go
        vertices3d(comp)=repmat(max(vertices3d(comp),[],1),length(components1{i}),1);

    end
    vertices(:,:,:,simstep)=vertices3d; % Transform vertices back into 4D

    end % of if length(components1{1})>0

end % of for simstep

end
