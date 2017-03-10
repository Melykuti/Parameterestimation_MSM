function [ncolors, wells_before_contamination, wells_after_contamination, open_edges, experiment_edgelist, max_vertices, max_edges, originalrows, shiftedrows, true_values, totalvertices, totaledgelocations] = createsynthdata(nrows, ncols, simstep_max, shape, lambda, mu, method)
%function createsynthdata(nrows, ncols, simstep_max, shape, lambda, mu, method)
%
% CREATESYNTHDATA.m provides an easy way to create (and potentially save to
% file) a random dataset. It can be used to create a synthetic dataset for
% testing MSM.m or for visualization.
%
% Possible calls of this function
% 
% [ncolors, wells_before_contamination, wells_after_contamination, open_edges, experiment_edgelist, max_vertices, max_edges, originalrows, shiftedrows, true_values, totalvertices, totaledgelocations] = createsynthdata(10, 10, 1, 0, [0.1 0.05 0.2], 0.1, 1)
% For figures:
% createsynthdata(8, 12, 1, 0, [0.1 0.05 0.2], 0.1, 1)
% With critical mu:
% createsynthdata(20, 16, 1, 0, [0.05 0.1 0.2], 2*sin(pi/18), 1)
%
% to check Table 6:
% createsynthdata(10, 10, 1, 0, [0.02 0.04 0.05], 0.03, 1)
% createsynthdata(25, 25, 1, 0, [0.1 0.05 0.07], 0.06, 1)
% createsynthdata(100, 100, 1, 0, [0.07 0.05 0.04], 0.03, 1)
% createsynthdata(300, 300, 1, 0, [0.05 0.06 0.03], 0.02, 1)
% createsynthdata(500, 500, 1, 0, [0.03 0.04 0.05], 0.02, 1)
% createsynthdata(661, 661, 1, 0, [0.02 0.04 0.05], 0.03, 1)
% createsynthdata(800, 800, 1, 0, [0.02 0.04 0.05], 0.03, 1)
%
%
% INPUTS:
%    nrows, ncols, simstep_max, shape, lambda, mu, method
%
% OUTPUTS:
%    ncolors, wells_before_contamination, wells_after_contamination,
%    open_edges, experiment_edgelist, max_vertices, max_edges, originalrows,
%    shiftedrows, true_values, totalvertices, totaledgelocations
%
% VARIABLES SAVED (optionally):
%    The union of INPUTS and OUTPUTS.
%
% EXPLANATION OF VARIABLES
%
% nrows             Number of rows in the grid
% ncols             Number of columns in the grid
% simstep_max       Number of independent realizations of the grid
%
% shape             Shape of the grid:
%
% shape==0 -- Odd numbered rows are shifted to the right by half a unit
%  1 2 3 ... ncols-1 ncols (shiftedrows)
% 1 2 3 ... ncols-1 ncols  (originalrows)
%  1 2 3 ... ncols-1 ncols (shiftedrows)
% ...
%
% shape==1 -- Even numbered rows are shifted to the right by half a unit
% 1 2 3 ... ncols-1 ncols  (originalrows)
%  1 2 3 ... ncols-1 ncols (shiftedrows)
% 1 2 3 ... ncols-1 ncols  (originalrows)
% ...
%
% lambda            Row vector of seeding rates for the ncolors colors
% mu                Contamination rate (probability of any given
%                   undirected edge being open)
% method            method==1 uses Bernoulli variables,
%                   method==2 uses random permutations to generate random realizations.
% filename          String, name of file to save the data. If filename=''
%                   is set, then data is not saved.
%
% Outputs:
%
% ncolors           Number of colors, equals length(lambda).
% wells_before_contamination & inivertices
%                   Size: (nrows, ncols, 1+ncolors, simstep_max), essentially
%                   simstep_max independent random realizations of the seeding,
%                   before contamination would have taken place. Primary use
%                   is in visualization. See also vertices.
% wells_after_contamination & vertices
%                   Size: (nrows, ncols, 1+ncolors, simstep_max), essentially
%                   simstep_max independent random realizations of the dataset.
%                   vertices(:,:,1,simstep) are indicators of the area of interest. 
%                   (These are identical for simstep=1:simstep_max.)
%                   vertices(i,j,2:ncolors+1) represent the presence 
%                   of individual colors in the dataset in the given vertex (i,j).
% open_edges & edges
%                   Size (nrows,ncols,3,simstep_max) matrix whose entries
%                   (i,j,k,simstep) are, for every fixed simstep, indicator
%                   variables of open edges between well (i,j) and its neighbor
%                   to the (k==1) right, (k==2) right down, (k==3) left down.
% experiment_edgelist
%                   List of the edges (connected neighbors) from the first of
%                   simstep_max simulations (simstep==1). It has size (4, no.
%                   of open edges). Each column stores an edge (v_1, v_2) as
%                   [rowindex(v_1); colindex(v_1); rowindex(v_2); colindex(v_2)].
% max_vertices      Matrix of size (nrows,ncols), with 0 and 1 entries, which
%                   are the indicators of the area of interest.
%                   If max_vertices(i,j)==1, then that vertex is going to be
%                   part of the calculations. If max_vertices(i,j)==0, then
%                   the corresponding matrix entries of colors will be
%                   ignored in any further calculations.
% max_edges         Matrix of size (nrows,ncols,3) that indicates if the edge
%                   can possibly be open between vertices (i,j) and three
%                   of its neighbors. It is 1 when both endpoints are in
%                   the area of interest (i.e. max_vertices(.,.)==1 for
%                   both endpoints), otherwise 0. The last dimension
%                   represents the three possible directions of edges. Edge
%                   direction: 1==right, 2==right down, 3==left down.
%                   E.g. if there could be contamination between vertex
%                   (i,j) and its right down neighbor, max_edges(i,j,2)==1.
%                   If not (i.e. the neighbor lies outside the area of
%                   interest), then max_edges (i,j,2)==0.
% originalrows      Row vector of row numbers of vertices that are in
%                   default position, 'further left'.
%                   See explanation of shape.
% shiftedrows       Row vector of row numbers of vertices that are 'further
%                   right' (shifted right by a half unit relative to originalrows).
%                   See explanation of shape.
% true_values       Row vector [lamda,mu]
% totalvertices     The total number of vertices in the lattice that belong to
%                   the area of interest. Defined as the sum of all entries
%                   of max_vertices.
% totaledgelocations The total number of edges, whether open or closed,
%                   whose both ends belong to the area of interest. Defined
%                   as the sum of all entries of max_edges.
%
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 25/4/2016, 2/11/2016, 30/1-9/2/2017

filename='707x707_testdata_m1.mat';
%filename='';

%rng('default'); % Reset random generator

ncolors=length(lambda);
experiment_matrix=cat(3, ones(nrows,ncols), zeros(nrows,ncols,ncolors));
%{
% There is an opportunity here to trim the area of interest by zeroing out
% certain entries of experiment_matrix(:,:,1)
experiment_matrix(:,:,1)=0;
%}
max_vertices=experiment_matrix(:,:,1);
[originalrows, shiftedrows, max_edges, totalvertices, totaledgelocations] = createsynthdata_determ(nrows, ncols, ncolors, shape, max_vertices);

if method==1
    rvcolors=rand(nrows,ncols,ncolors,simstep_max);
    rvedges=rand(nrows,ncols,3,simstep_max);
    rplabelled_colors=0; % unused, dummy variable
    rplabelled_edges=0; % unused, dummy variable
else % i.e. method==2
    rvcolors=0; % unused, dummy variable
    rvedges=0; % unused, dummy variable
    [rplabelled_colors, rplabelled_edges] = createsynthdata_random_m2(nrows, ncols, ncolors, experiment_matrix(:,:,1), max_edges, totalvertices, totaledgelocations, simstep_max);
end

% Start simulation(s)
if method==1
    [inivertices, vertices, edges, experiment_edgelist] = simulation_m1(nrows, ncols, ncolors, shape, originalrows, shiftedrows, lambda, mu, simstep_max, max_vertices, max_edges, rvcolors, rvedges);
else % i.e. method==2
    [inivertices, vertices, edges, experiment_edgelist] = simulation_m2(nrows, ncols, ncolors, shape, originalrows, shiftedrows, lambda, mu, simstep_max, max_vertices, max_edges, totalvertices, totaledgelocations, rplabelled_colors, rplabelled_edges);
end

wells_before_contamination=inivertices;
wells_after_contamination=vertices;
open_edges=edges;
true_values=[reshape(lambda,[1 ncolors]), mu];

if length(filename)>0 % then save the data
    % This corresponds to Felix Beck's nomenclature.
    
    % save(sprintf('dat%ix%i_3colors_m2.mat', nrows, ncols), ...
    save(filename,...
        'nrows', 'ncols', 'ncolors', 'simstep_max', 'shape',...
        'wells_before_contamination', 'wells_after_contamination',...
        'open_edges', 'experiment_edgelist', 'max_vertices', 'max_edges',...
        'originalrows', 'shiftedrows', 'lambda', 'mu', 'method',...
        'true_values', 'totalvertices', 'totaledgelocations'); 
    
    % Alternative names
    % 'inivertices', 'vertices', 'edges'
end
end