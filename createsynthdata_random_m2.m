function [rplabelled_colors, rplabelled_edges] = createsynthdata_random_m2(nrows, ncols, ncolors, max_vertices, max_edges, totalvertices, totaledgelocations, simstep_max)

% CREATESYNTHDATA_RANDOM_M2.m is creating the randomness according to
% Method 2 (random permutations) that will be kept fixed and used with
% different (lambda, mu) values to create realizations.
%
% max_vertices      Matrix of size (nrows,ncols), with 0 and 1 entries, which
%                   are the indicators of the area of interest. It's defined by
%                   experiment_matrix(:,:,1).
%                   If max_vertices(i,j,1)==1, then that vertex is going to be
%                   part of the calculations. If max_vertices(i,j,1)==0, then the
%                   corresponding matrix entries of colors will be ignored in any
%                   further calculations.
% j1, j2            Indices of the 1s in the (nrows, ncols) matrix max_vertices,
%                   i.e. where a seed is even possible.
% max_edges         Matrix of size (nrows,ncols,3) with 0 and 1 entries, which
%                   are the indicators of the edges whose both endpoints are in
%                   the area of interest.
% j                 Linear index of the 1s in the (nrows, ncols, 3) matrix
%                   max_edges, i.e. where there is an edge to be chosen open or
%                   closed.
% rplabelled_colors Matrix of size (nrows, ncols, ncolors, simstep_max) whose
%                   entries are 0 outside the area of interest (i.e. where
%                   max_vertices==0) and a random permutation of (1, 2, ...,
%                   totalvertices) for every fixed i and simstep in
%                   rplabelled_colors(:,:,i,simstep). It is the source of
%                   randomness for the seeding of vertices.
% rplabelled_edges  Matrix of size (nrows, ncols, 3, simstep_max) whose entries
%                   are 0 where there can be no edge (i.e. where max_edges==0)
%                   and a random permutation of (1, 2, ..., totaledgelocations)
%                   for every fixed simstep in rplabelled_edges(:,:,:,simstep).
%                   It is the source of randomness for the edges.
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 3/2/2017


% Preparing vertices

[j1, j2]=find(max_vertices); % indices of vertices where seeds are possible

rplabelled_colors=zeros(nrows, ncols, ncolors, simstep_max);

% j1(:), j2(:) enforce j1 and j2 to be column vectors; only needed when
% matrix consists of one row because otherwise j1, j2 are column vectors
% even without (:).
for simstep=1:simstep_max
    for i=1:ncolors
%{
        % More transparent, original version
        rpcolors(:,i,simstep)=randperm(totalvertices);
        rplabelled_colors(sub2ind([nrows, ncols, ncolors, simstep_max],j1,j2,i,simstep)) = rpcolors(:,i,simstep);
%}
        rplabelled_colors(sub2ind([nrows, ncols, ncolors, simstep_max], j1(:), j2(:), repmat(i, [totalvertices,1]), repmat(simstep, [totalvertices,1]))) = randperm(totalvertices);
    end
end


% Preparing edges

j=find(max_edges);

rplabelled_edges=zeros(nrows, ncols, 3, simstep_max);

for simstep=1:simstep_max
    rplabelled_edges(j+nrows*ncols*3*(simstep-1)) = randperm(totaledgelocations);
end

end
