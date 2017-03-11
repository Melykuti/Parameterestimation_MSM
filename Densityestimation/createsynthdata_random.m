function [experiment_seeds, experiment_edges, experiment_matrix, experiment_edgelist] = createsynthdata_random(nrows, ncols, shape, lambda, mu, method, simstep_max, shiftedrows, originalrows, max_edges, experiment_matrix_prototype, totalwells, totaledgelocations)

% CREATESYNTHDATA_RANDOM.m is using the deterministic framework created by createsynthdata_determ.m (shiftedrows, originalrows, max_edges, totalwells [not needed for method==1], totaledgelocations [not needed for method==1]) to generate the randomness and create the synthetic data.

% method==1 uses the Bernoulli variables,
% method==2 uses random permutations

% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 25/4/2016, 2-3/11/2016

ncolors=length(lambda);

if method==1
    rvRGB=rand(nrows,ncols,ncolors,simstep_max);
    rvedges=rand(nrows,ncols,3,simstep_max);

    [experiment_seeds, experiment_edges, experiment_matrix, experiment_edgelist] = simulation1(ncolors, nrows, ncols, shape, lambda, mu, simstep_max, rvRGB, experiment_matrix_prototype, rvedges, max_edges, shiftedrows, originalrows);

else % ie. method==2

% Preparing wells

indexed_wells=experiment_matrix_prototype(:,:,1);
counter=0;
for i=1:nrows*ncols
    counter=counter+indexed_wells(i);
    indexed_wells(i)=indexed_wells(i)*counter; % Only insert counter if not 0 (ie. 1, it belongs to area of interest)
end

rpwells=zeros(totalwells,ncolors,simstep_max); % This is the source of randomness for wells.
for simstep=1:simstep_max
    for i=1:ncolors
        rpwells(:,i,simstep)=randperm(totalwells);
    end
end

% Preparing edges

indexed_edges=max_edges;
counter=0;
for i=1:nrows*ncols*3
    counter=counter+max_edges(i);
    indexed_edges(i)=max_edges(i)*counter; % Only insert counter if not 0 (ie. 1)
end

rpedges=zeros(totaledgelocations,simstep_max); % This is the source of randomness for edges.
for simstep=1:simstep_max
    rpedges(:,simstep)=randperm(totaledgelocations);
end

[experiment_seeds, experiment_edges, experiment_matrix, experiment_edgelist] = simulation2(ncolors, nrows, ncols, shape, lambda, mu, simstep_max, experiment_matrix_prototype, indexed_wells, totalwells,  rpwells, max_edges, indexed_edges, totaledgelocations, rpedges, shiftedrows, originalrows);

end
%save(sprintf('dat%ix%i_highcont.mat',nrows,ncols), 'experiment_seeds', 'experiment_edges', 'experiment_matrix', 'experiment_edgelist', 'shape', 'lambda', 'mu', 'method');

%save(sprintf('data%i.mat',nrows), 'experiment_seeds', 'experiment_edges', 'experiment_matrix', 'experiment_edgelist', 'shape', 'lambda', 'mu', 'method');
end
