function [shiftedrows, originalrows, max_edges, experiment_matrix, totalwells, totaledgelocations] = createsynthdata_determ(ncolors, nrows, ncols, shape)

% CREATESYNTHDATA_DETERM.m creates the deterministic framework needed for the synthetic data: shiftedrows, originalrows, max_edges, totalwells [not needed for method==1], totaledgelocations [not needed for method==1]

% originalrows    Rows in wells that are 'further left'.
% shiftedrows     Rows in wells that are shifted right by a half unit.

% method==1 uses the Bernoulli variables,
% method==2 uses random permutations

% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 25/4/2016, 2-3/11/2016

% Create originalrows & shiftedrows
if shape==0;
    originalrows=2:2:nrows;
    shiftedrows=1:2:nrows;
else    %shape ==1
    originalrows=1:2:nrows;
    shiftedrows=2:2:nrows;
end


% We create three matrices max_edges(:,:,1:3) for possible open edges (ie
% possible contamination between two neighbors), needed for the calculations of the moments:
% Edge direction: 1=right, 2=right down, 3=left down. If an edge would
% connect a well outside and one inside the area of interest, it must be set to closed (0). We first allow every edge to be open randomly, then zero those which cannot be:
% (originalrows, 1): no left down; (shiftedrows, end): no right, no right down; (originalrows, end): no right; (end,:): no right down, no left down
max_edges=ones(nrows,ncols,3);
max_edges(originalrows,1,3)=0;
max_edges(shiftedrows,end,1)=0;
max_edges(shiftedrows,end,2)=0;
max_edges(originalrows,end,1)=0;
max_edges(end,:,[2 3])=0;

experiment_matrix=cat(3, ones(nrows,ncols,1), zeros(nrows,ncols,ncolors));
%experiment_matrix(:,:,1)=1;

% We comment out a part relevant only for irregular areas of interest from real data
%{
% Checking if the neighbors for each well are part of area of interest and
% if not, close the corresponding edge.
for i=1:nrows
    for j=1:ncols
        if experiment_matrix(i,j,1)==1 % 'Area of interest'
            if max(i==originalrows)==1 % If row is part of originalrows
                if j<ncols && experiment_matrix(i,j+1,1)==0 % If right neighbor is not in area of interest
                    max_edges(i,j,1)=0; % right edge cannot exist
                end
                if i<nrows && experiment_matrix(i+1,j,1)==0 % If right down neighbor is not in area of interest
                    max_edges(i,j,2)=0; % right down edge cannot exist
                end
                if i<nrows && j>1 && experiment_matrix(i+1,j-1,1)==0 % If down left neighbor is not in area of interest
                    max_edges(i,j,3)=0; % left down edge cannot exist
                end
            else  % Analogue if row is part of shiftedrows
                if j<ncols && experiment_matrix(i,j+1,1)==0
                    max_edges(i,j,1)=0;
                end
                if j<ncols && i<nrows && experiment_matrix(i+1,j+1,1)==0
                    max_edges(i,j,2)=0;
                end
                if i<nrows && experiment_matrix(i+1,j,1)==0
                    max_edges(i,j,3)=0;
                end
            end
        else
            max_edges(i,j,:)=0; % If (i,j) is not in area of interst, all three edges are zero
        end
    end
end
%}

% In fact, these are needed only for method==2
totalwells=sum(sum(experiment_matrix(:,:,1)));
totaledgelocations=sum(sum(sum(max_edges)));

%save(sprintf('dat%ix%i_highcont.mat',nrows,ncols), 'experiment_seeds', 'experiment_edges', 'experiment_matrix', 'experiment_edgelist', 'shape', 'lambda', 'mu', 'method');

%save(sprintf('data%i.mat',nrows), 'experiment_seeds', 'experiment_edges', 'experiment_matrix', 'experiment_edgelist', 'shape', 'lambda', 'mu', 'method');
end
