function [av, avnb] = simcalcs(nrows, ncols, ncolors, originalrows, shiftedrows, experiment_matrix, totalvertices, totaledgelocations)

% SIMCALCS.m calculates the moments we consider in the parameter
% estimation: average densities of colors and of neighbors with the same
% color.
% (Multiple colors per vertex as an additional moment is not used
% and its implementation for fixed ncolors==3 is commented out.)
%
% EXPLANATION OF VARIABLES
%
% nrows             Number of rows
% ncols             Number of columns
% ncolors           Number of colors. E.g. RGB corresponds to ncolors==3.
% simstep_max       Number of simulations
% experiment_matrix Size: (nrows,ncols,1+ncolors,simstep_max). Each
%                   simulation step has one matrix of size
%                   (nrows,ncols,1+ncolors). (:,:,1) represents the
%                   area of interest. (:,:,2:ncolors+1) represent the vertices of
%                   the data set that can be filled with individual colors
%                   (e.g. if vertex (i,j) of the data set has color type 2,
%                   then synthetic_matrix(i,j,3)==1, else 0).
%                   See: synthetic_matrix in msm.m.
% av                Column vector of length ncolors. Entry i is the relative
%                   frequency of vertices with color i present in the dataset.
% avnb              Column vector of length ncolors. Entry i is the relative
%                   frequency of neighbor pairs with both having color i in the
%                   dataset.
%{
In this version not used:
% RGsim             Matrix with 1 where there is both R and G color, else 0
% RBsim             Matrix with 1 where there is both R and B color, else 0
% GBsim             Matrix with 1 where there is both G and B color, else 0
% avRG              Average #R&G per vertex (divided by total number of vertices & divided by number of simulations)
% avRB              Average #R&B per vertex (divided by total number of vertices & divided by number of simulations)
% avGB              Average #G&B per vertex (divided by total number of vertices & divided by number of simulations)
% avtwocol          [RG,RB,GB]-vector of average over all vertices of instances of multiple colors present in one vertex
%}
%
% totalvertices     The number of vertices in the lattice that belong to
%                   the area of interest.
% totaledgelocations The number of edges, whether open or closed, whose
%                   both ends belong to the area of interest.
% originalrows      Rows in experiment_matrix that are 'further left'. See
%                   explanation of shape
% shiftedrows       Rows in experiment_matrix that are 'further right'. See
%                   explanation of shape
% nb_r              vector of number of same colored neighbors in all
%                   simulations along right directed edge
% nb_dr             vector of number of same colored neighbors in all
%                   simulations along right down directed edge
% nb_dl             vector of number of same colored neighbors in all
%                   simulations along left down directed edge
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 16/7/2015, 3/2/2017

simstep_max=size(experiment_matrix,4); % By original definition of simstep_max

av=zeros(ncolors, simstep_max);
avnb=zeros(ncolors, simstep_max);

for i=1:ncolors

% The number of right directed edges with color i in both end vertices
    nb_r=sum(sum(experiment_matrix(:, 1:(ncols-1), 1+i, :) == experiment_matrix(:, 2:ncols, 1+i, :) & experiment_matrix(:, 1:(ncols-1), 1+i, :) == 1,1),2);

% The number of possible down right directed edges is one less if it starts from a right-shifted row
    nb_dr=sum(sum(experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:(ncols-1), 1+i, :) == experiment_matrix(shiftedrows(shiftedrows~=nrows)+1, 2:ncols, 1+i, :) ...
        & experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:(ncols-1), 1+i, :) == 1, 1),2)+...
        sum(sum(experiment_matrix(originalrows(originalrows~=nrows), 1:ncols, 1+i, :) == experiment_matrix(originalrows(originalrows~=nrows)+1, 1:ncols, 1+i, :) ...
        & experiment_matrix(originalrows(originalrows~=nrows), 1:ncols, 1+i, :) == 1, 1),2);

% The number of possible down left directed edges is one less if it starts from a not shifted row
    nb_dl=sum(sum(experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:ncols, 1+i, :) == experiment_matrix(shiftedrows(shiftedrows~=nrows)+1, 1:ncols, 1+i, :) ...
        & experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:ncols, 1+i, :) == 1, 1),2)+...
        sum(sum(experiment_matrix(originalrows(originalrows~=nrows), 2:ncols, 1+i, :) == experiment_matrix(originalrows(originalrows~=nrows)+1, 1:(ncols-1), 1+i, :) ...
        & experiment_matrix(originalrows(originalrows~=nrows), 2:ncols, 1+i, :) == 1, 1),2);

    av(i,:)=sum(sum(experiment_matrix(:,:,1+i,:),1),2)/totalvertices; % The right-hand side is size (1,1,1,simstep_max), is converted automatically.
    avnb(i,:)=(nb_r+nb_dr+nb_dl)/totaledgelocations; % The right-hand side is size (1,1,1,simstep_max), is converted automatically.

end % of for i

av=mean(av,2);
avnb=mean(avnb,2);

%{
Old implementation for 3 colors with no replacement (this capability has been lost and would need two for loops to go through all pairs of colours)
% Calculate number of vertices with (at least) two colors by adding
% together all vertices of all simulations that are RG,GB or RB. Divide by
% total number of vertices and by total number of simulations to get average
% per vertex.
RGsim=experiment_matrix(:,:,2,:)==experiment_matrix(:,:,3,:) & experiment_matrix(:,:,2,:)==1;
RBsim=experiment_matrix(:,:,2,:)==experiment_matrix(:,:,4,:) & experiment_matrix(:,:,2,:)==1;
GBsim=experiment_matrix(:,:,3,:)==experiment_matrix(:,:,4,:) & experiment_matrix(:,:,3,:)==1;

avRG=sum(sum(sum(RGsim(:,:,1,:)))./sum(sum(experiment_matrix(:,:,1,:))))/simstep_max;
avRB=sum(sum(sum(RBsim(:,:,1,:)))./sum(sum(experiment_matrix(:,:,1,:))))/simstep_max;
avGB=sum(sum(sum(GBsim(:,:,1,:)))./sum(sum(experiment_matrix(:,:,1,:))))/simstep_max;

avtwocol=[avRG,avRB,avGB];
%}
end