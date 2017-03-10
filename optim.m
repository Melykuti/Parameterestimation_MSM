function  epsilon_new = optim(param, nrows, ncols, ncolors, shape, originalrows, shiftedrows, simstep_max, dsRGB, dstwocol, dsnb, max_vertices, max_edges, totalvertices, totaledgelocations, rvcolors, rvedges, rplabelled_colors, rplabelled_edges, method)

% OPTIM.m computes the objective function with given values of lambda and
% mu for the optimization of parameter values. Itself does not do
% optimization. The matrices for the simulation (i.e. vertices) are
% received from SIMULATION_M1.m or SIMULATION_M2.m and the moments needed
% for the optimization from SIMCALCS.m.
%
% EXPLANATION OF VARIABLES
%
% param             Row vector [lamda,mu]
% simstep_max       Number of simulations for a fixed parameter 
% vertices          Size: (nrows, ncols, 1+ncolors, simstep_max), essentially
%                   simstep_max independent random realizations of the dataset.
%                   They store the synthetic dataset, similarly to how
%                   experiment_matrix stores the experimental data.
%                   vertices(:,:,1,simstep) are indicators of the area of interest. 
%                   (These are identical for simstep=1:simstep_max.)
%                   vertices(i,j,2:ncolors+1) represent the presence 
%                   of individual colors in the dataset in the given vertex (i,j).
% max_vertices      Matrix of size (nrows,ncols), with 0 and 1 entries, which
%                   are the indicators of the area of interest. It's defined by
%                   experiment_matrix(:,:,1).
%                   If max_vertices(i,j)==1, then that vertex is going to be
%                   part of the calculations. If max_vertices(i,j)==0, then the
%                   corresponding matrix entries of colors will be ignored in any
%                   further calculations.
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
% dsRGB             Column vector of length ncolors. Entry i is the relative
%                   frequency of vertices with color i present in the experimental
%                   dataset.
% dsnb              Column vector of length ncolors. Entry i is the relative
%                   frequency of neighbor pairs with both having color i in the
%                   experimental dataset.
% dstwocol          Currently unused. Each of its entries would mean relative
%                   frequency of cooccurrence of two colors in a vertex in the
%                   experimental dataset.
% avRGB             Column vector of length ncolors. Entry i is the relative
%                   frequency of vertices with color i present in the synthetic
%                   dataset averaged over simstep_max realizations.
% avnb              Column vector of length ncolors. Entry i is the relative
%                   frequency of neighbor pairs with both having color i in the
%                   synthetic dataset averaged over simstep_max realizations.
% avtwocol          Currently unused. Each of its entries would mean relative
%                   frequency of cooccurrence of two colors in a vertex in the
%                   synthetic dataset averaged over simstep_max realizations.
% d_RGB             Column vector of normalized squared distances of averages
%                   dsRGB and avRGB.
% d_nb              Column vector of normalized squared distances of averages
%                   dsnb and avnb.
% d_twocol          Currently unused. It would be the column vector of
%                   normalized squared distances of averages dstwocol and avtwocol.
%
% absd              Column vector consisting of d_RGB and d_nb (and potentially d_twocol).
% epsilon_new       Sum of all distances in absd, which is to be minimized.
%
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 16/7/2015, 2/2/2017

lambda=param(1:ncolors);
mu=param(end);

% The simulation
if method==1

    [~, vertices, ~, ~] = simulation_m1(nrows, ncols, ncolors, shape, originalrows, shiftedrows, lambda, mu, simstep_max, max_vertices, max_edges, rvcolors, rvedges);

else % i.e. method==2

    [~, vertices, ~, ~] = simulation_m2(nrows, ncols, ncolors, shape, originalrows, shiftedrows, lambda, mu, simstep_max, max_vertices, max_edges, totalvertices, totaledgelocations, rplabelled_colors, rplabelled_edges);

end

% Get averages for simulations from simcalcs
[avRGB, avnb] = simcalcs(nrows, ncols, ncolors, originalrows, shiftedrows, vertices, totalvertices, totaledgelocations);
%avtwocol=0; % This is neither used nor implemented in simcalcs, we specify a dummy zero.

% Presetting for absolute distances between data set and simulations
d_RGB=zeros(ncolors,1);
d_nb=zeros(ncolors,1);
%d_twocol=0;

% Squared distances between data set and simulations normalised by moments of experimental dataset.
for i=1:ncolors
    if dsRGB(i)~=0
        d_RGB(i)=((dsRGB(i)-avRGB(i))/dsRGB(i))^2;
    else
        d_RGB(i)=(dsRGB(i)-avRGB(i))^2;
    end
    if dsnb(i)~=0
        d_nb(i)=((dsnb(i)-avnb(i))/dsnb(i))^2;
    else
        d_nb(i)=(dsnb(i)-avnb(i))^2;
    end
end

%{
Cooccurrence of two colors in same vertex would be a for loop over (ncolors choose 2):
    if dstwocol(i)~=0
        d_twocol(i)=0*((dstwocol(i)-avtwocol(i))/dstwocol(i))^2;
    else
        d_twocol(i)=0*(dstwocol(i)-avtwocol(i))^2;
    end
%}

absd=[d_RGB; d_nb];
epsilon_new=sum(absd); % the output is not absd but this scalar

end
