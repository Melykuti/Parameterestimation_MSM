function densityestimation(total_nrows, total_ncols, nrows, ncols, shape, simstep_max)
% DENSITYESTIMATION.m computes the densities of coloured vertices and of
% adjacent pairs of coloured vertices as a function of lambda (seeding
% rate) and mu (percolation parameter) when number of colours is 1. For
% each parameter vector (lambda, mu), we create a total_nrows x total_ncols
% lattice, and cut out the middle nrows x ncols sublattice to compute the
% densities.
%
% simstep_max       Number of simulations per parameter vector
%
% There are 4 lines to comment out or uncomment to toggle between saving
% data and not saving. These four are offset by one space character.
%
% createsynthdata(300, 300, 0, [0.1 0.02 0.05], 0.1, 1)
% densityestimation(300, 300, 100, 100, 0, 1)
% densityestimation(100, 100, 33, 33, 0, 1)
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 16/7/2015, 2-10/11/2016

tic

data_file_name='densities_run_1.mat';
mu_crit=2*sin(pi/18);
method=1;

% This is the part that ensures looping through a grid of (lambda, mu) parameter pairs
lambda=0; % This is here just to ensure length(lambda)==1, needed in the definition of ncolors.
%lambda_vector=linspace(0.001,0.9,5);
% In dataset_A.mat:
lambda_vector=logspace(log10(5e-4),0,11);
lambda_vector(end)=[]; % We know what happens for lambda==1.
% In dataset_B.mat:
%lambda_vector=logspace(log10(5e-4*(1/5e-4)^-0.6),log10(5e-4*(1/5e-4)^0.1),8);
mu_vector=logspace(-4,log10(0.5),9);
mu_vector=[mu_vector(1:8) mu_crit mu_vector(9)];

[lambda_array mu_array]=meshgrid(lambda_vector, mu_vector);
[nmus nlambdas]=size(lambda_array);
nparams=nmus*nlambdas;
% To go as for i in mu, for j in lambda, and not the other way around so that computation goes through lowest mu first
lambda_array=lambda_array'; mu_array=mu_array';

ncolors=length(lambda);

% Data storage, initialisation
lambda_st=zeros(0, ncolors);
mu_st=zeros(0,1);
total_nrows_st=zeros(0,1);
total_ncols_st=zeros(0,1);
nrows_st=zeros(0,1);
ncols_st=zeros(0,1);
shape_st=zeros(0,1);
simstep_max_st=zeros(0,1);
av_st=zeros(0, ncolors);
avnb_st=zeros(0, ncolors);
avseeds_st=zeros(0, ncolors);
var_st=zeros(0, ncolors);
varnb_st=zeros(0, ncolors);
varseeds_st=zeros(0, ncolors);
 save(data_file_name, '-v7', 'lambda_st', 'mu_st', 'total_nrows_st', 'total_ncols_st', 'nrows_st', 'ncols_st', 'shape_st', 'simstep_max_st', 'av_st', 'avnb_st', 'avseeds_st', 'var_st', 'varnb_st', 'varseeds_st');
% End of data storage, initialisation

for param=1:nparams
%for param=nparams:-1:41
%for param=49:56

lambda=lambda_array(param)
mu=mu_array(param)

av=zeros(simstep_max, ncolors);
avnb=zeros(simstep_max, ncolors);
avseeds=zeros(simstep_max, ncolors);

[shiftedrows, originalrows, max_edges, experiment_matrix_prototype, totalwells, totaledgelocations] = createsynthdata_determ(ncolors, total_nrows, total_ncols, shape);

shiftedrows_orig=shiftedrows;
originalrows_orig=originalrows;

row_1=floor((total_nrows-nrows)/2);
col_1=floor((total_ncols-ncols)/2);

% originalrows and shiftedrows for the central part
originalrows=originalrows(originalrows>=row_1+1 & originalrows<=row_1+nrows)-row_1;
shiftedrows=shiftedrows(shiftedrows>=row_1+1 & shiftedrows<=row_1+nrows)-row_1;


[experiment_seeds, experiment_edges, experiment_matrix, experiment_edgelist] = createsynthdata_random(total_nrows, total_ncols, shape, lambda, mu, method, simstep_max, shiftedrows_orig, originalrows_orig, max_edges, experiment_matrix_prototype, totalwells, totaledgelocations);
%experiment_seeds
%experiment_matrix

% We crop the experiment_matrix and use the central nrows x ncols submatrix.
experiment_matrix=experiment_matrix(row_1+1:row_1+nrows, col_1+1:col_1+ncols, :, :);
experiment_seeds =experiment_seeds(row_1+1:row_1+nrows, col_1+1:col_1+ncols, :, :);

% Count how many edges the central submatrix can have (totaledgelocations)
max_edges=ones(nrows,ncols,3);
max_edges(originalrows,1,3)=0;
max_edges(shiftedrows,end,1)=0;
max_edges(shiftedrows,end,2)=0;
max_edges(originalrows,end,1)=0;
max_edges(end,:,[2 3])=0;
totaledgelocations=sum(sum(sum(max_edges)));


for i=1:ncolors

 nb_r=sum(sum(experiment_matrix(:, 1:(ncols-1), 1+i, :) == experiment_matrix(:, 2:ncols, 1+i, :) & experiment_matrix(:, 1:(ncols-1), 1+i, :) == 1,1),2);

% the number of possible down right arrows is one less if it starts from a right-shifted row
 nb_dr=sum(sum(experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:(ncols-1), 1+i, :) == experiment_matrix(shiftedrows(shiftedrows~=nrows)+1, 2:ncols, 1+i, :) ...
        & experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:(ncols-1), 1+i, :) == 1, 1),2)+...
        sum(sum(experiment_matrix(originalrows(originalrows~=nrows), 1:ncols, 1+i, :) == experiment_matrix(originalrows(originalrows~=nrows)+1, 1:ncols, 1+i, :) ...
        & experiment_matrix(originalrows(originalrows~=nrows), 1:ncols, 1+i, :) == 1, 1),2);

% the number of possible down left arrows is one less if it starts from a not shifted row
 nb_dl=sum(sum(experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:ncols, 1+i, :) == experiment_matrix(shiftedrows(shiftedrows~=nrows)+1, 1:ncols, 1+i, :) ...
        & experiment_matrix(shiftedrows(shiftedrows~=nrows), 1:ncols, 1+i, :) == 1, 1),2)+...
        sum(sum(experiment_matrix(originalrows(originalrows~=nrows), 2:ncols, 1+i, :) == experiment_matrix(originalrows(originalrows~=nrows)+1, 1:(ncols-1), 1+i, :) ...
        & experiment_matrix(originalrows(originalrows~=nrows), 2:ncols, 1+i, :) == 1, 1),2);

% sum(sum(sum(experiment_matrix(:,:,1+i,:))))
% sum(nb_r+nb_dr+nb_dl)

%for simstep=1:simstep_max
 av(:,i)=sum(sum(experiment_matrix(:,:,1+i,:),1),2)/(ncols*nrows);
 avnb(:,i)=(nb_r+nb_dr+nb_dl)/totaledgelocations;
 avseeds(:,i)=sum(sum(experiment_seeds(:,:,1+i,:),1),2)/(ncols*nrows);
%end

end % of for i

vary=var(av,1)
varnb=var(avnb,1)
varseeds=var(avseeds,1)

av=mean(av,1);
avnb=mean(avnb,1);
avseeds=mean(avseeds,1);

 load(data_file_name);
lambda_st=[lambda_st; lambda];
mu_st=[mu_st; mu];
total_nrows_st=[total_nrows_st; total_nrows];
total_ncols_st=[total_ncols_st; total_ncols];
nrows_st=[nrows_st; nrows];
ncols_st=[ncols_st; ncols];
shape_st=[shape_st; shape];
simstep_max_st=[simstep_max_st; simstep_max];
av_st=[av_st; av];
avnb_st=[avnb_st; avnb];
avseeds_st=[avseeds_st; avseeds];
var_st=[var_st; vary];
varnb_st=[varnb_st; varnb];
varseeds_st=[varseeds_st; varseeds];
 save(data_file_name, '-v7', 'lambda_st', 'mu_st', 'total_nrows_st', 'total_ncols_st', 'nrows_st', 'ncols_st', 'shape_st', 'simstep_max_st', 'av_st', 'avnb_st', 'avseeds_st', 'var_st', 'varnb_st', 'varseeds_st');

disp(sprintf('Done with %i out of %i parameter combinations.',param,nparams))
end % of for param

running_time_sec=toc
 save(data_file_name, '-append', 'running_time_sec', 'nmus', 'nlambdas');
