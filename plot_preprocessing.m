function [nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(experiment_matrix, shape)

% PLOT_PREPROCESSING.m extracts important parameters from the dataset.
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 26/4/2016, 7/2/2017

nrows=size(experiment_matrix,1);
ncols=size(experiment_matrix,2);
ncolors=size(experiment_matrix,3)-1;

if shape==0;
    originalrows=2:2:nrows;
    shiftedrows=1:2:nrows;
else % shape==1
    originalrows=1:2:nrows;
    shiftedrows=2:2:nrows;
end

end
