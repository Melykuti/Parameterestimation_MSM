function [solutions, result, experiment_matrix] = msm(simstep_max, loops, mu_max, method)

% MSM.m is the main function for the parameter estimation by the method of
% simulated moments. It loads the individual colour channel matrices and
% the shape of the dataset from a file. The moments we need for the
% estimation are called from simcalcs.m. Initial guesses for the estimators
% are chosen, the objective function for the optimization is computed in
% optim.m and then optimized with fminsearchbnd. The best estimator is
% saved in a file.
%
% This works with any number of colors (i.e. color channels).
% It saves the following variables to a file: 'estimators', 'time',
% 'errors', 'input_values', 'optimizations'
%
%
% EXPLANATION OF VARIABLES
%
% simstep_max        Number of simulations for each parameter vector tested
% loops              Number of optimizations with different initial values
% mu_max             Maximum initial mu value for optimizations. Choose at
%                    most 0.2 as fminsearchbnd is set to look for mu <=0.2.
% method             method==1 uses Bernoulli variables,
%                    method==2 uses random permutations to generate random realizations.
%
% nrows              Number of rows in the grid
% ncols              Number of columns in the grid
% ncolors            Number of colors, equals length(lambda).
% shape              Shape of the grid:
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
% originalrows      Row vector of row numbers of vertices that are in
%                   default position, 'further left'.
%                   See explanation of shape.
% shiftedrows       Row vector of row numbers of vertices that are 'further
%                   right' (shifted right by a half unit relative to originalrows).
%                   See explanation of shape.
%                   Note:
%                   The main purpose is to know the relative positioning of
%                   odd and even rows. Those that are 'further right' are
%                   called shiftedrows. It is irrelevant whether shiftedrows
%                   were actually shifted to the right or if originalrows
%                   were shifted to the left. The set theoretic union of
%                   originalrows and shiftedrows is 1:nrows.
%
% lambda            Row vector of seeding rates for the ncolors colors
% mu                Contamination rate (probability of any given
%                   undirected edge being open)
% lambda_max        maximum possible lambda, chosen also as maximum initial lambda
%                   Coordinate i is defined by sum of all vertices of color i divided
%                   by total number of vertices. (Obviously, there could not have been
%                   more seeds than how many times the color is observed after
%                   contamination).
%
% experiment_matrix Size: (nrows, ncols, 1+ncolors), essentially 1+ncolors
%                   matrices of size (nrows,ncols) stacked with entries of
%                   either 0 or 1. They store the observations of the experimental
%                   the dataset on the grid - . experiment_matrix(:,:,1) represents
%                   area of interest. E.g. it can happen that some parts of a
%                   column or row are not within the limits of our original
%                   image. If experiment_matrix(i,j,1)==1, then that vertex
%                   is going to be part of the calculations. If
%                   experiment_matrix(i,j,1)==0, the corresponding
%                   matrix entries of colors will be ignored in any further
%                   calculations.
%                   experiment_matrix(i,j,2:ncolors+1) represent the presence 
%                   of individual colors in the dataset in the given vertex
%                   (e.g. if vertex (i,j) of the dataset has color type 2,
%                   then experiment_matrix(i,j,3)=1, else 0).
% max_edges         Matrix of size (nrows,ncols,3) that indicates if the edge
%                   can possibly be open between vertices (i,j) and three
%                   of its neighbors. It is 1 when both endpoints are in
%                   the area of interest (i.e. experiment_matrix(.,.,1)==1
%                   for both endpoints), otherwise 0. The last dimension
%                   represents the three possible directions of edges. Edge
%                   direction: 1==right, 2==right down, 3==left down.
%                   E.g. if there could be contamination between vertex
%                   (i,j) and its right down neighbor, max_edges(i,j,2)==1.
%                   If not (i.e. the neighbor lies outside the area of
%                   interest), then max_edges (i,j,2)==0.
% totalvertices     The total number of vertices in the lattice that belong to
%                   the area of interest. Defined as the sum of all entries of
%                   experiment_matrix(:,:,1).
% totaledgelocations The total number of edges, whether open or closed,
%                   whose both ends belong to the area of interest. Defined
%                   as the sum of all entries of max_edges.
%
% rvcolors          Independent uniform random variables on [0,1] of size
%                   (nrows,ncols,ncolors,simstep_max), used to define seeding
%                   in method 1.
% rvedges           Independent uniform random variables on [0,1] of size
%                   (nrows,ncols,3,simstep_max), used to define which edges
%                   are open in method 1.
% rplabelled_colors (nrows, ncols, ncolors, simstep_max)-matrix whose entries are 0
%                   outside the area of interest (i.e. where experiment_matrix(:,:,1)==0)
%                   and a random permutation of (1, 2, ..., totalvertices) for every
%                   fixed i and simstep in rplabelled_colors(:,:,i,simstep). It is used
%                   to define seeding in method 2.
% rplabelled_edges  (nrows,ncols,3,simstep_max)-matrix whose entries are 0 where
%                   there can be no edge (i.e. where max_edges==0) and a random
%                   permutation of (1, 2, ..., totaledgelocations) for every fixed
%                   simstep in rplabelled_edges(:,:,:,simstep). It is used to
%                   define which edges are open in method 2.
%
% dsRGB             Column vector of length ncolors. Entry i is the relative
%                   frequency of vertices with color i present in the experimental
%                   dataset.
% dsnb              Column vector of length ncolors. Entry i is the relative
%                   frequency of neighbor pairs with both having color i in the
%                   experimental dataset.
% dstwocol          Currently unused. Each of its entries would mean relative
%                   frequency of cooccurrence of two colors in a vertex in the
%                   experimental dataset.
%
% solutions         Matrix of size (loops, 2*(ncolors+1)+2). For every row l, it contains
%                   the initial [lambda, mu], the corresponding squared error between
%                   moments of the experimental dataset and of the synthetic data (the
%                   objective function value), the optimal [lambda, mu] for this initial
%                   state, and the corresponding squared error.
% optimizations     Cell array of size (1+loops, 1+2*(ncolors+1)+2). It contains the
%                   matrix solutions prepended with the numbers 1:loops in the
%                   first column and with labels in the first row.
% estimator         Row vector, the optimal [lambda, mu] values from the row of
%                   solutions which has the lowest squared error in solutions(:,end).
%                   This is the solution of the problem.
% result            Cell array of size (2,1) where {1} is a label, {2} is the
%                   row vector estimator.
% estimators        Cell array of size (2,2) where {1,1} and {2,1} are the labels
%                   'lambda' and 'mu', respectively. {2,1} and {2,2} are the optimal
%                   lambda and mu, respectively.
% time              Cell array of size (1,2) where {1} is a label, {2} is the
%                   running time of the parameter estimation.
% errors            Cell array of size (2,2) where {1,1} and {2,1} are labels,
%                   {1,2} and {2,2} are squared errors between moments of the
%                   experimental dataset and the synthetic data, ({1,2}) with the
%                   trivial [lambda, mu] with mu==0 and lambda the relative
%                   frequencies of colors (i.e. sum of all colored vertices / total
%                   number of vertices) and ({2,2}) with the optimal [lambda, mu].
% input_values      Cell array of size (4,2) where the first column is labels,
%                   the second column stores the inputs of MSM.m.
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 16/7/2015, 3/2/2017


% Load the color matrices and shape of the dataset whose parameters we
% estimate, i.e. the variables 'experiment_matrix' and 'shape'.
filename=input('Filename for dataset: ','s');
load(filename);

dots=strfind(filename,'.'); % Find the dot in the filename
if length(dots)>0
    outputsuggestion=sprintf('%s_estimate',filename(1:max(dots)-1));
else
    outputsuggestion=sprintf('%s_estimate',filename);
end
outputname=input(sprintf('Enter filename for output (or press Return for the default %s.mat): ',outputsuggestion),'s');

% This renaming needed for our benchmark synthetic data files:
experiment_matrix=wells_after_contamination;

tic % Starting stopwatch

% Setting sizes of simulation matrices
nrows=size(experiment_matrix,1);
ncols=size(experiment_matrix,2);
ncolors=size(experiment_matrix,3)-1;

[originalrows, shiftedrows, max_edges, totalvertices, totaledgelocations] = createsynthdata_determ(nrows, ncols, ncolors, shape, experiment_matrix(:,:,1));

% rng('shuffle'); % Shuffle random generator
% rng('default'); % Reset random generator to default initial state.

% Generating the randomness which is kept fixed throughout and used for all
% parameter values considered in the optimization
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

% Calculations of the moments for experimental dataset
% In earlier version: [dsRGB, dstwocol, dsnb] =
[dsRGB, dsnb] = simcalcs(nrows, ncols, ncolors, originalrows, shiftedrows, experiment_matrix, totalvertices, totaledgelocations);
dstwocol=0; % dstwocol is neither used nor implemented in simcalcs, we specify a dummy zero.

solutions=zeros(loops,2*(ncolors+1)+2); % This will store the solutions from every initial (lambda, mu) value

% Maximal initial lambda: for each color, sum of all colored vertices / total number of vertices
lambda_max=zeros(1,ncolors);
lambda_max(:)=sum(sum(experiment_matrix(:,:,2:end),1),2)/totalvertices;

for l=1:loops % For number of optimizations

    % Settings for initial estimator(s) of lambda and mu
    if l==1 % If first optimization
        lambda=lambda_max; % lambda gets the largest value possible
        mu=0; % mu gets the smallest value possible
    else            % If loops>1, mu is going up, lambda is going down; lambda is always >0
        lambda=lambda-lambda_max/loops; % New initial lambda estimators (assumption: lambda_max/loops<=lambda<=lambda_max)
        mu=mu+mu_max/(loops-1);  % New initial mu estimators (assumption: 0<=mu<=mu_max)
    end

    solutions(l,1:ncolors+1)=[lambda, mu]; % Estimators before optimization

    % Calculating the max error with initial values of lambda and mu to test if
    % the optimized estimator really gives smaller error

    solutions(l,ncolors+2) = optim([lambda,mu], nrows, ncols, ncolors, shape, originalrows, shiftedrows, simstep_max, dsRGB, dstwocol, dsnb, experiment_matrix(:,:,1), max_edges, totalvertices, totaledgelocations, rvcolors, rvedges, rplabelled_colors, rplabelled_edges, method);

    % fminsearchbnd behaves similarly to fminsearch, except you can add bound
    % constraints.
    % source and description: http://de.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon/content/FMINSEARCHBND/fminsearchbnd.m
    % author: John D'Errico

    % Options for fminsearchbnd optimization routine (function optim returns
    % the scalar squared error for lambda and mu)
    options = optimset('MaxFunEvals',20000,'MaxIter',20000,'TolFun',1e-7,'TolX',1e-7);
    % We bound search interval for mu from above by 0.2.
    [param,fval]=fminsearchbnd(@(param) optim(param, nrows, ncols, ncolors, shape, originalrows, shiftedrows, simstep_max, dsRGB, dstwocol, dsnb, experiment_matrix(:,:,1), max_edges, totalvertices, totaledgelocations, rvcolors, rvedges, rplabelled_colors, rplabelled_edges, method), [lambda,mu], zeros(1,ncolors+1), [lambda_max, 0.2], options);

    solutions(l,ncolors+3:2*(ncolors+1)+2)=[param,fval]; % Estimators after optimization and objective function value

end     % of for l

[~, rowmin]=min(solutions(:,end)); % Searching for best estimator (lowest fval from optim) in loops number of independent solutions

% Preparing output
estimator=solutions(rowmin(1), ncolors+3:2*(ncolors+1)+1); % Take first entry of rowmin in case there are multiple equally good solutions
result=[{'estimators'}; estimator]; % Output of the best estimator [lambda_1, ..., lambda_ncolors, mu]
absd=solutions(rowmin(1), end);

%{
To be deleted.
% Returning an optimized simulation as an output with visualization, using optim_output because the optimization function optim has limited output options
% (Note that if simstep_max>1, the image that will be shown as an output is not representative for the calculated errors/distances because the image only shows the first simulation!)
lambda=estimator(1:ncolors);
mu=estimator(1+ncolors);
[wells, absd_repeat]=optim_output([lambda,mu], nrows, ncols, ncolors, shape, originalrows, shiftedrows, simstep_max, dsRGB, dstwocol, dsnb, synthetic_matrix, max_edges, totalvertices, totaledgelocations, rvcolors, rvedges, rplabelled_colors, rplabelled_edges, method);
%}

optimizations=cell(size(solutions,1)+1,size(solutions,2)+1); % Create variable with all estimators and errors before and after optimizations
titles=['optim number'; 'lambda_ini  '];
for i=2:ncolors
   titles=[titles;'            ']; 
end
titles=[titles; 'mu_ini      '; 'error_ini   '; 'lambda      '];
for i=2:ncolors
   titles=[titles;'            ']; 
end
titles=[titles; 'mu          '; 'error       '];

optimizations(1,:)=cellstr(titles);
for i=1:size(solutions,1)
    optimizations(i+1,:)=num2cell([i,solutions(i,:)]);
end

% Display some variables for control purposes
optimizations
absd
result{2}     % Output of final estimators


% Create a file with parameters of the function call and results of the parameter estimation

% Create output for the best estimator
estimators=cell(2,2);
estimators{1,1}='lambda';
estimators{2,1}='mu';
estimators{1,2}=estimator(1:ncolors);
estimators{2,2}=estimator(ncolors+1);

% Create output for elapsed time
eltime=toc;
time=cell(1,2);
time_titles='Elapsed time in sec';
time(1)=cellstr(time_titles);
%time(2)=num2cell(eltime);
time{2}=eltime

% Create output for trivial error and optimized error
errors{1,1}='Trivial error';
errors{2,1}='Optimized error';
errors{1,2}=solutions(1,ncolors+2);
errors{2,2}=absd;

% Create output for user inputs
input_values=cell(4,2);
input_values{1,1}='Number of simulations n_s';
input_values{2,1}='Number of optimizations n_opt';
input_values{3,1}='Max initial mu value for optimizations mu_max';
input_values{4,1}='Method';
input_values(:,2)=num2cell([simstep_max; loops; mu_max; method]);

% User dialog whether to save the data into file

%choice='Yes'; % Uncomment this and comment out the two menu calls to shortcut to saving the data without prompt.

if length(outputname)==0 % If user decides to use the default name outputsuggestion

     % Check if the user would like to save the file
     choice = questdlg(['Do you want to save the results in the file ', outputsuggestion, '.mat?'],'Save output?','Yes','No','Yes');
     switch choice
         case 'Yes'
             save (sprintf('%s.mat',outputsuggestion), 'estimators', 'time', 'errors', 'input_values', 'optimizations');
             disp(['The file ', outputsuggestion, '.mat was created successfully.'])
         case 'No'
             disp('No file was created.')
     end

else % If user entered a file name

    % Check if the user would like to save the file
    choice = questdlg(['Do you want to save the results in the file ', outputname, '?'],'Save output?','Yes','No','Yes');
    switch choice
        case 'Yes'
            save (sprintf('%s',outputname), 'estimators', 'time', 'errors', 'input_values', 'optimizations');
            disp(['The file ', outputname, ' was created successfully.'])
        case 'No'
            disp('No file was created.')
    end

end % end of if length(outputname)==0

end
