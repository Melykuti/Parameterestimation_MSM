function table_text=latex_table5

% LATEX_TABLE5.m turns saved data into LaTeX code for table 5 in the paper.
%
% In the resulting LaTeX code:
% 1) Replace every occurrence of e- and e+ in exponential notations with \mathrm{e-} and \mathrm{e+}.
%
% We could drop the column of n_p if the table is too wide.
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 7-13/2/2017


header=latex_table_header;
footer=latex_table_footer;
block='';
for n=1:4
    blocktext=latex_block(n);
    block=[block blocktext];
end
table_text=[header block footer];


    fileID=fopen('Text/table.tex','w');
    fprintf(fileID,table_text);
    fclose(fileID);
%}
end


function blocktext=latex_block(n)
nblockrows = 3*(n==1) + 2*(n==2) + 4*(n==3) + 2*(n==4); % how many rows each block consists of
blocktext='\\hline\n'; % the \hline is at top of each block, and not part of header
for k=1:nblockrows
    blockrow=latex_blockrow(n, nblockrows, k);
    blocktext=[blocktext blockrow];
end
end

function blockrow=latex_blockrow(n, nblockrows, k)
% k is the row index within block n
filename0=latex_finddataset(n);
load(filename0);
[nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(wells_after_contamination, shape);
[originalrows, shiftedrows, max_edges, totalvertices, totaledgelocations] = createsynthdata_determ(nrows, ncols, ncolors, shape, wells_after_contamination(:,:,1));
r=sprintf('%i\\\\times%i & %i & ', nrows, ncols, totaledgelocations);

filename1=latex_findestimate(n, k, 1);
filename2=latex_findestimate(n, k, 2);
load(filename1);

%r=[r num2str(input_values{1,2}) ' & ' num2str(input_values{2,2}) ' & ' num2str(input_values{3,2}) ' & '];
r=[r num2str(input_values{1,2}) ' & ' num2str(input_values{2,2}) ' & '];
simstep_max=input_values{1,2};
[dsRGB, dsnb] = simcalcs(nrows, ncols, ncolors, originalrows, shiftedrows, wells_after_contamination, totalvertices, totaledgelocations);
dstwocol=0;
alpha_triv(1) = errors{1,2};
alpha(1) = errors{2,2};

for method=1:2
rng('default'); % Reset random generator to default initial state.
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
    [rplabelled_colors, rplabelled_edges] = createsynthdata_random_m2(nrows, ncols, ncolors, wells_before_contamination(:,:,1), max_edges, totalvertices, totaledgelocations, simstep_max);
end

alpha_0(method) = optim(true_values, nrows, ncols, ncolors, shape, originalrows, shiftedrows, simstep_max, dsRGB, dstwocol, dsnb, wells_before_contamination(:,:,1), max_edges, totalvertices, totaledgelocations, rvcolors, rvedges, rplabelled_colors, rplabelled_edges, method);
%alpha_triv(method) = optim([lambda_max, 0], nrows, ncols, ncolors, shape, originalrows, shiftedrows, simstep_max, dsRGB, dstwocol, dsnb, wells_before_contamination(:,:,1), max_edges, totalvertices, totaledgelocations, rvcolors, rvedges, rplabelled_colors, rplabelled_edges, method);
end

load(filename2);
alpha_triv(2) = errors{1,2};
alpha(2) = errors{2,2};
r=[r sprintf('%.3g & %.3g & %.3g & %.3g & %.3g & %.3g', alpha_0(1), alpha_0(2), alpha_triv(1), alpha_triv(2), alpha(1), alpha(2))];

if n==4 & k==nblockrows
    r=[r '\n'];
else
    r=[r '\\\\\n'];
end
blockrow=r;
end


function header=latex_table_header
% Header for Table 5
    header='\\begin{table}\n\\begin{center}\n\\begin{tabular}{>{$}c<{$} >{$}c<{$} >{$}c<{$} >{$}c<{$} | >{$}c<{$}   >{$}c<{$}  >{$}c<{$}  >{$}c<{$}  >{$}c<{$}  >{$}c<{$}} \n n_I & n_p & n_s & n_{\\textrm{opt}} & \\alpha_{\\theta_0}^{(\\textrm{M}1)} & \\alpha_{\\theta_0}^{(\\textrm{M}2)} & \\alpha_{\\textrm{triv}}^{(\\textrm{M}1)} & \\alpha_{\\textrm{triv}}^{(\\textrm{M}2)} & \\alpha_{\\hat{\\theta}^{(\\textrm{M}1)}_{n_s,n_I}} & \\alpha_{\\hat{\\theta}^{(\\textrm{M}2)}_{n_s,n_I}} \\\\\n';
end


function footer=latex_table_footer
% Footer for Table 5
footer='\\end{tabular}\n';
caption='\\caption{A comparison of the values of the objective functions for the true value~$\\theta_0$, for the trivial estimator and for the MSM estimator. The four synthetic datasets used are the same as in Tables~\\ref{tb:table1}--\\ref{tb:table4}.}';
%caption=['\\caption' sprintf('{%i estimates for a synthetic dataset with $n_I=%i\\\\times %i=%i$ vertices ($n_p=%i$) and ', 2*nblocks, nrows, ncols, totalvertices, totaledgelocations) '$\\theta_0=(' num2str(lambda(1)) ', ' num2str(lambda(2)) ', ' num2str(lambda(3)) ', ' num2str(mu) ')$.}\n'];
tablelabel='\\label{tb:table5}\n';
footerend='\\end{center}\n\\end{table}\n';

footer=[footer caption tablelabel footerend];
end


function filename=latex_finddataset(n)
filepath='../';
%filepath='Synthetic_datasets_for_estimation/';

switch n
case 1
filename='25x25.mat';
case 2
filename='100x100.mat';
case 3
filename='300x300.mat';
case 4
filename='500x500.mat';
end

filename=[filepath filename];
end


function filename=latex_findestimate(n, k, method)
filepath='Data/';
%filepath='Data_estimates/';

switch n
case 1

    switch k

        case 1
            filename=sprintf('25x25_10_10_0.1_m%i_estimate.mat',method);
        case 2
            filename=sprintf('25x25_50_10_0.1_m%i_estimate.mat',method);
        case 3
            filename=sprintf('25x25_100_10_0.1_m%i_estimate.mat',method);
%{
        case 1
            filename=sprintf('25x25_10_10_0.1_m%i_estimators.mat',method);
        case 2
            filename=sprintf('25x25_50_10_0.1_m%i_estimators.mat',method);
        case 3
            filename=sprintf('25x25_100_10_0.1_m%i_estimators.mat',method);
%}
    end

case 2

switch k

case 1
filename=sprintf('100x100_20_10_0.05_m%i_estimate.mat',method);
case 2
filename=sprintf('100x100_40_10_0.05_m%i_estimate.mat',method);
%{
case 1
filename=sprintf('100x100_20_10_0.05_m%i_estimators.mat',method);
case 2
filename=sprintf('100x100_40_10_0.05_m%i_estimators.mat',method);
%}
end

case 3

switch k

case 1
filename=sprintf('300x300_2_8_0.05_m%i_estimate.mat',method);
case 2
filename=sprintf('300x300_4_4_0.05_m%i_estimate.mat',method);
case 3
filename=sprintf('300x300_8_2_0.05_m%i_estimate.mat',method);
case 4
filename=sprintf('300x300_5_5_0.05_m%i_estimate.mat',method);
%{
case 1
filename=sprintf('300x300_1_6_0.03_m%i_estimators.mat',method);
case 2
filename=sprintf('300x300_6_3_0.03_m%i_estimators.mat',method);
%}
end

case 4
    
    switch k
        case 1
            filename=sprintf('500x500_1_1_0.04_m%i_estimate.mat',method);
        case 2
            filename=sprintf('500x500_5_5_0.04_m%i_estimate.mat',method);
    end

end

filename=[filepath filename];
end
