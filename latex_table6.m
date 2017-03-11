function table_text=latex_table6

% LATEX_TABLE6.m turns saved data into LaTeX code for table 6 in the paper.
%
% In the resulting LaTeX code:
% 1) Replace every occurrence of e- and e+ in exponential notations with \mathrm{e-} and \mathrm{e+}.
%
%
% (version 13/2/2017)
% tic, table_text=latex_table6, toc
% Elapsed time is 1437.633318 seconds.
%
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 7-14/2/2017


header=latex_table_header;
footer=latex_table_footer;
body='\\hline\n'; % the \hline is at top of each block, and not part of header
for n=1:6
    row_text=latex_row(n);
    body=[body row_text];
end
table_text=[header body footer];


    fileID=fopen('Text/table.tex','w');
    fprintf(fileID,table_text);
    fclose(fileID);
%}
end


function r=latex_row(n)

filename0=latex_finddataset(n);
load(filename0);
[nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(wells_after_contamination, shape);
[originalrows, shiftedrows, max_edges, totalvertices, totaledgelocations] = createsynthdata_determ(nrows, ncols, ncolors, shape, wells_after_contamination(:,:,1));
simstep_max=10;
r=sprintf('%i\\\\times%i & %i & %i & %i & ', nrows, ncols, totalvertices, totaledgelocations, simstep_max);

[dsRGB, dsnb] = simcalcs(nrows, ncols, ncolors, originalrows, shiftedrows, wells_after_contamination, totalvertices, totaledgelocations);
dstwocol=0;

alpha_0=zeros(2,2);

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

[alpha_0(method,1), alpha_0(method,2)] = optim_norma(true_values, nrows, ncols, ncolors, shape, originalrows, shiftedrows, simstep_max, dsRGB, dstwocol, dsnb, wells_before_contamination(:,:,1), max_edges, totalvertices, totaledgelocations, rvcolors, rvedges, rplabelled_colors, rplabelled_edges, method);
end

alpha_0
%{
load(filename2);
alpha_triv(2) = errors{1,2};
alpha(2) = errors{2,2};
%}
r=[r sprintf('%.3g & %.3g & %.3g & %.3g', alpha_0(1,1), alpha_0(2,1), alpha_0(1,2), alpha_0(2,2))];

if n==7
    r=[r '\n'];
else
    r=[r '\\\\\n'];
end
end


function header=latex_table_header
% Header for Table 6
    header='\\begin{table}\n\\begin{center}\n\\begin{tabular}{>{$}c<{$} >{$}c<{$} >{$}c<{$} >{$}c<{$} | >{$}c<{$}  >{$}c<{$}  >{$}c<{$}  >{$}c<{$}}\n\\textrm{Size} & n_I & n_p & n_s & \\alpha_{\\theta_0}^{(\\textrm{M}1)} & \\alpha_{\\theta_0}^{(\\textrm{M}2)} & \\widetilde{\\alpha}_{\\theta_0}^{(\\textrm{M}1)} & \\widetilde{\\alpha}_{\\theta_0}^{(\\textrm{M}2)} \\\\\n';
end


function footer=latex_table_footer
% Footer for Table 6
footer='\\end{tabular}\n';
caption='\\caption{Realisations of the objective function~$\\alpha$ for the true parameter value~$\\theta_0$ for different synthetic dataset sizes and of the not normalised variant of the objective function $\\widetilde{\\alpha}(\\eta)=\\eta\\T \\eta$. Here $\\theta_0=(0.03,0.04,0.05,0.02)$ across fresh synthetic datasets.}';
tablelabel='\\label{tb:table6}\n';
footerend='\\end{center}\n\\end{table}\n';

footer=[footer caption tablelabel footerend];
end


function filename=latex_finddataset(n)
%{
filepath='../';
switch n
case 1
filename='25x25.mat';
case 2
filename='100x100.mat';
case 3
filename='300x300.mat';
case 4
filename='500x500.mat';
case 5
filename='661x661.mat';
end
%}
filepath='Synthetic_datasets_for_estimation/';
switch n
case 1
filename='25x25_testdata_m1.mat';
case 2
filename='100x100_testdata_m1.mat';
case 3
filename='300x300_testdata_m1.mat';
case 4
filename='500x500_testdata_m1.mat';
case 5
filename='707x707_testdata_m1.mat';
case 6
filename='1000x1000_testdata_m1.mat';
end
%}
filename=[filepath filename];
end
