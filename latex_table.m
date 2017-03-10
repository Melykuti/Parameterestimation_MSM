function table_text=latex_table(n)

% LATEX_TABLE.m turns saved data into LaTeX code for the tables 1-4 in the paper.
%
% In the resulting LaTeX code:
% 1) Replace every occurrence of e- and e+ in exponential notations with \mathrm{e-} and \mathrm{e+}.
% 2) In captions, replace leading numbers with text as appropriate, check plural.
%
% n==1 for Table 1, 25x25.mat
% n==2 for Table 2, 100x100.mat
% n==3 for Table 3, 300x300.mat
% n==4 for Table 4, 500x500.mat
%
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 7-9/2/2017

    nblocks= 3*(n==1) + 2*(n==2) + 4*(n==3) + 2*(n==4); % how many blocks the table consists of
    header=latex_table_header;
    footer=latex_table_footer(n, nblocks);
    block='';
    for k=1:nblocks
        blocktext=latex_block(n, nblocks, k);
        block=[block blocktext];
    end
    table_text=[header block footer];


    fileID=fopen('Text/table.tex','w');
    fprintf(fileID,table_text);
    fclose(fileID);
%}
end


function blocktext=latex_block(n, nblocks, k)
% k is the row index of the block within the table
    r=cell(4,1); % array of rows
    filename0=latex_finddataset(n);
    filename1=latex_findestimate(n, k, 1); filename2=latex_findestimate(n, k, 2);
    load(filename0);
    load(filename1);

    r{1}=[num2str(input_values{1,2}) '&' num2str(input_values{2,2}) '&' num2str(input_values{3,2}) '&'];
    r{2}='&&& '; r{3}=r{2}; r{4}=r{2};

    v=[estimators{1,2} estimators{2,2}];
    perc=100*abs((v-true_values)./true_values);

    for i=1:4
        r{i}=[r{i} num2str(true_values(i)) ' & ']; % column of \theta_0

        if i<=3 % column of \hat{\theta}^{(\textrm{M}1)}_{n_s,n_I}
            r{i}=[r{i} sprintf('%5.4f & ', estimators{1,2}(i))];
        else
            r{i}=[r{i} sprintf('%5.4f & ', estimators{2,2})];
        end

        r{i}=[r{i} sprintf('%.2f', perc(i)) '\\%% & ']; % column of d^{(\textrm{M}1)}

        switch i % column of \alpha_{\hat{\theta}^{(\textrm{M}1)}_{n_s,n_I}}
            case 1
                r{i}=[r{i} sprintf('%.3g & ', errors{2,2})];
            case 4
                r{i}=[r{i} '(' sprintf('%.3g',time{2}) '\\,\\textrm{s})&'];
            otherwise
                r{i}=[r{i} '& '];
        end

%{
        if i==1 % column of \alpha_{\hat{\theta}^{(\textrm{M}1)}_{n_s,n_I}}
            r{i}=[r{i} sprintf('%.3g & ', errors{2,2})];
        else
            %{
            switch i
                case 2
                    r{i}=[r{i} '& '];
                case 3
                    r{i}=[r{i} '\\textrm{time (s)}& '];
                case 4
                    r{i}=[r{i} sprintf('%.3g & ',time{2})];
            end
            %}
            %r{i}=[r{i} '& '];
        end
%}
    end

    load(filename2);

    v=[estimators{1,2} estimators{2,2}];
    perc=100*abs((v-true_values)./true_values);

    for i=1:4
        if i<=3 % column of \hat{\theta}^{(\textrm{M}2)}_{n_s,n_I}
            r{i}=[r{i} sprintf('%5.4f & ', estimators{1,2}(i))];
        else
            r{i}=[r{i} sprintf('%5.4f & ', estimators{2,2})];
        end

        r{i}=[r{i} sprintf('%.2f', perc(i)) '\\%% & ']; % column of d^{(\textrm{M}2)}

        switch i % column of \alpha_{\hat{\theta}^{(\textrm{M}2)}_{n_s,n_I}}
            case 1
                r{i}=[r{i} sprintf('%.3g', errors{2,2})]; % no & after last entry
            case 4
                r{i}=[r{i} '(' sprintf('%.3g',time{2}) '\\,\\textrm{s})'];
        end

%{
        if i==1 % column of \alpha_{\hat{\theta}^{(\textrm{M}2)}_{n_s,n_I}}
            r{i}=[r{i} sprintf('%.3g', errors{2,2})]; % no & after last entry
        else
            switch i
                case 3
                    r{i}=[r{i} '\\textrm{time (s)}'];
                case 4
                    r{i}=[r{i} sprintf('%.3g',time{2})];
            end
        end
%}

        if k==nblocks & i==4
            r{i}=[r{i} '\n'];
        else
            r{i}=[r{i} '\\\\\n'];
        end
    end
    
    if k<nblocks
       r{4}=[r{4} '\\hline\n'];
    end

    blocktext=[r{1} r{2} r{3} r{4}];
end


function header=latex_table_header
% Header for Tables 1-4
header='\\begin{table}\n\\begin{center}\n\\begin{tabular}{ >{$}c<{$} >{$}c<{$} >{$}c<{$} >{$}c<{$} | >{$}r<{$} >{$}r<{$} >{$}r<{$} >{$}r<{$} >{$}r<{$} >{$}r<{$}}\nn_s & n_{\\textrm{opt}} & \\mu_{\\textrm{max}} & \\theta_0 & \\hat{\\theta}^{(\\textrm{M}1)}_{n_s,n_I} & d^{(\\textrm{M}1)} & \\alpha_{\\hat{\\theta}^{(\\textrm{M}1)}_{n_s,n_I}} & \\hat{\\theta}^{(\\textrm{M}2)}_{n_s,n_I} & d^{(\\textrm{M}2)} & \\alpha_{\\hat{\\theta}^{(\\textrm{M}2)}_{n_s,n_I}}\\\\ \\hline\n';
end


function footer=latex_table_footer(n, nblocks)

filename=latex_finddataset(n);
load(filename);
[nrows, ncols, ncolors, originalrows, shiftedrows] = plot_preprocessing(wells_after_contamination, shape);
[originalrows, shiftedrows, max_edges, totalvertices, totaledgelocations] = createsynthdata_determ(nrows, ncols, ncolors, shape, wells_after_contamination(:,:,1));
lambda=true_values(1:3);
mu=true_values(4);
density_before=sum(sum(wells_before_contamination(:,:,2:end),1),2)/totalvertices
density_after=sum(sum(wells_after_contamination(:,:,2:end),1),2)/totalvertices
density_edges=sum(sum(sum(open_edges)))/totaledgelocations

% Footer for Tables 1-4
footer='\\end{tabular}\n';
%%caption=['\\caption' sprintf('{%i estimates for a synthetic dataset with $n_I=%i\\\\times %i=%i$ vertices ($n_p=%i$) and $\\\\theta_0=(%f, %f, %f, %f)$.}\n', 2*nblocks, nrows, ncols, totalvertices, totaledgelocations, lambda(1), lambda(2), lambda(3), mu)];
%caption=['\\caption' sprintf('{%i estimates for a synthetic dataset with $n_I=%i\\\\times %i=%i$ vertices ($n_p=%i$) and ', 2*nblocks, nrows, ncols, totalvertices, totaledgelocations) '$\\theta_0=(' num2str(lambda(1)) ', ' num2str(lambda(2)) ', ' num2str(lambda(3)) ', ' num2str(mu) ')$.}\n'];
caption=['\\caption' sprintf('{%i estimates for a synthetic dataset with $n_I=%i\\\\times %i=%i$ vertices ($n_p=%i$) and ', 2*nblocks, nrows, ncols, totalvertices, totaledgelocations) '$\\theta_0=(' num2str(lambda(1)) ', ' num2str(lambda(2)) ', ' num2str(lambda(3)) ', ' num2str(mu) ')$. ' sprintf('In this synthetic dataset, the relative frequency of the incidence of the three colours in the seeding is $\\\\bar{\\\\mathcal{X}}=(%.3g, %.3g, %.3g)$, while in the contamination-impacted observed data, it is $\\\\bar{\\\\mathcal{Y}}=(%.3g, %.3g, %.3g)$. The relative frequency of adjacent vertices having an open edge between them is $%.3g$.', density_before(1), density_before(2), density_before(3), density_after(1), density_after(2), density_after(3), density_edges) '}\n'];
tablelabel=['\\label' sprintf('{tb:table%i}', n) '\n'];
footerend='\\end{center}\n\\end{table}\n';

footer=[footer caption tablelabel footerend];
end


function filename=latex_finddataset(n)
filepath='../Synthetic_datasets_for_estimation/';
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
filepath='Data_estimates/';

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
