function data_visualise_publication(n)

% Visualising all data from
% dataset_A.mat (Dataset A)
% dataset_B.mat (Dataset B)
% for use in publication.
%
% n==0 or 1 or 2 or 5 to create Figure A.1 or A.2 or A.3 or A.4 of Appendix
% of arXiv:1604.08908v3
%
% Save image with e.g.
% >> print('data_0.png','-dpng', '-r300')
% >> set(gcf, 'PaperSize', [9 6.5]); print('dat8x12_ver1.pdf','-dpdf','-bestfit')
% >> print('data_0.eps','-depsc')
%
% Bence Melykuti (University of Freiburg, Germany)
% 15/12/2016, 7-11/3/2017

mu_crit=2*sin(pi/18)

if n==1 | n==2
    load('dataset_A.mat'); v=[nmus nlambdas];
elseif n>=3
    load('dataset_B.mat'); v=[8 7];
else % i.e. n==0
    load('dataset_A.mat'); v=[nmus nlambdas];
    lambda_broad=lambda_st; mu_broad=mu_st;
    load('dataset_B.mat'); v=[8 7];
end

if n==0
fig1=figure;
fig1.Units='inches';
fig1.PaperUnits=fig1.Units; % possible for both Units and PaperUnits: 'inches' | 'centimeters' | 'points'
fig1.PaperSize=fig1.OuterPosition(3:4);
end

if n<=4 & n>0
fig1=figure;
f1pos=fig1.OuterPosition;
fig1.OuterPosition(4)=f1pos(4)*2;
fig1.Units='inches';
fig1.PaperUnits=fig1.Units; % possible for both Units and PaperUnits: 'inches' | 'centimeters' | 'points'
fig1.PaperSize=fig1.OuterPosition(3:4);
end

if n==5
fig1=figure;
f1pos=fig1.OuterPosition;
fig1.OuterPosition(3)=f1pos(3)*2;
fig1.OuterPosition(4)=f1pos(4)*2;
fig1.Units='inches';
fig1.PaperUnits=fig1.Units; % possible for both Units and PaperUnits: 'inches' | 'centimeters' | 'points'
fig1.PaperSize=fig1.OuterPosition(3:4);
end

switch n
    case 0
        loglog(lambda_broad,mu_broad,'bx',lambda_st,mu_st,'ko',[min(lambda_st/5) 1],[mu_crit mu_crit],'r');
        %title('Sampled parameter combinations');
        legend('Dataset A','Dataset B','Critical value_ p_c'); xlabel('lambda'); ylabel('mu');
        % Without the underscore after value, the saved pdf would have a
        % large space between p and its subscript c, I don't know why.
        axis([min(lambda_st/5) 1 min(mu_broad/5) 10]);
        print('app_fig1.pdf','-dpdf')
    case 1
%{
        subplot(2,1,1);data_surface(lambda_st,mu_st,av_st,v,'Empirical mean of_ Y_i',1);
        subplot(2,1,2);data_contour(lambda_st,mu_st,av_st,v,'...and its level curves',1);
%}
        subplot(2,1,1);data_surface(lambda_st,mu_st,av_st,v,'Y_i',1);
        subplot(2,1,2);data_contour(lambda_st,mu_st,av_st,v,'',1);
        print('app_fig2.pdf','-dpdf')
    case 2
%{
        subplot(2,1,1);data_surface(lambda_st,mu_st,avnb_st,v,'Empirical mean of_ Y_i Y_j',1);
        subplot(2,1,2);data_contour(lambda_st,mu_st,avnb_st,v,'...and its level curves',1);
%}
        subplot(2,1,1);data_surface(lambda_st,mu_st,avnb_st,v,'Y_i Y_j',1);
        subplot(2,1,2);data_contour(lambda_st,mu_st,avnb_st,v,'',1);
        print('app_fig3.pdf','-dpdf')
    case 3
        subplot(2,1,1);data_surface(lambda_st,mu_st,av_st,v,'Empirical mean of_ Y_i',0);
        subplot(2,1,2);data_contour(lambda_st,mu_st,av_st,v,'...and its level curves',0);
    case 4
        subplot(2,1,1);data_surface(lambda_st,mu_st,avnb_st,v,'Empirical mean of_ Y_i Y_j',0);
        subplot(2,1,2);data_contour(lambda_st,mu_st,avnb_st,v,'...and its level curves',0);
    case 5
%{        
        subplot(2,2,1);data_surface(lambda_st,mu_st,av_st,v,'Empirical mean of_ Y_i',0);
        subplot(2,2,3);data_contour(lambda_st,mu_st,av_st,v,'...and its level curves',0);
        subplot(2,2,2);data_surface(lambda_st,mu_st,avnb_st,v,'Empirical mean of_ Y_i Y_j',0);
        subplot(2,2,4);data_contour(lambda_st,mu_st,avnb_st,v,'...and its level curves',0);
%}
        subplot(2,2,1);data_surface(lambda_st,mu_st,av_st,v,'Y_i',0);
        subplot(2,2,3);data_contour(lambda_st,mu_st,av_st,v,'',0);
        subplot(2,2,2);data_surface(lambda_st,mu_st,avnb_st,v,'Y_i Y_j',0);
        subplot(2,2,4);data_contour(lambda_st,mu_st,avnb_st,v,'',0);
        print('app_fig4.pdf','-dpdf')
end
end

function data_surface(lambda_st,mu_st,data_values,v,caption,criticalline)
 surf(reshape(lambda_st, v),reshape(mu_st, v),reshape(data_values, v));
 title(caption);
 xlabel('lambda');ylabel('mu');
 if criticalline==1
  mu_crit=2*sin(pi/18);
  line([0 0.5 0.5],[mu_crit mu_crit mu_crit],[0 0 1],'Color','red','LineWidth',2);
 end
end

function data_contour(lambda_st,mu_st,data_values,v,caption,criticalline)
 contour3(reshape(lambda_st, v),reshape(mu_st, v),reshape(data_values, v));
 title(caption);
 xlabel('lambda');ylabel('mu');view(0,90);
 if criticalline==1
  mu_crit=2*sin(pi/18);
  line([0 0.5 0.5],[mu_crit mu_crit mu_crit],[0 0 1],'Color','red','LineWidth',2);
 end
end