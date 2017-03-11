function data_visualise

% Visualising all data from
% dataset_A.mat (Dataset A)
% dataset_B.mat (Dataset B)
%
% Bence Melykuti (University of Freiburg, Germany)
% 15/12/2016, 11/3/2017

mu_crit=2*sin(pi/18)
load('dataset_A.mat'); v=[nmus nlambdas];
lambda_broad=lambda_st;
mu_broad=mu_st;

data_surface(lambda_st,mu_st,av_st,v,'Empirical mean of Y_i',1);
data_contour(lambda_st,mu_st,av_st,v,'Empirical mean of Y_i',1);
data_surface(lambda_st,mu_st,avnb_st,v,'Empirical mean of Y_i Y_j',1);
data_contour(lambda_st,mu_st,avnb_st,v,'Empirical mean of Y_i Y_j',1);

load('dataset_B.mat'); v=[8 7];
data_surface(lambda_st,mu_st,av_st,v,'Empirical mean of Y_i',0);
data_contour(lambda_st,mu_st,av_st,v,'Empirical mean of Y_i',0);
data_surface(lambda_st,mu_st,avnb_st,v,'Empirical mean of Y_i Y_j',0);
data_contour(lambda_st,mu_st,avnb_st,v,'Empirical mean of Y_i Y_j',0);

figure; loglog(lambda_broad,mu_broad,'bx',lambda_st,mu_st,'go',[min(lambda_st/5) 1],[mu_crit mu_crit],'r');
title('Sampled parameter combinations'); legend('Dataset A','Dataset B'); xlabel('lambda'); ylabel('mu');
axis([min(lambda_st/5) 1 min(mu_broad/5) 10]);
end

function data_surface(lambda_st,mu_st,data_values,v,caption,criticalline)
 figure;
 surf(reshape(lambda_st, v),reshape(mu_st, v),reshape(data_values, v));
 %surf(reshape(log10(lambda_st), v),reshape(log10(mu_st), v),reshape(data_values, v));
 title(caption);xlabel('lambda');ylabel('mu');
 if criticalline==1
  mu_crit=2*sin(pi/18);
  line([0 0.5 0.5],[mu_crit mu_crit mu_crit],[0 0 1],'Color','red','LineWidth',2);
 % line(log10([min(lambda_st(:)) 0.5 0.5]),log10([mu_crit mu_crit mu_crit]),[0 0 1],'Color','red','LineWidth',2);
 end
end

function data_contour(lambda_st,mu_st,data_values,v,caption,criticalline)
 figure;
 contour3(reshape(lambda_st, v),reshape(mu_st, v),reshape(data_values, v));
 %contour3(reshape(log10(lambda_st), v),reshape(log10(mu_st), v),reshape(data_values, v));
 title(caption);xlabel('lambda');ylabel('mu');view(0,90);
 if criticalline==1
  mu_crit=2*sin(pi/18);
  line([0 0.5 0.5],[mu_crit mu_crit mu_crit],[0 0 1],'Color','red','LineWidth',2);
 % line(log10([min(lambda_st(:)) 0.5 0.5]),log10([mu_crit mu_crit mu_crit]),[0 0 1],'Color','red','LineWidth',2);
 end
end