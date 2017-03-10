function plot_edges(nrows, ncols, originalrows, shiftedrows, d, d_sqrt32, edgelist, color)

% PLOT_EDGES.m is called by PLOT_PANEL.m to draw edges.
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 26/4/2016, 7/2/2017

edgesx=d_sqrt32*([edgelist(1,:); edgelist(3,:)])+0.5;
edgesy=zeros(2,0);
for i=1:size(edgelist,2)
    edgesy=[edgesy d*([edgelist(2,i); edgelist(4,i)]+0.5*[any(shiftedrows==edgelist(1,i)); any(shiftedrows==edgelist(3,i))])+0.5];
end
plot(edgesx,edgesy,'LineWidth',1.5,'Color',color)
end
