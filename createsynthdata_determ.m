function [originalrows, shiftedrows, max_edges, totalvertices, totaledgelocations] = createsynthdata_determ(nrows, ncols, ncolors, shape, max_vertices)

% CREATESYNTHDATA_DETERM.m creates the deterministic framework needed for
% the synthetic data: shiftedrows, originalrows, max_edges, totalvertices
% [not needed for method==1], totaledgelocations [not needed for method==1]
%
% method==1 uses the Bernoulli variables,
% method==2 uses random permutations
%
% EXPLANATION OF VARIABLES
%
% originalrows      Row vector of row numbers of vertices that are in
%                   default position, 'further left'.
% shiftedrows       Row vector of row numbers of vertices that are 'further
%                   right' (shifted right by a half unit relative to originalrows).
% max_edges         Matrix of size (nrows,ncols,3) that indicates if the edge
%                   can possibly be open between vertices (i,j) and three
%                   of its neighbors. It is 1 when both endpoints are in
%                   the area of interest (i.e. max_vertices(.,.)==1
%                   for both endpoints), otherwise 0. The last dimension
%                   represents the three possible directions of edges. Edge
%                   direction: 1==right, 2==right down, 3==left down.
%                   E.g. if there could be contamination between vertex
%                   (i,j) and its right down neighbor, max_edges(i,j,2)==1.
%                   If not (i.e. the neighbor lies outside the area of
%                   interest), then max_edges (i,j,2)==0.
%
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 2/2/2017

% Create originalrows & shiftedrows
if shape==0;
    originalrows=2:2:nrows;
    shiftedrows=1:2:nrows;
else % shape==1
    originalrows=1:2:nrows;
    shiftedrows=2:2:nrows;
end


% Edge direction: 1=right, 2=right down, 3=left down
%                   
% If an edge were to connect a vertex outside and one inside the area of
% interest, it must be set to closed (0). We first allow every edge to be
% open, then zero those which cannot physically be:
% (originalrows, 1): no left down;
% (shiftedrows, end): no right, no right down;
% (originalrows, end): no right;
% (end,:): no right down, no left down
max_edges=ones(nrows,ncols,3);
max_edges(originalrows,1,3)=0;
max_edges(shiftedrows,end,[1 2])=0;
max_edges(originalrows,end,1)=0;
max_edges(end,:,[2 3])=0;

% This part is relevant only if the real data represents an irregular area
% of interest. Otherwise it can be commented out. BEGIN

% Checking if the neighbors for each well are part of area of interest and
% if not, close the corresponding edge.
for i=1:nrows
    for j=1:ncols
        if max_vertices(i,j)==1 % 'Area of interest'
            if any(i==originalrows) % If row is part of originalrows
                if j<ncols && max_vertices(i,j+1)==0 % If right neighbor is not in area of interest
                    max_edges(i,j,1)=0; % then right edge cannot exist
                end
                if i<nrows && max_vertices(i+1,j)==0 % If right down neighbor is not in area of interest
                    max_edges(i,j,2)=0; % then right down edge cannot exist
                end
                if i<nrows && j>1 && max_vertices(i+1,j-1)==0 % If left down neighbor is not in area of interest
                    max_edges(i,j,3)=0; % then left down edge cannot exist
                end
            else  % Analogue if row is part of shiftedrows
                if j<ncols && max_vertices(i,j+1)==0
                    max_edges(i,j,1)=0;
                end
                if j<ncols && i<nrows && max_vertices(i+1,j+1)==0
                    max_edges(i,j,2)=0;
                end
                if i<nrows && max_vertices(i+1,j)==0
                    max_edges(i,j,3)=0;
                end
            end
        else
            max_edges(i,j,:)=0; % If (i,j) is not in area of interest, all three edges are zero
        end
    end
end

% END End of possible commenting out.

totalvertices=sum(sum(max_vertices));
totaledgelocations=sum(sum(sum(max_edges)));

end
