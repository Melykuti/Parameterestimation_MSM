function [edgelist, components1]=contamination(edges, nrows, ncols, originalrows, shiftedrows)

% CONTAMINATION.m creates from the matrix edges the list of edges in a
% 4-row matrix, edgelist and a cell array of connected components by
% breadth-first search.
%
% EXPLANATION OF VARIABLES
%
% edges             Size (nrows,ncols,3) matrix whose entries (i,j,k) are
%                   indicator variables of open edges between well (i,j) and
%                   its neighbor to the (k==1) right, (k==2) right down,
%                   (k==3) left down.
% edgelist          List of the edges (connected neighbors) of size (4, no.
%                   of open edges). Each column stores an edge (v_1, v_2) as
%                   [rowindex(v_1); colindex(v_1); rowindex(v_2); colindex(v_2)].
% originalrows      Rows that are 'further left'. See explanation of shape
%                   in MSM.m.
% shiftedrows       Rows that are 'further right'. See explanation of shape
%                   in MSM.m.
% edgelist1         edgelist converted to size (2, no. of edges) by
%                   sub2ind, i.e. each vertex is indexed by a linear index.
%                   Each column stores a [linearindex(v_1); linearindex(v_2)].
% components1{i}    Matrix which contains in one row the elements of the
%                   i'th connected component (linear indexing). Its purpose
%                   is that each vertex within a connected component will
%                   ultimately have the same color.
% edgelistq1        Queue of edges. Created as a copy of edgelist1, the
%                   visited edges are removed from it until it is empty.
% vtov              Vertices to visit. A row vector of vertices which have
%                   been discovered as neighbors of an already visited
%                   element of the current components1{i}.
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 12/8/2015, 7/2/2017


% We use i=i(:); j=j(:); to enforce column vector form. This is needed when originalrows or shiftedrows consists of one row only and the other is nonempty, i.e. when there are 2 or 3 rows altogether, otherwise i and j in the right-down arrow and left-down arrow would be row vectors.
% In the definitions of edgelist, we enforce originalrows(i) and shiftedrows(i) to be in row vector form by using i'. This is necessary when length(originalrows) or length(shiftedrows) is 1 because i is a column vector.

edgelist=zeros(4,0);

% For right arrows
[i, j]=find(edges(:,:,1)); % Find row and column indices of all nonzero elements of edges(:,:,1)
i=i(:); j=j(:);
% 1st coordinate: row indices of 1st vertices of edges
% 2nd coordinate: column indices of the same vertices
% 3rd coordinate: row indices of 2nd vertices of edges
% 4th coordinate: column indices of the same vertices
edgelist=[edgelist [i';j';i';j'+1]];


% Analogue for right-down arrows
[i, j]=find(edges(shiftedrows(:),:,2)); i=i(:); j=j(:);
edgelist=[edgelist [shiftedrows(i');j';shiftedrows(i')+1;j'+1]];
[i, j]=find(edges(originalrows(:),:,2)); i=i(:); j=j(:);
edgelist=[edgelist [originalrows(i');j';originalrows(i')+1;j']];

% Analogue for left-down arrows
[i, j]=find(edges(shiftedrows(:),:,3)); i=i(:); j=j(:);
edgelist=[edgelist [shiftedrows(i');j';shiftedrows(i')+1;j']];
[i, j]=find(edges(originalrows(:),:,3)); i=i(:); j=j(:);
edgelist=[edgelist [originalrows(i');j';originalrows(i')+1;j'-1]];

% Converting the representation of wells in edgelist from (row, col) to (col-1)*nrows+row (from 2 coordinates to 1).
% The breadth-first search (BFS) implemented this way is faster.
edgelist1=[sub2ind([nrows ncols], edgelist(1,:), edgelist(2,:)); sub2ind([nrows ncols], edgelist(3,:), edgelist(4,:))];


% We start with the first edge in edgelistq1. Put its endpoints into both
% components1{i} and vtov. Remove this edge from edgelistq1.
%
% While vtov not empty: Take first element from vtov. k2: all edges from
% edgelistq1 which have vtov(1) in the second row. Restrict k2 to k which
% are those where the neighbor of vtov(1) in the first row is not in
% components1{i} yet. Put these neighbors into both components1{i} and
% vtov. Remove the k2 edges from edgelistq1 as these have been traversed.
% Repeat with k1: all edges from edgelistq1 which have vtov(1) in the first
% row. Remove vtov(1) from vtov. When vtov is empty, then the component has
% been fully explored. If there is still an edge left in edgelistq1, then
% it starts a new component.

components1=cell(1);
edgelistq1=edgelist1;
i=1;

while size(edgelistq1,2)>0
    % Start the new component with endvertices of first edge of edgelistq1
    components1{i}=[edgelistq1(1,1) edgelistq1(2,1)]; vtov=components1{i};
    edgelistq1(:,1)=[];
    while length(vtov)>0
        k2=find(edgelistq1(2,:)==vtov(1)); % Elements in edgelistq1 where the vertex vtov(1) is the second vertex
        k=[];
        for l=1:length(k2) % Look through the neighbors k2 and only keep the ones that are not in components1{i}
            if any(edgelistq1(1,k2(l))==components1{i})==0
                k=[k k2(l)]; % k has the indices of all neighbors of vtov(1) that are not yet in components1{i} when vtov(1) is in the second row of edgelistq1
            end
        end
        components1{i}=[components1{i} edgelistq1(1,k)];
        vtov=[vtov edgelistq1(1,k)];
        edgelistq1(:,k2)=[];
        
        % Do this also when vtov(1) is in the first row of edgelistq1
        k1=find(edgelistq1(1,:)==vtov(1)); % Elements in edgelistq1 where the vertex vtov(1) is the first vertex
        k=[];
        for l=1:length(k1) % Look through the neighbors k1 and only keep the ones that are not in components1{i}
            if any(edgelistq1(2,k1(l))==components1{i})==0
                k=[k k1(l)]; % k has the indices of all neighbors of vtov(1) that are not yet in components1{i} when vtov(1) is in the first row of edgelistq1
            end
        end
        components1{i}=[components1{i} edgelistq1(2,k)];
        vtov=[vtov edgelistq1(2,k)];
        edgelistq1(:,k1)=[];
        
        vtov(1)=[];
    end % of while length(vtov)>0
    i=i+1;
end % of while size(edgelistq1,2)>0
end
