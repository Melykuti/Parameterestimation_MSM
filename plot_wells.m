function plot_wells(wells, nrows, ncols, ncolors, originalrows, shiftedrows, d, d_sqrt32, edgecolor, wells_on)

% PLOT_WELLS.m is called by PLOT_PANEL.m to draw the vertices of the grid,
% i.e. the wells or cavities.
%
% wells_on==0 displays tiny thin outlines in color edgecolor (i.e. white) with black interior for vertices
% wells_on==1 displays balls of color with identical outline
% wells_on==2 displays thin outlines in color edgecolor (i.e. gray) with black interior for empty
% wells and balls of color with identical outline
%
% Felix Beck, Bence Melykuti (University of Freiburg, Germany)
% 26/4/2016, 7/2/2017

switch wells_on
    case 0
        %color=reshape(wells(originalrows,:,2:4),[3 1 1]);
        for j=1:ncols
        for i=1:length(originalrows)
            rectangle('Position',[d_sqrt32*originalrows(i)+0.25,d*j+0.25,0.5,0.5], 'Curvature',[1,1], 'FaceColor','k', 'EdgeColor',edgecolor);
        end
        for i=1:length(shiftedrows)
            rectangle('Position',[d_sqrt32*shiftedrows(i)+0.25,d*(j+0.5)+0.25,0.5,0.5], 'Curvature',[1,1], 'FaceColor','k', 'EdgeColor',edgecolor);
        end
        end
    case 1
        [indwellsx indwellsy]=find(any(wells(originalrows(:),:,2:1+ncolors),3));
        for i=1:length(indwellsx)
            color=reshape(wells(originalrows(indwellsx(i)),indwellsy(i),2:1+ncolors),[ncolors 1 1]);
            rectangle('Position',[d_sqrt32*originalrows(indwellsx(i)),d*indwellsy(i),1,1], 'Curvature',[1,1], 'FaceColor',plot_colorshift(color), 'EdgeColor',plot_colorshift(color));
        end

        [indwellsx indwellsy]=find(any(wells(shiftedrows(:),:,2:1+ncolors),3));
        for i=1:length(indwellsx)
            color=reshape(wells(shiftedrows(indwellsx(i)),indwellsy(i),2:1+ncolors),[ncolors 1 1]);
            rectangle('Position',[d_sqrt32*shiftedrows(indwellsx(i)),d*(indwellsy(i)+0.5),1,1], 'Curvature',[1,1], 'FaceColor',plot_colorshift(color), 'EdgeColor',plot_colorshift(color));
        end
    case 2
        % balls of color with identical outline:
        [indwellsx indwellsy]=find(any(wells(originalrows(:),:,2:1+ncolors),3));
        for i=1:length(indwellsx)
            color=reshape(wells(originalrows(indwellsx(i)),indwellsy(i),2:1+ncolors),[ncolors 1 1]);
            rectangle('Position',[d_sqrt32*originalrows(indwellsx(i)),d*indwellsy(i),1,1], 'Curvature',[1,1], 'FaceColor',plot_colorshift(color), 'EdgeColor',plot_colorshift(color));
        end

        [indwellsx indwellsy]=find(any(wells(shiftedrows(:),:,2:1+ncolors),3));
        for i=1:length(indwellsx)
            color=reshape(wells(shiftedrows(indwellsx(i)),indwellsy(i),2:1+ncolors),[ncolors 1 1]);
            rectangle('Position',[d_sqrt32*shiftedrows(indwellsx(i)),d*(indwellsy(i)+0.5),1,1], 'Curvature',[1,1], 'FaceColor',plot_colorshift(color), 'EdgeColor',plot_colorshift(color));
        end
        
        % thin white outlines with black interior for empty wells:
        [indwellsx indwellsy]=find(all(1-wells(originalrows(:),:,2:1+ncolors),3));
        for i=1:length(indwellsx)
            %color=reshape(wells(originalrows(indwellsx(i)),indwellsy(i),2:4),[3 1 1]);
            rectangle('Position',[d_sqrt32*originalrows(indwellsx(i)),d*indwellsy(i),1,1], 'Curvature',[1,1], 'FaceColor','k', 'EdgeColor',edgecolor);
        end

        [indwellsx indwellsy]=find(all(1-wells(shiftedrows(:),:,2:1+ncolors),3));
        for i=1:length(indwellsx)
            %color=reshape(wells(shiftedrows(indwellsx(i)),indwellsy(i),2:4),[3 1 1]);
            rectangle('Position',[d_sqrt32*shiftedrows(indwellsx(i)),d*(indwellsy(i)+0.5),1,1], 'Curvature',[1,1], 'FaceColor','k', 'EdgeColor',edgecolor);
        end
end

end
