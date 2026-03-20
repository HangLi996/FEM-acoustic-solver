function plotMesh(h, nodes, elements, conf, B, limiter, col, v, loc)
%   Input:  h       ... figure handle
%           nodes   ... nodal coordinates
%           elements... elements
%           B       ... boundary nodes and elements
%           limiter ... show nodes, element or boundaries
%           col     ... color of the lines
%           v       ... view angle
%           loc     ... location of the plot window
%   Output: ~
%
%   Plot mesh.
%
%   by Lennart Moheit 18/08/2014
n=limiter(1);
e=limiter(2);
b=limiter(3);

figure(h)
movegui(h, loc)
hold all

%% Restriction to boundary nodes
%{
    nodes.all = [];
    nodes.nall = [];
    elements.vert = [];
    elements.all = [];
    for qq = B.innerGamma
        nodes.all = [nodes.all; B.nodes{qq}];
        nodes.nall = [nodes.nall; B.n{qq}];
        elements.vert = [elements.vert; B.vert{qq}];
        elements.all = [elements.all; B.elements{qq}];
    end
    elements.vert = changem(elements.vert, 1:length(nodes.nall), nodes.nall);
    elements.all = changem(elements.all, 1:length(nodes.nall), nodes.nall);
    elements.type = elements.type ./ elements.type * 8;
%}

%%
N.n = size(nodes.all, 1);
N.el = size(elements.vert, 1);

dim = getDimension(elements.type);

%xlim([ min(nodes.all(:,1))-0.1*abs(max(nodes.all(:,1))-min(nodes.all(:,1))) max(nodes.all(:,1))+0.1*abs(max(nodes.all(:,1))-min(nodes.all(:,1))) ]);  
%ylim([ min(nodes.all(:,2))-0.1*abs(max(nodes.all(:,2))-min(nodes.all(:,2))) max(nodes.all(:,2))+0.1*abs(max(nodes.all(:,2))-min(nodes.all(:,2))) ]);
%zlim([ min(nodes.all(:,3))-0.1*abs(max(nodes.all(:,3))-min(nodes.all(:,3))) max(nodes.all(:,3))+0.1*abs(max(nodes.all(:,3))-min(nodes.all(:,3))) ]);
%% Nodes
plot3(nodes.all(:,1), nodes.all(:,2), nodes.all(:,3),'kx')
if n==true
    text(nodes.all(:,1), nodes.all(:,2), nodes.all(:,3), cellstr(num2str([1:N.n]')), 'VerticalAlignment','bottom', 'HorizontalAlignment','right')
end
%% Lines
    % in gmsh line coordinates as well as nodes are given, so it could be a better
    % solution to read this data
    N.lines=0;
    for el=1:N.el % for each element
        elnodes.n = elements.all(el, :); % respective nodes per element el
            elnodes.n = elnodes.n(find(elnodes.n~=0)); % cut zero columns for mixed meshes and elements with more vertices
            N.nperel = size(elnodes.n, 2); % Number of nodes per element
        elnodes.nodes = nodes.all(nodes.nind(elnodes.n, 2), :); % respective nodal coordinates per el
        %elnodes.type = elements.type(el); % elements.type
        switch unique(elements.type) % http://geuz.org/gmsh/doc/texinfo/gmsh.html
            case 1 % 2-node line
                line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 2 % 3-node triangle
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(1,1)],[elnodes.nodes(3,2),elnodes.nodes(1,2)],[elnodes.nodes(3,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(1,1)],[elnodes.nodes(3,2),elnodes.nodes(1,2)],[elnodes.nodes(3,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                end
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 3 % 4-node quadrangle
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(4,1)],[elnodes.nodes(2,2),elnodes.nodes(4,2)],[elnodes.nodes(2,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(3,1)],[elnodes.nodes(4,2),elnodes.nodes(3,2)],[elnodes.nodes(4,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(1,1)],[elnodes.nodes(3,2),elnodes.nodes(1,2)],[elnodes.nodes(3,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(4,1)],[elnodes.nodes(3,2),elnodes.nodes(4,2)],[elnodes.nodes(4,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(1,1)],[elnodes.nodes(4,2),elnodes.nodes(1,2)],[elnodes.nodes(3,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                end
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)), mean(elnodes.nodes(:, 2)), mean(elnodes.nodes(:, 3)), num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 4 % 4-node tetrahedron
 % COMSOL
% GMSH?
                line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                line([elnodes.nodes(3,1),elnodes.nodes(1,1)],[elnodes.nodes(3,2),elnodes.nodes(1,2)],[elnodes.nodes(3,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                line([elnodes.nodes(3,1),elnodes.nodes(4,1)],[elnodes.nodes(3,2),elnodes.nodes(4,2)],[elnodes.nodes(3,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                line([elnodes.nodes(4,1),elnodes.nodes(1,1)],[elnodes.nodes(4,2),elnodes.nodes(1,2)],[elnodes.nodes(4,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                line([elnodes.nodes(2,1),elnodes.nodes(4,1)],[elnodes.nodes(2,2),elnodes.nodes(4,2)],[elnodes.nodes(2,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                line([elnodes.nodes(2,1),elnodes.nodes(4,1)],[elnodes.nodes(2,2),elnodes.nodes(4,2)],[elnodes.nodes(2,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),mean(elnodes.nodes(:, 3)),num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 5 % 8-node hexahedron
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                    disp('! NO COMSOL DATA FOR ELEMENT TYPE 5');
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(4,1)],[elnodes.nodes(3,2),elnodes.nodes(4,2)],[elnodes.nodes(3,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(1,1)],[elnodes.nodes(4,2),elnodes.nodes(1,2)],[elnodes.nodes(4,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(5,1),elnodes.nodes(6,1)],[elnodes.nodes(5,2),elnodes.nodes(6,2)],[elnodes.nodes(5,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(6,1),elnodes.nodes(7,1)],[elnodes.nodes(6,2),elnodes.nodes(7,2)],[elnodes.nodes(6,3),elnodes.nodes(7,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(7,1),elnodes.nodes(8,1)],[elnodes.nodes(7,2),elnodes.nodes(8,2)],[elnodes.nodes(7,3),elnodes.nodes(8,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(8,1),elnodes.nodes(5,1)],[elnodes.nodes(8,2),elnodes.nodes(5,2)],[elnodes.nodes(8,3),elnodes.nodes(5,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(1,1),elnodes.nodes(5,1)],[elnodes.nodes(1,2),elnodes.nodes(5,2)],[elnodes.nodes(1,3),elnodes.nodes(5,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(6,1)],[elnodes.nodes(2,2),elnodes.nodes(6,2)],[elnodes.nodes(2,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(7,1)],[elnodes.nodes(3,2),elnodes.nodes(7,2)],[elnodes.nodes(3,3),elnodes.nodes(7,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(8,1)],[elnodes.nodes(4,2),elnodes.nodes(8,2)],[elnodes.nodes(4,3),elnodes.nodes(8,3)],'Color',col); N.lines=N.lines+1;
                end
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),mean(elnodes.nodes(:, 3)),num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 8 % 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
                % COMSOL
% GMSH?
                line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)), mean(elnodes.nodes(:, 2)), num2str(el), 'Color', [0.7 0.7 0.7]);  
                end
            case 9 % 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(5,1)],[elnodes.nodes(3,2),elnodes.nodes(5,2)],[elnodes.nodes(3,3),elnodes.nodes(5,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(5,1),elnodes.nodes(6,1)],[elnodes.nodes(5,2),elnodes.nodes(6,2)],[elnodes.nodes(5,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(6,1),elnodes.nodes(4,1)],[elnodes.nodes(6,2),elnodes.nodes(4,2)],[elnodes.nodes(6,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(1,1)],[elnodes.nodes(4,2),elnodes.nodes(1,2)],[elnodes.nodes(4,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                    line([elnodes.nodes(1,1),elnodes.nodes(4,1)],[elnodes.nodes(1,2),elnodes.nodes(4,2)],[elnodes.nodes(1,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(2,1)],[elnodes.nodes(4,2),elnodes.nodes(2,2)],[elnodes.nodes(4,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(5,1)],[elnodes.nodes(2,2),elnodes.nodes(5,2)],[elnodes.nodes(2,3),elnodes.nodes(5,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(5,1),elnodes.nodes(3,1)],[elnodes.nodes(5,2),elnodes.nodes(3,2)],[elnodes.nodes(5,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(6,1)],[elnodes.nodes(3,2),elnodes.nodes(6,2)],[elnodes.nodes(3,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(6,1),elnodes.nodes(1,1)],[elnodes.nodes(6,2),elnodes.nodes(1,2)],[elnodes.nodes(6,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                end
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),mean(elnodes.nodes(:, 3)),num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 10 % 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(6,1)],[elnodes.nodes(3,2),elnodes.nodes(6,2)],[elnodes.nodes(3,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(6,1),elnodes.nodes(9,1)],[elnodes.nodes(6,2),elnodes.nodes(9,2)],[elnodes.nodes(6,3),elnodes.nodes(9,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(9,1),elnodes.nodes(7,1)],[elnodes.nodes(9,2),elnodes.nodes(7,2)],[elnodes.nodes(9,3),elnodes.nodes(7,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(7,1),elnodes.nodes(4,1)],[elnodes.nodes(7,2),elnodes.nodes(4,2)],[elnodes.nodes(7,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(1,1)],[elnodes.nodes(4,2),elnodes.nodes(1,2)],[elnodes.nodes(4,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                    line([elnodes.nodes(1,1),elnodes.nodes(5,1)],[elnodes.nodes(1,2),elnodes.nodes(5,2)],[elnodes.nodes(1,3),elnodes.nodes(5,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(5,1),elnodes.nodes(2,1)],[elnodes.nodes(5,2),elnodes.nodes(2,2)],[elnodes.nodes(5,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(6,1)],[elnodes.nodes(2,2),elnodes.nodes(6,2)],[elnodes.nodes(2,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(6,1),elnodes.nodes(3,1)],[elnodes.nodes(6,2),elnodes.nodes(3,2)],[elnodes.nodes(6,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(7,1)],[elnodes.nodes(3,2),elnodes.nodes(7,2)],[elnodes.nodes(3,3),elnodes.nodes(7,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(7,1),elnodes.nodes(4,1)],[elnodes.nodes(7,2),elnodes.nodes(4,2)],[elnodes.nodes(7,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(8,1)],[elnodes.nodes(4,2),elnodes.nodes(8,2)],[elnodes.nodes(4,3),elnodes.nodes(8,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(8,1),elnodes.nodes(1,1)],[elnodes.nodes(8,2),elnodes.nodes(1,2)],[elnodes.nodes(8,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                end
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),mean(elnodes.nodes(:, 3)),num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 11 % 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                    disp('! NO COMSOL DATA FOR ELEMENT TYPE 11');
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
% THESE GIVES JUST LINES BETWEEN THE VERTICES IGNORING THE ADDITIONAL
% POINTS ON THE EDGES!! IS THAT ENOUGH??
                    line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(4,1)],[elnodes.nodes(2,2),elnodes.nodes(4,2)],[elnodes.nodes(2,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(1,1)],[elnodes.nodes(4,2),elnodes.nodes(1,2)],[elnodes.nodes(4,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(1,1),elnodes.nodes(3,1)],[elnodes.nodes(1,2),elnodes.nodes(3,2)],[elnodes.nodes(1,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(3,1)],[elnodes.nodes(4,2),elnodes.nodes(3,2)],[elnodes.nodes(4,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                end
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),mean(elnodes.nodes(:, 3)),num2str(el),'Color',[0.7 0.7 0.7]);
                end                
            case 12 % 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                    disp('! NO COMSOL DATA FOR ELEMENT TYPE 12');
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                    
% DO WE REALLY NEED ALL THE SMALL LINES? TO SHOW THE MESH, WE JUST NEED
% STRAIGHT LINES BETWEEN THE VERTICES IGNORING THE POINTS ON THE EDGES
                    line([elnodes.nodes(1,1),elnodes.nodes(9,1)],[elnodes.nodes(1,2),elnodes.nodes(9,2)],[elnodes.nodes(1,3),elnodes.nodes(9,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(9,1),elnodes.nodes(2,1)],[elnodes.nodes(9,2),elnodes.nodes(2,2)],[elnodes.nodes(9,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(13,1)],[elnodes.nodes(2,2),elnodes.nodes(13,2)],[elnodes.nodes(2,3),elnodes.nodes(13,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(13,1),elnodes.nodes(6,1)],[elnodes.nodes(13,2),elnodes.nodes(6,2)],[elnodes.nodes(13,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(6,1),elnodes.nodes(17,1)],[elnodes.nodes(6,2),elnodes.nodes(17,2)],[elnodes.nodes(6,3),elnodes.nodes(17,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(17,1),elnodes.nodes(5,1)],[elnodes.nodes(17,2),elnodes.nodes(5,2)],[elnodes.nodes(17,3),elnodes.nodes(5,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(5,1),elnodes.nodes(11,1)],[elnodes.nodes(5,2),elnodes.nodes(11,2)],[elnodes.nodes(5,3),elnodes.nodes(11,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(11,1),elnodes.nodes(1,1)],[elnodes.nodes(11,2),elnodes.nodes(1,2)],[elnodes.nodes(11,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(3,1),elnodes.nodes(14,1)],[elnodes.nodes(3,2),elnodes.nodes(14,2)],[elnodes.nodes(3,3),elnodes.nodes(14,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(14,1),elnodes.nodes(4,1)],[elnodes.nodes(14,2),elnodes.nodes(4,2)],[elnodes.nodes(14,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(4,1),elnodes.nodes(16,1)],[elnodes.nodes(4,2),elnodes.nodes(16,2)],[elnodes.nodes(4,3),elnodes.nodes(16,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(16,1),elnodes.nodes(8,1)],[elnodes.nodes(16,2),elnodes.nodes(8,2)],[elnodes.nodes(16,3),elnodes.nodes(8,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(8,1),elnodes.nodes(20,1)],[elnodes.nodes(8,2),elnodes.nodes(20,2)],[elnodes.nodes(8,3),elnodes.nodes(20,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(20,1),elnodes.nodes(7,1)],[elnodes.nodes(20,2),elnodes.nodes(7,2)],[elnodes.nodes(20,3),elnodes.nodes(7,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(7,1),elnodes.nodes(15,1)],[elnodes.nodes(7,2),elnodes.nodes(15,2)],[elnodes.nodes(7,3),elnodes.nodes(15,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(15,1),elnodes.nodes(3,1)],[elnodes.nodes(15,2),elnodes.nodes(3,2)],[elnodes.nodes(15,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(1,1),elnodes.nodes(10,1)],[elnodes.nodes(1,2),elnodes.nodes(10,2)],[elnodes.nodes(1,3),elnodes.nodes(10,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(10,1),elnodes.nodes(4,1)],[elnodes.nodes(10,2),elnodes.nodes(4,2)],[elnodes.nodes(10,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(5,1),elnodes.nodes(18,1)],[elnodes.nodes(5,2),elnodes.nodes(18,2)],[elnodes.nodes(5,3),elnodes.nodes(18,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(18,1),elnodes.nodes(8,1)],[elnodes.nodes(18,2),elnodes.nodes(8,2)],[elnodes.nodes(18,3),elnodes.nodes(8,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(6,1),elnodes.nodes(19,1)],[elnodes.nodes(6,2),elnodes.nodes(19,2)],[elnodes.nodes(6,3),elnodes.nodes(19,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(19,1),elnodes.nodes(7,1)],[elnodes.nodes(19,2),elnodes.nodes(7,2)],[elnodes.nodes(19,3),elnodes.nodes(7,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(2,1),elnodes.nodes(12,1)],[elnodes.nodes(2,2),elnodes.nodes(12,2)],[elnodes.nodes(2,3),elnodes.nodes(12,3)],'Color',col); N.lines=N.lines+1;
                    line([elnodes.nodes(12,1),elnodes.nodes(3,1)],[elnodes.nodes(12,2),elnodes.nodes(3,2)],[elnodes.nodes(12,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
                end
                %
                if e == true
                    text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),mean(elnodes.nodes(:, 3)),num2str(el),'Color',[0.7 0.7 0.7]);
                end
            case 16 % 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)
%                 line([elnodes.nodes(1,1),elnodes.nodes(2,1)],[elnodes.nodes(1,2),elnodes.nodes(2,2)],[elnodes.nodes(1,3),elnodes.nodes(2,3)],'Color',col); N.lines=N.lines+1;
%                 line([elnodes.nodes(2,1),elnodes.nodes(3,1)],[elnodes.nodes(2,2),elnodes.nodes(3,2)],[elnodes.nodes(2,3),elnodes.nodes(3,3)],'Color',col); N.lines=N.lines+1;
%                 line([elnodes.nodes(3,1),elnodes.nodes(5,1)],[elnodes.nodes(3,2),elnodes.nodes(5,2)],[elnodes.nodes(3,3),elnodes.nodes(5,3)],'Color',col); N.lines=N.lines+1;
%                 line([elnodes.nodes(5,1),elnodes.nodes(8,1)],[elnodes.nodes(5,2),elnodes.nodes(8,2)],[elnodes.nodes(5,3),elnodes.nodes(8,3)],'Color',col); N.lines=N.lines+1;
%                 line([elnodes.nodes(8,1),elnodes.nodes(7,1)],[elnodes.nodes(8,2),elnodes.nodes(7,2)],[elnodes.nodes(8,3),elnodes.nodes(7,3)],'Color',col); N.lines=N.lines+1;
%                 line([elnodes.nodes(7,1),elnodes.nodes(6,1)],[elnodes.nodes(7,2),elnodes.nodes(6,2)],[elnodes.nodes(7,3),elnodes.nodes(6,3)],'Color',col); N.lines=N.lines+1;
%                 line([elnodes.nodes(6,1),elnodes.nodes(4,1)],[elnodes.nodes(6,2),elnodes.nodes(4,2)],[elnodes.nodes(6,3),elnodes.nodes(4,3)],'Color',col); N.lines=N.lines+1;
%                 line([elnodes.nodes(4,1),elnodes.nodes(1,1)],[elnodes.nodes(4,2),elnodes.nodes(1,2)],[elnodes.nodes(4,3),elnodes.nodes(1,3)],'Color',col); N.lines=N.lines+1;
%                 %
%                 if e == true
%                     text(mean(elnodes.nodes(:, 1)),mean(elnodes.nodes(:, 2)),mean(elnodes.nodes(:, 3)),num2str(el),'Color',[0.7 0.7 0.7]);
%                 end
            otherwise disp(['! WARNING: Unable to plot element No. ' num2str(el) ' yet.'])
        end
    end

%% Boundary Lines
if b == true
    % Boundary normal vectors
    if iscell(B.elcenterpoint)
        cmap = cool(size(B.elcenterpoint,2));
        for qq = 1:size(B.elcenterpoint,2)
        	quiver3(B.elcenterpoint{qq}(:, 1), B.elcenterpoint{qq}(:, 2), B.elcenterpoint{qq}(:, 3), B.elnormalvector{qq}(:, 1), B.elnormalvector{qq}(:, 2), B.elnormalvector{qq}(:, 3), 0.2, 'Color', cmap(qq,:))
            for pp = 1:size(B.elements{qq}, 1)
                switch dim
                    case 1 % 1D
                        disp('! WARNING: No 1D support yet.')
                    case 2 % 2D
                        text(mean([nodes.all(nodes.nind(B.elements{qq}(pp,1), 2),1), nodes.all(nodes.nind(B.elements{qq}(pp,end), 2),1)]),...
                            mean([nodes.all(nodes.nind(B.elements{qq}(pp,1), 2),2), nodes.all(nodes.nind(B.elements{qq}(pp,end), 2),2)]),...
                            ['\Gamma_{' num2str(pp) '}'],'Color', cmap(qq,:));
                    case 3 % 3D
                        text(mean([nodes.all(nodes.nind(B.elements{qq}(pp,1), 2),1), nodes.all(nodes.nind(B.elements{qq}(pp,2), 2),1), nodes.all(nodes.nind(B.elements{qq}(pp,3), 2),1)]),...
                             mean([nodes.all(nodes.nind(B.elements{qq}(pp,1), 2),2), nodes.all(nodes.nind(B.elements{qq}(pp,2), 2),2), nodes.all(nodes.nind(B.elements{qq}(pp,3), 2),2)]),...
                             mean([nodes.all(nodes.nind(B.elements{qq}(pp,1), 2),3), nodes.all(nodes.nind(B.elements{qq}(pp,2), 2),3), nodes.all(nodes.nind(B.elements{qq}(pp,3), 2),3)]),...
                            ['\Gamma_{' num2str(pp) '}'],'Color', cmap(qq,:));
                    otherwise disp('! WARNING: Invalid dimension for the domain.')
                end
            end
        end
    else
        quiver3(B.elcenterpoint(:, 1), B.elcenterpoint(:, 2), B.elcenterpoint(:, 3), B.elnormalvector(:, 1), B.elnormalvector(:, 2), B.elnormalvector(:, 3), 0.2, 'Color', 'm')
        tmp.outerBelements = B.elements{B.outerGamma};
        for qq = 1:size(B.elements{B.outerGamma}, 1)
            switch dim
                case 1 % 1D
                    disp('! WARNING: No 1D support yet.')
                case 2 % 2D
                    text(mean([nodes.all(nodes.nind(tmp.outerBelements(qq,1), 2),1), nodes.all(nodes.nind(tmp.outerBelements(qq,end), 2),1)]),...
                        mean([nodes.all(nodes.nind(tmp.outerBelements(qq,1), 2),2), nodes.all(nodes.nind(tmp.outerBelements(qq,end), 2),2)]),...
                        ['\Gamma_{' num2str(qq) '}'],'Color', 'm');
                case 3 % 3D
                    text(mean([nodes.all(nodes.nind(B.elements{1}(qq,1), 2),1), nodes.all(nodes.nind(B.elements{1}(qq,2), 2),1), nodes.all(nodes.nind(B.elements{1}(qq,3), 2),1)]),...
                         mean([nodes.all(nodes.nind(B.elements{1}(qq,1), 2),2), nodes.all(nodes.nind(B.elements{1}(qq,2), 2),2), nodes.all(nodes.nind(B.elements{1}(qq,3), 2),2)]),...
                         mean([nodes.all(nodes.nind(B.elements{1}(qq,1), 2),3), nodes.all(nodes.nind(B.elements{1}(qq,2), 2),3), nodes.all(nodes.nind(B.elements{1}(qq,3), 2),3)]),...
                        ['\Gamma_{' num2str(qq) '}'],'Color', 'm');
                otherwise disp('! WARNING: Invalid dimension for the domain.')
            end
        end
    end
end
hold off
view(v)
end