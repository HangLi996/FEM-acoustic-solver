function [nodes, elements, dim, N] = msh2mat(file)
%   Input:  file    ... file struct
%   Output: nodes   ... nodes struct
%           elements... elements struct
%           dim     ... dimension of the domain
%
%   Import of Gmsh .msh files
%   Gmsh 2.8 documentation:
%       http://geuz.org/gmsh/doc/texinfo/gmsh.html
%   Gmsh .msh-file structure
%       http://stderr.org/doc/gmsh/html/gmsh_10.html#SEC63
%
%   Written by Lennart Moheit (Lennart.Moheit@tum.de)
%   06/02/2016
%
    fileID = fopen(file.mesh);
    formatSpec = '%s';
    C = textscan(fileID, formatSpec, 'Delimiter', '\t');
    fclose(fileID);

    %% Check how many elements and nodes are given in the .msh-file
    for qq = 1:size(C{1}, 1)
        % Find nodal table
        if strcmp(cellstr(C{1}{qq}), '$Nodes')
            pos.n = qq + 1;   % keep row in mind
            N.n = str2double(C{1}{pos.n});
            if ~strcmp(cellstr(C{1}{pos.n + 1 + N.n}), '$EndNodes')
                disp('! WARNING: Something is wrong with the number of nodes in the .msh-file')
            end
        % Find elements table
        elseif strcmp(cellstr(C{1}{qq}), '$Elements')
            pos.el = qq + 1;   % keep row in mind
            N.el = str2double(C{1}{pos.el});
            if ~strcmp(cellstr(C{1}{pos.el + 1 + N.el}), '$EndElements')
                disp('! WARNING: Something is wrong with the number of elements in the .msh-file')
            end
        end
    end

    %% Get table of nodes
    nodes.nall = zeros(N.n, 1);
    nodes.all = zeros(N.n, 3);
    for qq = 1:N.n
        tmp.nodesnall = str2num(C{1}{pos.n + qq});
        tmp.nodesall = str2num(C{1}{pos.n + qq});
        nodes.nall(qq, :) = tmp.nodesnall(1);
        nodes.all(qq, :) = tmp.nodesall(2:4);
    end

    %% Get table of elements
    tmp.elements = zeros(N.el, 5+27); % maximum number of columns (5 columns + 27 nodes)
    tmp.eldim = [];
    elements.type = []; elements.order = []; 
    elements.vert = []; elements.all = [];

    % Classify elements according to http://stderr.org/doc/gmsh/html/gmsh_10.html#SEC63
    for qq = 1:N.el
        tmp.elements(qq, :) = [str2num(C{1}{pos.el + qq}), zeros(1, size(tmp.elements, 2) - size(str2num(C{1}{pos.el + qq}), 2))]; % fill empty columns with zeros
        switch tmp.elements(qq, 2)
            % 0D - Point
            case 15 % 1-node point

            % 1D - Line
            case 1 % 2-node line (for 1D elements)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 1];
                elements.shape = 'line';
            case 8 % 3-node second order line (2 nodes associated with the vertices and 1 with the edge)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 2];
                elements.shape = 'line';
                
            % 2D - Surface
            case 2 % 3-node triangle
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 1];
                elements.shape = 'tri';
            case 3 % 4-node quadrangle
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 1];
                elements.shape = 'quad';
            case 9 % 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 2];
                elements.shape = 'tri';
            case 10 % 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 2];
                elements.shape = 'quad';
            case 16 % 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 2];
                elements.shape = 'quad';
            case 21 % 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges and 1 with the face)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 3];
                elements.shape = 'tri';

            % 3D - Volume
            case 4 % 4-node tetrahedron
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 1];
            case 5 % 8-node hexahedron
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 1];
            case 6 % 6-node prism
                disp(['! WARNING: Element type "' num2str(tmp.elements(qq, 2)) '" not supported yet.'])
            case 7 % 5-node pyramid
                disp(['! WARNING: Element type "' num2str(tmp.elements(qq, 2)) '" not supported yet.'])         
            case 11 % 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 2];
            case 12 % 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 2];
            case 13 % 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)
                disp(['! WARNING: Element type "' num2str(tmp.elements(qq, 2)) '" not supported yet.'])
            case 14 % 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)
                disp(['! WARNING: Element type "' num2str(tmp.elements(qq, 2)) '" not supported yet.'])        
            case 17 % 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)
                elements.type = [elements.type; tmp.elements(qq, 2)];
                elements.all = [elements.all; tmp.elements(qq, 6:end)]; tmp.eldim = [tmp.eldim; getDimension(tmp.elements(qq, 2))];
                elements.order = [elements.order; 2];
            case 18 % 15-node second order prism (6 nodes associated with the vertices and 9 with the edges)
                disp(['! WARNING: Element type "' num2str(tmp.elements(qq, 2)) '" not supported yet.'])        
            case 19 % 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)
                disp(['! WARNING: Element type "' num2str(tmp.elements(qq, 2)) '" not supported yet.'])        
            otherwise disp('! WARNING: Unkown element type.')
        end
    end
    %
    tmp.verts ={{1 2};      % 2-node line.
                {1:3};      % 3-node triangle.
                {1:4};      % 4-node quadrangle. 
                {1:4};      % 4-node tetrahedron. 
                {1:8};      % 8-node hexahedron. 
                {'NA'};     % 6-node prism. 
                {'NA'};     % 5-node pyramid. 
                {1 2};      % 3-node second order line (2 nodes associated with the vertices and 1 with the edge). 
                {1:3};    % 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges). 
                {1:4};  % 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face). 
                {1:4};  % 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges). 
                {1:8}; % 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume). 
                {'NA'};     % 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces). 
                {'NA'};     % 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face). 
                {1};        % 1-node point. 
                {1:4};  % 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges). 
                {1:8}; % 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges). 
                {'NA'};     % 15-node second order prism (6 nodes associated with the vertices and 9 with the edges). 
                {'NA'}};    % 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges). 
    
    dim = getDimension(elements.type);
    elements.order = unique(elements.order); if length(elements.order) > 1; disp('! WARNING: Different element orders in msh2mat.'); end;
    % Remove boundary elements of dimension dim-1
    elements.all = elements.all(tmp.eldim == max(tmp.eldim), :);
    elements.type = elements.type(tmp.eldim == max(tmp.eldim));
    % Cut zero columns (for mixed meshs)
    elements.all = elements.all(:, any(elements.all, 1));    
    N.el = size(elements.all, 1);
    for qq = 1:N.el
        elements.vert(qq, :) = elements.all(qq, cell2mat(tmp.verts{elements.type(qq)}));
    end
    elements.corners = size(elements.vert, 2);
    
    %% Get rid of nodes/dof that are not part of the mesh
    tmp.meshnodes = unique(elements.all);
    nodes.all = nodes.all(tmp.meshnodes, :);
    nodes.nall = nodes.nall(tmp.meshnodes, :);
    N.n = length(tmp.meshnodes); % New total number of dof
    N.i = size(elements.all, 2);
    % Create a new table which gives the indices of the nodes of the mesh    
    nodes.nind = zeros(length(nodes.nall), 2);
    nodes.nind(nodes.nall, 1) = nodes.nall;    
    tmp.nonzeroind = find(ismember(nodes.nind(:, 1), nodes.nall));
    nodes.nind(tmp.nonzeroind, 2) = 1:length(nodes.nall);
    nodes.ninmesh = nodes.nind(find(nodes.nind(:, 2)), 2);
    % Create another elements.all table using straight indices from 1:N.n
    elements.nall = changem(elements.all, nodes.nind(:, 2), nodes.nind(:, 1));
end