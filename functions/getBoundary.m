function [B, N] = getBoundary(file, conf, N, nodes, elements)
% function [B, BEM, inner, N] = getBoundary(file, N, nodes, elements)
%   Input:  nodes           ... nodal matrix
%           elements        ... element matrix
%   Output: B.n             ... nodes at the boundary
%           B.nodes         ... nodal coordinates of boundary nodes
%           B.lines         ... boundary lines
%           B.elements      ... boundary elements (lines for 2D or surfaces for 3D)
%           B.type        ... element type for boundary elements
%           B.centerpoint   ... centerpoint of boundary lines and surfaces
%       analogously for inner nodes and elements
%           inner.n
%           inner.nodes
%           inner.lines
%           inner.surfaces
%           inner.elements
%
%   Find boundaries of the given domain by finding those lines or surfaces, that do not appear
%   more than once as for inner elements.
%
%   by Lennart Moheit 06/06/2016
% 检查file.path字段是否存在
if ~isfield(file, 'path')
    file.path = './';
end
if exist([file.path 'Boundary.mat'], 'file')
    load([file.path '/Boundary.mat']);
else
    B.order = elements.order;
    B.type = zeros(size(elements.type, 1), 1);    
    dim = getDimension(elements.type);
    switch dim
        %______________________________________________________________________
        %% 1D
        case 1 % 0D boundary point for 1D lines
            disp('! WARNING: No 1D support in getBoundary.')
        %%_________________________________________________________________
        %     ____  ____  
        %    |___ \|  _ \ 
        %      __) | | | |
        %     / __/| |_| |
        %    |_____|____/ 
        %%_________________________________________________________________
        case 2 % 1D boundary lines for 2D surfaces
            %% elements.type 2 - 3-node triangle
            lines = [];
            elindices = find(elements.type == 2);
            if elindices
                lines = [lines; ...
                    elements.all(elindices, 1:2); ...
                    elements.all(elindices, 2:3); ...
                    [elements.all(elindices, 3), elements.all(elindices, 1)]];
                B.type = 1; % 2-node line (for 1D elements)
            end
            %% elements.type 3 - 4-node quadrangle
            elindices = find(elements.type == 3);
            if elindices
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                lines = [lines; ...
                    elements.all(elindices, 1:2); ...
                    [elements.all(elindices, 2), elements.all(elindices, 4)]
                    [elements.all(elindices, 4), elements.all(elindices, 3)]
                    [elements.all(elindices, 3), elements.all(elindices, 1)]];
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                lines = [lines; ...
                    elements.all(elindices, 1:2); ...
                    elements.all(elindices, 2:3); ...
                    elements.all(elindices, 3:4); ...
                    [elements.all(elindices, 4), elements.all(elindices, 1)]];    
                end
                B.type = 1; % 2-node line (for 1D elements)
            end
            %% elements.type 9 - 6-node second order triangle (3 nodes associated with the vertices and 3 with the lines)
            elindices = find(elements.type == 9);
            if elindices
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                lines = [lines; ...
                    elements.all(elindices, 1:2); ...
                    elements.all(elindices, 2:3); ...
                    [elements.all(elindices, 3), elements.all(elindices, 5)]; ...
                    [elements.all(elindices, 5), elements.all(elindices, 6)]; ...
                    [elements.all(elindices, 6), elements.all(elindices, 4)]; ...
                    [elements.all(elindices, 4), elements.all(elindices, 1)]];
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                lines = [lines; ...
                    [elements.all(elindices, 1), elements.all(elindices, 4)]; ...
                    [elements.all(elindices, 4), elements.all(elindices, 2)]; ...
                    [elements.all(elindices, 2), elements.all(elindices, 5)]; ...
                    [elements.all(elindices, 5), elements.all(elindices, 3)]; ...
                    [elements.all(elindices, 3), elements.all(elindices, 6)]; ...
                    [elements.all(elindices, 6), elements.all(elindices, 1)]];
                end
                B.type = 8; % 3-node second order line (2 nodes associated with the vertices and 1 with the edge)
            end
            %% elements.type 10 - 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the lines and 1 with the face)
            elindices = find(elements.type == 10);
            if elindices
                % COMSOL
                if strcmp(conf.mesher, 'COMSOL')
                lines = [lines; ...
                    elements.all(elindices, 1:2); ...
                    elements.all(elindices, 2:3); ...
                    [elements.all(elindices, 3), elements.all(elindices, 6)]; ...
                    [elements.all(elindices, 6), elements.all(elindices, 9)]; ...
                    [elements.all(elindices, 9), elements.all(elindices, 8)]; ...
                    [elements.all(elindices, 8), elements.all(elindices, 7)]; ...
                    [elements.all(elindices, 7), elements.all(elindices, 4)]; ...
                    [elements.all(elindices, 4), elements.all(elindices, 1)]];
                % GMSH
                elseif strcmp(conf.mesher, 'GMSH')
                lines = [lines; ...
                    [elements.all(elindices, 1), elements.all(elindices, 5)]; ...
                    [elements.all(elindices, 5), elements.all(elindices, 2)]; ...
                    [elements.all(elindices, 2), elements.all(elindices, 6)]; ...
                    [elements.all(elindices, 6), elements.all(elindices, 3)]; ...
                    [elements.all(elindices, 3), elements.all(elindices, 7)]; ...
                    [elements.all(elindices, 7), elements.all(elindices, 4)]; ...
                    [elements.all(elindices, 4), elements.all(elindices, 8)]; ...
                    [elements.all(elindices, 8), elements.all(elindices, 1)]];
                end
                B.type = 8; % 3-node second order line (2 nodes associated with the vertices and 1 with the edge)
            end
            %% elements.type 16 - 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)
            elindices = find(elements.type == 16);
            if elindices
                lines = [lines; ...
                    elements.all(elindices, 1:2); ...
                    elements.all(elindices, 2:3); ...
                    [elements.all(elindices, 3), elements.all(elindices, 5)]; ...
                    [elements.all(elindices, 5), elements.all(elindices, 8)]; ...
                    [elements.all(elindices, 8), elements.all(elindices, 7)]; ...
                    [elements.all(elindices, 7), elements.all(elindices, 6)]; ...
                    [elements.all(elindices, 6), elements.all(elindices, 4)]; ...
                    [elements.all(elindices, 4), elements.all(elindices, 1)]];
                B.type = 8; % 3-node second order line (2 nodes associated with the vertices and 1 with the edge)
            end
            %% elements.type 21 - 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges and 1 with the face)
            elindices = find(elements.type == 21);
            if elindices
                lines = [lines; ...
                    elements.all(elindices, 1:2); ...
                    elements.all(elindices, 2:3); ...
                    elements.all(elindices, 3:4); ...
                    [elements.all(elindices, 4), elements.all(elindices, 7)]; ...
                    [elements.all(elindices, 7), elements.all(elindices, 9)]; ...
                    [elements.all(elindices, 9), elements.all(elindices, 10)]; ...
                    [elements.all(elindices, 10), elements.all(elindices, 8)]; ...
                    [elements.all(elindices, 8), elements.all(elindices, 5)]; ...
                    [elements.all(elindices, 5), elements.all(elindices, 1)]];
                B.type = NaN; % 4-node third order line (2 nodes associated with the vertices and 2 with the edge)
            end

            % Anything else?
            if any(ismember(elements.type, [2 3 9 10 16 21])==0)
                disp('! WARNING: At least one choosen element type is not implemented yet in getBoundary.')
            end

            %% Get unique boundary lines
            lines = sort(lines, 2);       
            [tmp.b, tmp.m, tmp.n] = unique(lines,'rows');
                % tmp.bounds=tmp.b(tmp.n,:); % = boundary
            % This is how it works:
                % tmp.m gives the indices of bounds, which are unique
                % find(ismember(~)==0) gives those indices of bounds, which are NOT unique
                % bounds(~) gives each two nodes (= one boundary), which are NOT unique
                % ismember(bounds,~,'rows') gives indices for row numbers with ...
                %   ... boundaries that are NOT unique
                % find(~==0) gives those rows that are ABSOLUTELY unique
                % bounds(~, :) gives a matrix of ABSOLUTELY unique boundaries
            tmp.Blines = lines(find(ismember(lines,lines(find(ismember(1:size(lines, 1), tmp.m)==0), :),'rows')==0), :);
                tmp.rawBlines = tmp.Blines;
            inner.lines = setdiff(lines, tmp.Blines, 'rows');
            inner.elements = 0; % No inner element surfaces possible in 2D domain
            
            %% INNER BOUNDARIES / HOLES
            dd = 0;
            % find complete boundary-threads
            tmp.Bnodes = unique(tmp.Blines);
            tmp.dominoesfound = false; % found all lines?
            while tmp.dominoesfound == false  
                dd = dd + 1;
                % Dominoes        
                tmp.dominoline = 1;
                tmp.dominostone = tmp.Blines(1, :);
                tmp.Blines = tmp.Blines(2:end, :);
                tmp.dominonew = tmp.dominostone(2);
                tmp.dominothread = [tmp.dominostone];
                B.Gamma{dd} = [tmp.dominostone];
                tmp.dominoloop = false;
                while tmp.dominoloop == false
                    [tmp.row, ~] = find(ismember(tmp.rawBlines, tmp.dominonew));               
                    tmp.dominoneighbours = tmp.rawBlines(tmp.row, :);
                    if size(tmp.dominoneighbours, 1) ~= 2
                        disp('! WARNING: At least one boundary node has less or more than two neighbours.')
                    end
                    tmp.dominostone = tmp.dominoneighbours(find(~ismember(tmp.dominoneighbours, tmp.dominostone, 'rows')), :);
                        B.Gamma{dd} = [B.Gamma{dd}; tmp.dominostone];
                    tmp.dominonew = tmp.dominostone(find(tmp.dominostone ~= tmp.dominonew));
                    if tmp.dominonew == tmp.dominothread(1)
                        tmp.dominoloop = true;
                    else
                        tmp.dominothread = [tmp.dominothread, tmp.dominonew];
                    end     
                    tmp.Blines = tmp.Blines(find(~ismember(tmp.Blines, tmp.dominostone, 'rows')), :);
                end
                B.domino{dd} = tmp.dominothread;
                %
                if size(tmp.Blines, 1) == 0
                    tmp.dominoesfound = true;
                elseif size(tmp.Blines, 1) == 2
                    disp('! WARNING: Invalid number of remaining, unassigned boundary nodes in getBoundary.')
                end            
            end
            for dd = 1:size(B.domino, 2)
                [~, ~, tmp.r{dd}] = cart2sph(nodes.all(nodes.nind(double(B.domino{dd}), 2), 1), nodes.all(nodes.nind(double(B.domino{dd}), 2), 2), nodes.all(nodes.nind(double(B.domino{dd}), 2), 3));
                tmp.r0(dd) = mean(double(tmp.r{dd}));
            end
            % Find domino line with largest radius -> this is assumed to be the
            % outer boundary line
                B.outerGamma = find(tmp.r0==max(tmp.r0));
                B.innerGamma = setdiff(1:size(B.Gamma, 2), B.outerGamma);
            B.lines{B.outerGamma} = double(B.Gamma{B.outerGamma});
            for qq = 1:length(B.innerGamma)
                B.lines{B.innerGamma(qq)} = double(B.Gamma{B.innerGamma(qq)});
            end
            
            %% BOUNDARY AND INNER NODES
            B.n{B.outerGamma} = unique(B.lines{B.outerGamma});
            B.nodes{B.outerGamma} = nodes.all(nodes.nind(B.n{B.outerGamma}, 2), :);
                B.ninner = [];
            for qq = 1:length(B.innerGamma)
                B.n{B.innerGamma(qq)} = unique(B.lines{B.innerGamma(qq)});
                B.nodes{B.innerGamma(qq)} = nodes.all(nodes.nind(B.n{B.innerGamma(qq)}, 2), :);
                B.ninner = [B.ninner; B.n{B.innerGamma(qq)}];
            end
            B.ninner = sort(B.ninner);
            B.nall = [B.n{B.outerGamma}; B.ninner];
            B.nall = unique(B.nall);
                inner.n = setdiff(nodes.nall, B.nall)';
            B.nodesall = nodes.all(nodes.nind(B.nall, 2), :);
                inner.nodes = nodes.all(nodes.nind(inner.n, 2), :);
                
            %% BOUNDARY ELEMENTS AND VERTS
            if B.order == 1 % linear FE elements -> linear boundary elements
                B.elements{B.outerGamma} = B.lines{B.outerGamma};                
                for qq = 1:length(B.innerGamma)
                    B.elements{B.innerGamma(qq)} = B.lines{B.innerGamma(qq)};
                end
            elseif B.order == 2 % quad FE elements
                tmp.Blines = B.domino{B.outerGamma}.';
                tmp.Blines = kron(tmp.Blines, [1; 0]);
                tmp.Blines(6:4:length(tmp.Blines)) = tmp.Blines(5:4:length(tmp.Blines));
                tmp.Blines = tmp.Blines(tmp.Blines~=0);
                tmp.Blines(end+1) = tmp.Blines(1);
                tmp.Blines = reshape(tmp.Blines, 3, size(tmp.Blines, 1)/3).'; % get quadratic boundary lines (3 nodes)
                B.elements{B.outerGamma} = tmp.Blines;              
                for qq = 1:length(B.innerGamma)
                    tmp.Blines = B.domino{B.innerGamma(qq)}.';
                    tmp.Blines = kron(tmp.Blines, [1; 0]);
                    tmp.Blines(6:4:length(tmp.Blines)) = tmp.Blines(5:4:length(tmp.Blines));
                    tmp.Blines = tmp.Blines(tmp.Blines~=0);
                    tmp.Blines(end+1) = tmp.Blines(1);
                    tmp.Blines = reshape(tmp.Blines, 3, size(tmp.Blines, 1)/3).'; % get quadratic boundary lines (3 nodes)
                    B.elements{B.innerGamma(qq)} = tmp.Blines;
                end
            else disp('! WARNING: Invalid order of boundary elements.')
            end  
            % ... to get consistent order
            B.vert{B.outerGamma} = B.elements{B.outerGamma}(:, [1 end]);
            for qq = 1:length(B.innerGamma)
                B.vert{B.innerGamma(qq)} = B.elements{B.innerGamma(qq)}(:, [1 end]);
            end
            B.corners = size(B.vert{B.outerGamma}, 2);
%                 B.nvert = unique(B.vert);
%                 N.Bnvert = length(B.nvert);  

            %% CENTERPOINTS AND NORMALVECTORS
            tmp.outerBvert = B.vert{B.outerGamma};
            tmp.outerBelements = B.elements{B.outerGamma};
            inner.leftelementnodes = zeros(size(B.vert, 1), size(elements.all, 2));
            for el = 1:size(tmp.outerBvert, 1)
                % Find centerpoint of boundary lines
                B.elcenterpoint(el, :) = [mean([nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), 1), nodes.all(nodes.nind(tmp.outerBvert(el, 2), 2), 1)]), ...
                    mean([nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), 2), nodes.all(nodes.nind(tmp.outerBvert(el, 2), 2), 2)]), ...
                    mean([nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), 3), nodes.all(nodes.nind(tmp.outerBvert(el, 2), 2), 3)])];

                % Get correct outward normal vector and rearrange nodal orientation per element
                    % Find the respective complete surface element for the boundary line element
                    tmp.Bel(el) = find(sum(ismember(elements.all, tmp.outerBvert(el, :)), 2) == 2); % ...==2 for lines
                    % Respective left inner nodes for boundary elements
                    tmp.leftelementnodes(el, :) = setdiff(elements.all(tmp.Bel(el), :), tmp.outerBvert(el, :));
                        inner.leftelementnodes(el, :) = [setdiff(tmp.leftelementnodes(el, :), B.nall), zeros(1, size(elements.all, 2) - length(setdiff(tmp.leftelementnodes(el, :), B.nall)))];
                    %
                    tmp.boundaryline(el, :) = nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), :) - nodes.all(nodes.nind(tmp.outerBvert(el, 2), 2), :);
                    % normalvector = [x1, -y1] / sqrt(x1^2 + y1^2); while z1 keeps constant for 2D
                        tmp.normalvectorscale = 20 * sqrt(tmp.boundaryline(el,1)^2 + tmp.boundaryline(el,2)^2);
                    B.elnormalvector(el, :) = [[tmp.boundaryline(el,2), -tmp.boundaryline(el,1)] ./ tmp.normalvectorscale, tmp.boundaryline(el,3)];
                        % ... just another vector from the inner node of the element to
                        % the first node on the boundary element           
                        if isempty(inner.leftelementnodes) || sum(any(inner.leftelementnodes))==0
                            %tmp.Ovec = [0 0 0];
                            tmp.Ovec = mean(nodes.all(elements.nall(tmp.Bel(el), :), :));
                            tmp.vectorp = nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), :) - tmp.Ovec;
                        else
                            % Old version has some problems if there are
                            % no inner.leftelementnodes (e.g. for triangles
                            % with all element nodes on the boundary surface)
                            %   tmp.vectorp = nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), :) - mean(nodes.all(nodes.nind(inner.leftelementnodes(el, any(inner.leftelementnodes(el, :), 1)), 2), :), 1);
                            % New version taking the mean value of all
                            % element nodes as inner point:
                            tmp.vectorp = nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), :) - mean(nodes.all(nodes.nind(elements.all(tmp.Bel(el), :), 2), :), 1);
                        end                        
                
                    % negative scalar product for vectorp .* normalvector for wrong orientated elements
                    if sum(B.elnormalvector(el, :) .* tmp.vectorp) < 0 % rearrange the nodal order
%                         disp(['WRONG orientation for B.el# ' num2str(el)])
                        tmp.outerBvert(el, :) = flip(tmp.outerBvert(el, :));
                        tmp.outerBelements(el, :) = flip(tmp.outerBelements(el, :));
                    else % this one is fine, no need to worry!
%                         disp(['RIGHT orientation for B.el# ' num2str(el)])
                    end                    
                % Now calculate new normalvectors after rearrangement of boundary elements
                tmp.boundaryline(el,:) = nodes.all(nodes.nind(tmp.outerBvert(el, 1), 2), :) - nodes.all(nodes.nind(tmp.outerBvert(el, 2), 2), :);
                % normalvector = [x1, -y1] / sqrt(x1^2 + y1^2); while z1 keeps constant for 2D
                    tmp.normalvectorscale = 20 * sqrt(tmp.boundaryline(el,1)^2 + tmp.boundaryline(el,2)^2);
                B.elnormalvector(el, :) = [[tmp.boundaryline(el,2), -tmp.boundaryline(el,1)] ./ tmp.normalvectorscale, tmp.boundaryline(el,3)];                
            end
            B.vert{B.outerGamma} = tmp.outerBvert;
            B.elements{B.outerGamma} = tmp.outerBelements;
            B.higherordernodes = setdiff(B.elements{B.outerGamma}, B.vert{B.outerGamma}); % higher order nodes on boundary
            
            %% NEIGHBOR/ADJACENT ELEMENTS
            N.ofboundaries = numel(B.elements); % Number of closed boundaries including inner and outer boundaries
            qq = 0;
                % Nodes that are already in use for a certain element. 
                % If the first element has nodes 1 and 2, and the second
                % element has nodes 2 and 3, then element number 2 gets an
                % additional node index for the second node, N.n+1, because it 
                % is already in use and we need a new index for the vector of
                % nodal velocities (mutual boundary nodes)
                tmp.n_all = [];
            for bb = 1:N.ofboundaries % for all boundaries ...
                tmp.Bel = B.elements{bb};                
                for el = 1:size(tmp.Bel, 1) % for each boundary element of the current closed boundary (bb)
                    elindices = tmp.Bel(el, :);                    
                    % For each closed boundary line we walk through all its
                    % elements el, that's why we use qq as boundary element index
%                     qq = qq + 1;                    
                    %______________________________________________________________
                    % Whose neighbour element number is the current element qq
                    if el == 1 % fist boundary element
                        B.adjacent{bb}(el).el = [size(tmp.Bel, 1); ...
                                             2];
                    elseif el == size(tmp.Bel, 1) % last boundary element
                        B.adjacent{bb}(el).el = [size(tmp.Bel, 1) - 1; ...
                                             1];
                    else
                        B.adjacent{bb}(el).el = [el - 1; ...
                                             el + 1];
                    end
                    
                    % Adjacent elements
                    B.adjacent{bb}(el).elements = tmp.Bel(B.adjacent{bb}(el).el, :);
                    % Find those adjacent elements, whose element surface normal vectors are not parallel
                    B.adjacent{bb}(el).nonparallel = any(cross(B.elnormalvector(B.adjacent{bb}(el).el, :), repmat(B.elnormalvector(el, :), [size(B.adjacent{bb}(el).el, 1), 1])), 2);                 
                    % Find mutual nodes for current and adjacent elements
                        [tmp.nind, ind] = ismember(B.adjacent{bb}(el).elements, tmp.Bel(el, :));
                        [~, ind] = sort(ind(tmp.nind));
                        tmp.nind = find(tmp.nind);
                    B.adjacent{bb}(el).nind = tmp.nind(ind);
                    B.adjacent{bb}(el).n = B.adjacent{bb}(el).elements(B.adjacent{bb}(el).nind);
                    
%% SINNLOS?
%                     % Which shape function Ni have the mutual nodes at adjacent boundary elements?
%                         [tmp.Niind, ind] = sort(B.adjacent{bb}(el).nind);
%                         tmp.Niind = changem(1:length(tmp.Niind), tmp.Niind);
%                     B.adjacent{bb}(el).Niind = tmp.Niind(ind);
%%
                    
                    % Additional lines have to be added to the nodal vectors of
                    % boundary admittance and velocity
                    if qq == 1
                        B.adjacent{bb}(1).Yvn = [elindices;
                                             elindices];
                    else                            
                        B.adjacent{bb}(el).Yvn = zeros(2, length(elindices));
                        tmp.n_new = ~ismember(elindices, [B.adjacent{bb}(:).Yvn]);                            
                        tmp.Yvn = zeros(1, length(elindices));
                        tmp.Yvn(tmp.n_new) = elindices(tmp.n_new);
                        if any(~tmp.n_new)
                            tmp.Yvn(~tmp.n_new) = max([nodes.nall; unique([B.adjacent{bb}(:).Yvn])]) + (1:length(find(~tmp.n_new)));
                        end
                        B.adjacent{bb}(el).Yvn = [elindices; tmp.Yvn];
                    end                    
                end
            end
            
            B.shape = 'line';
            N.iB = size([B.elements{B.outerGamma}], 2); % Number of nodes per boundary element
        %%_________________________________________________________________
        %     _____ ____  
        %    |___ /|  _ \ 
        %      |_ \| | | |
        %     ___) | |_| |
        %    |____/|____/ 
        %%_________________________________________________________________   
        case 3 % 2D boundary surfaces for 3D volumes           
            lines = [];
            surfaces = [];
            for el = 1:size(elements.all, 1)
                elnodes = elements.all(el, :); % nodes of the current element el
                % Get lines for each element (ellines)
                switch elements.type(el)
                    case 4 % 4-node tetrahedron
                        if strcmp(conf.mesher, 'COMSOL')
                            % TODO same as gmsh?
                            disp('! WARNING: COMSOL not implemented yet.')
                        elseif strcmp(conf.mesher, 'GMSH')
                            ellines = [elnodes(1:2); ...
                                       elnodes(2:3); ...
                                       [elnodes(3), elnodes(1)]; ...
                                       [elnodes(2), elnodes(4)]; ...
                                       [elnodes(3), elnodes(4)]; ...
                                       [elnodes(4), elnodes(1)]];
                            elsurfaces = [elnodes(1:3); ...
                                         [elnodes(1), elnodes(2), elnodes(4)]; ...
                                         [elnodes(1), elnodes(3), elnodes(4)]; ...
                                         [elnodes(2), elnodes(3), elnodes(4)]];
                        end
                        B.type = 2; % 3-node triangle
                        B.shape = 'tri';
                    case 5 % 8-node hexahedron
                        % COMSOL
                        if strcmp(conf.mesher, 'COMSOL')
                        ellines = [elnodes(1:2); ...
                                   elnodes(2:3); ...
                                   elnodes(3:4); ...              
                                   [elnodes(4), elnodes(1)]; ...
                                   elnodes(5:6); ...
                                   elnodes(6:7); ...
                                   elnodes(7:8); ...              
                                   [elnodes(8), elnodes(5)]; ...
                                   [elnodes(1), elnodes(5)]; ...
                                   [elnodes(2), elnodes(6)]; ...
                                   [elnodes(3), elnodes(7)]; ...              
                                   [elnodes(4), elnodes(8)]];
                        % GMSH
                        elseif strcmp(conf.mesher, 'GMSH')
                        ellines = [elnodes(1:2); ...
                                   elnodes(2:3); ...
                                   elnodes(3:4); ...              
                                   [elnodes(4), elnodes(1)]; ...
                                   elnodes(5:6); ...
                                   elnodes(6:7); ...
                                   elnodes(7:8); ...              
                                   [elnodes(8), elnodes(5)]; ...
                                   [elnodes(1), elnodes(5)]; ...
                                   [elnodes(2), elnodes(6)]; ...
                                   [elnodes(3), elnodes(7)]; ...              
                                   [elnodes(4), elnodes(8)]];
                        elsurfaces = [elnodes(4:-1:1); ...
                                      elnodes(5:8); ...
                                      [elnodes(1), elnodes(2), elnodes(6), elnodes(5)]; ...
                                      [elnodes(1), elnodes(5), elnodes(8), elnodes(4)]; ...
                                      [elnodes(4), elnodes(8), elnodes(7), elnodes(3)]; ...
                                      [elnodes(2), elnodes(3), elnodes(7), elnodes(6)]];
                        end
                        B.type = 3; % 4-node quadrangle
                        B.shape = 'quad';
                    case 11 % 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the lines)
                        if strcmp(conf.mesher, 'COMSOL')
                            % TODO same as gmsh?
                            disp('! WARNING: COMSOL not implemented yet.')
                        elseif strcmp(conf.mesher, 'GMSH')
                            ellines = [[elnodes(1),elnodes(7)]; ...
                                       [elnodes(7),elnodes(3)]; ...
                                       [elnodes(3),elnodes(6)]; ...
                                       [elnodes(6),elnodes(2)]; ...
                                       [elnodes(2),elnodes(5)]; ...
                                       [elnodes(5),elnodes(1)]; ...
                                       [elnodes(1),elnodes(8)]; ...
                                       [elnodes(8),elnodes(4)]; ...
                                       [elnodes(4),elnodes(9)]; ...
                                       [elnodes(9),elnodes(3)]; ...
                                       [elnodes(4),elnodes(10)]; ...
                                       [elnodes(10),elnodes(2)]];
                            elsurfaces = [[elnodes(1:3) elnodes(5:7)]; ...
                                          [elnodes(2) elnodes(4) elnodes(3) elnodes(10) elnodes(9) elnodes(6)]; ...
                                          [elnodes(1) elnodes(4) elnodes(2) elnodes(8) elnodes(10) elnodes(5)]; ...
                                          [elnodes(1) elnodes(3) elnodes(4) elnodes(7) elnodes(9) elnodes(8)]];
                        end
                        B.type = 9; % 6-node triangle
                        B.shape = 'tri';
                    case 12 % 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the lines, 6 with the faces and 1 with the volume)
                        if strcmp(conf.mesher, 'COMSOL')
                            % TODO same as gmsh?
                            disp('! WARNING: COMSOL not implemented yet.')
                        elseif strcmp(conf.mesher, 'GMSH')
                            ellines = [[elnodes(1),elnodes(10)]; ...
                                       [elnodes(10),elnodes(4)]; ...
                                       [elnodes(4),elnodes(14)]; ...
                                       [elnodes(14),elnodes(3)]; ...
                                       [elnodes(3),elnodes(12)]; ...
                                       [elnodes(12),elnodes(2)]; ...
                                       [elnodes(2),elnodes(9)]; ...
                                       [elnodes(9),elnodes(1)]; ...
                                       %
                                       [elnodes(1),elnodes(11)]; ...
                                       [elnodes(11),elnodes(5)]; ...
                                       [elnodes(5),elnodes(10)]; ...
                                       [elnodes(10),elnodes(8)]; ...
                                       [elnodes(8),elnodes(16)]; ...
                                       [elnodes(16),elnodes(4)]; ...
                                       %
                                       [elnodes(3),elnodes(15)]; ...
                                       [elnodes(15),elnodes(7)]; ...
                                       [elnodes(7),elnodes(19)]; ...
                                       [elnodes(19),elnodes(6)]; ...
                                       [elnodes(6),elnodes(13)]; ...
                                       [elnodes(13),elnodes(2)]; ...
                                       %
                                       [elnodes(5),elnodes(17)]; ...
                                       [elnodes(17),elnodes(6)]; ...
                                       [elnodes(7),elnodes(20)]; ...
                                       [elnodes(20),elnodes(8)]];
                            elsurfaces = [[elnodes(4:-1:1), elnodes(14) elnodes(12), elnodes(9), elnodes(10), elnodes(21)]; ...
                                          [elnodes(5:8) elnodes(17) elnodes(19:20) elnodes(18) elnodes(26)]; ...
                                          [elnodes(1) elnodes(2) elnodes(6) elnodes(5) ...
                                             elnodes(9) elnodes(13) elnodes(17) elnodes(11) elnodes(22)]; ...
                                          [elnodes(1), elnodes(4), elnodes(8), elnodes(5), ...
                                             elnodes(10) elnodes(16) elnodes(18) elnodes(11) elnodes(23)]; ...
                                          [elnodes(4), elnodes(8), elnodes(7), elnodes(3), ...
                                             elnodes(16) elnodes(20) elnodes(15) elnodes(14) elnodes(25)]; ...
                                          [elnodes(2), elnodes(3), elnodes(7), elnodes(6), ...
                                             elnodes(12) elnodes(15) elnodes(19) elnodes(13) elnodes(24)]];
                        end
                        B.type = 10; % 9-node quadrangle
                        B.shape = 'quad';
                    case 17 % 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the lines)
                        disp('! WARNING: At least one choosen element type is not implemented yet in getBoundary.')  
                    otherwise
                        disp('! WARNING: At least one choosen element type is not implemented yet in getBoundary.')
                end
                
                % specify number of corners
                switch elements.shape
                    case 'tri'
                        B.corners = 3;
                    case 'quad'
                        B.corners = 4;
                end
                
                % reduce to unique surfaces (not necessairily needed here)
                elsurfaces = getUniqueSurfaces(elsurfaces);
                surfaces = [surfaces; elsurfaces];
                
                % lines...
                lines = [lines; ellines];
            end % end for el
            
            % extract boundary surfaces
            tmp.Bsurfaces = getBoundarySurfaces(surfaces);
            % extract inner surfaces
            inner.surfaces = setdiff(surfaces, tmp.Bsurfaces, 'rows');
            
            %% DETECT BOUNDARIES
            % create a cell in B.Gamma for every boundary
            workdone = false;
            currSet = 1;
            B.Gamma{currSet} = [];
            while workdone == false
                B.Gamma{currSet} = tmp.Bsurfaces(1,:);
                tmp.Bsurfaces(1,:) = [];
                setdone = false;
                while setdone == false
                    [tmp.row, ~] = find(ismember(tmp.Bsurfaces, B.Gamma{currSet}));
                    tmp.row = unique(tmp.row);
                    B.Gamma{currSet} = [B.Gamma{currSet}; tmp.Bsurfaces(tmp.row,:)];
                    tmp.Bsurfaces(tmp.row,:) = [];
                    if isempty(tmp.row)
                        setdone = true;
                        continue;
                    end
                end
                if size(tmp.Bsurfaces, 1) == 0
                    workdone = true;
                else
                    currSet = currSet + 1;
                end
            end
            
            % calculate boundaries' mean distance to origin
            for currSet = 1:size(B.Gamma, 2)
                [~, ~, tmp.r{currSet}] = cart2sph(nodes.all(nodes.nind(double(B.Gamma{currSet}), 2), 1), nodes.all(nodes.nind(double(B.Gamma{currSet}), 2), 2), nodes.all(nodes.nind(double(B.Gamma{currSet}), 2), 3));
                tmp.r0(currSet) = mean(tmp.r{currSet});
            end
            
            % boundary with largest mean distance is assumed to be the
            % outer boundary
            B.outerGamma = find(tmp.r0==max(tmp.r0));
            B.innerGamma = setdiff(1:size(B.Gamma, 2), B.outerGamma);
            B.surfaces{B.outerGamma} = B.Gamma{B.outerGamma};
            for qq = 1:length(B.innerGamma)
                B.surfaces{B.innerGamma(qq)} = B.Gamma{B.innerGamma(qq)};
            end
            
            %% BOUNDARY AND INNER BOUNDARY NODES
            B.n{B.outerGamma} = unique(B.surfaces{B.outerGamma});
            B.nodes{B.outerGamma} = nodes.all(nodes.nind(B.n{B.outerGamma}, 2), :);
            B.ninner = [];
            for qq = 1:length(B.innerGamma)
                B.n{B.innerGamma(qq)} = unique(B.surfaces{B.innerGamma(qq)});
                B.nodes{B.innerGamma(qq)} = nodes.all(nodes.nind(B.n{B.innerGamma(qq)}, 2), :);
                B.ninner = [B.ninner; B.n{B.innerGamma(qq)}];
            end
            B.ninner = sort(B.ninner);
            B.nall = [B.n{B.outerGamma}; B.ninner];
            B.nall = unique(B.nall);
            inner.n = setdiff(nodes.nall, B.nall)';
            B.nodesall = nodes.all(nodes.nind(B.nall, 2), :);
            inner.nodes = nodes.all(nodes.nind(inner.n, 2), :);
                     
            %% BOUNDARY ELEMENTS AND VERTS
            B.elements{B.outerGamma} = B.surfaces{B.outerGamma};                
            for qq = 1:length(B.innerGamma)
                B.elements{B.innerGamma(qq)} = B.surfaces{B.innerGamma(qq)};
            end
            
            switch B.corners
                case 3
                    B.vert{B.outerGamma} = B.elements{B.outerGamma}(:,1:3);
                    for qq = 1:length(B.innerGamma)
                        B.vert{B.innerGamma(qq)} = B.elements{B.innerGamma(qq)}(:,1:3);
                    end
                case 4
                    B.vert{B.outerGamma} = B.elements{B.outerGamma}(1:4,:);
                    for qq = 1:length(B.innerGamma)
                        B.vert{B.innerGamma(qq)} = B.elements{B.innerGamma(qq)}(:,1:4);
                    end
            end

            %% CENTERPOINTS AND NORMALVECTORS
            for qq = 1:size(B.elements, 2)
                tmp.Bel = [];
                B.elnormalvector{qq} = [];
                B.elcenterpoint{qq} = [];
                for el = 1:size(B.elements{qq}, 1)
                    % Find element central point
                    B.elcenterpoint{qq} = [B.elcenterpoint{qq}; mean(nodes.all(nodes.nind(B.elements{qq}(el,:), 2), :))];
                    % Find the respective complete volume element for the boundary surface element
                    tmp.Bel = find(sum(ismember(elements.all, B.elements{qq}(el, :)), 2) == size(B.elements{qq}(el,:),2));
                    % Cross product n = a x b
                    tmp.vectora = nodes.all(nodes.nind(B.elements{qq}(el, 1), 2), :) - nodes.all(nodes.nind(B.elements{qq}(el, 2), 2), :);
                    tmp.vectorb = nodes.all(nodes.nind(B.elements{qq}(el, 3), 2), :) - nodes.all(nodes.nind(B.elements{qq}(el, 2), 2), :);
                    normalvector = cross(tmp.vectora, tmp.vectorb);
                    % vector pointing into the volume
                    tmp.vectorp = B.elcenterpoint{qq}(end,:) - mean(nodes.all(nodes.nind(elements.all(tmp.Bel, :), 2), :), 1);

                    % negative scalar product for vectorp .* normalvector 
                    % indicates wrong orientated elements
                    if dot(normalvector,tmp.vectorp) < 0 % rearrange the nodal order
                        normalvector = -1.*normalvector;
                    end
                    B.elnormalvector{qq}(el,:) = normalvector;
                end
            end
            
            %% NEIGHBOR/ADJACENT ELEMENTS
            currauxdof = max(B.nall)+1;
            for qq = 1:numel(B.elements) % boundary loop
                for ele = 1:size(B.elements{qq},1) % element loop
                    % find all nodes which connect to adjacent elements
                    B.adjacent{qq}(ele).n = unique(B.elements{qq}(ismember(B.elements{qq},B.elements{qq}(ele,:))));
                    % find adjacent elements
                    adjEleInd = find(sum(ismember(B.elements{qq},B.adjacent{qq}(ele).n'),2));
                    adjEleInd(adjEleInd == ele) = []; % remove current element from the list
                    B.adjacent{qq}(ele).elements = double(B.elements{qq}(adjEleInd,:));
                    B.adjacent{qq}(ele).el = double(adjEleInd);
                    % DEBUG:
                    % plotSurfaces(B.adjacent{qq}(ele).elements,nodes);
                    % coords = nodes.all(nodes.nind(B.elements{qq}(ele,:),2),:);
                    % plot3(coords(:,1), coords(:,2), coords(:,3),'ko','MarkerSize',10);
                    % 
                    % plotSurfaces(B.elements{qq}(B.adjacent{qq}(ele).el,:),nodes);
                    % coords = nodes.all(nodes.nind(B.elements{qq}(ele,:),2),:);
                    % plot3(coords(:,1), coords(:,2), coords(:,3),'ko','MarkerSize',10);
                    
                    % check for nonparallel elnormalvector in adjacent elements
                    B.adjacent{qq}(ele).nonparallel = ...
                        any(cross(B.elnormalvector{qq}(B.adjacent{qq}(ele).el, :), repmat(B.elnormalvector{qq}(ele, :), [size(B.adjacent{qq}(ele).el, 1), 1])), 2); 
                    % DEBUG:
                    % if sum(B.adjacent{qq}(ele).nonparallel) ~= 0
                    %     plotSurfaces(B.adjacent{qq}(ele).elements,nodes);
                    %     coords = nodes.all(nodes.nind(B.elements{qq}(ele,:),2),:);
                    %     plot3(coords(:,1), coords(:,2), coords(:,3),'ko','MarkerSize',10);
                    %     quiver3( ...
                    %         B.elcenterpoint{qq}(B.adjacent{qq}(ele).el(B.adjacent{qq}(ele).nonparallel), 1), ...
                    %         B.elcenterpoint{qq}(B.adjacent{qq}(ele).el(B.adjacent{qq}(ele).nonparallel), 2), ...
                    %         B.elcenterpoint{qq}(B.adjacent{qq}(ele).el(B.adjacent{qq}(ele).nonparallel), 3), ...
                    %         B.elnormalvector{qq}(B.adjacent{qq}(ele).el(B.adjacent{qq}(ele).nonparallel), 1), ...
                    %         B.elnormalvector{qq}(B.adjacent{qq}(ele).el(B.adjacent{qq}(ele).nonparallel), 2), ...
                    %         B.elnormalvector{qq}(B.adjacent{qq}(ele).el(B.adjacent{qq}(ele).nonparallel), 3), ...
                    %         0.2, 'Color', 'm')
                    % end
                    
                    % create auxiliary dofs for nodes which connect to adjacent elements
                    if ele == 1
                        B.adjacent{qq}(ele).Yvn = [B.elements{qq}(ele,:); ...
                                                   B.elements{qq}(ele,:)];
                    else
                        % check which dofs have to be replaced
                        auxdofs = double(ismember(B.elements{qq}(ele,:), [B.adjacent{qq}(:).Yvn]));
                        B.adjacent{qq}(ele).Yvn = [B.elements{qq}(ele,:); ...
                                                   B.elements{qq}(ele,:).*~auxdofs];
                        B.adjacent{qq}(ele).Yvn(2,auxdofs>0) = ...
                            (1:length(find(auxdofs))) + currauxdof;
                        currauxdof = currauxdof + sum(auxdofs);
                        % DEBUG:
                        % disp(B.adjacent{qq}(ele).Yvn);
                        % disp(currauxdof);
                    end
                end % element loop
            end % boundary loop
            %% COORDINATES
            for qq = 1:size(B.elements,2)
                B.nodes{qq} = nodes.all(nodes.nind(B.n{qq}, 2),:);
                inner.nodes = nodes.all(nodes.nind(inner.n, 2), :);
            end
            
        otherwise disp('! WARNING: Invalid dimension of the domain.')
    end % endswitch
    %__________________________________________________________________________
    
    % B.type extend to vector
    B.type = B.type * ones(size(B.elements{B.outerGamma}, 1), 1);
    % Numbers
    N.Bn = size(B.n{B.outerGamma}, 1);
    N.Bel = size(B.elements{B.outerGamma}, 1);
    % BEM transformation
    BEM.n = []; BEM.elements = [];
    for qq = 1:length(B.innerGamma)
        BEM.n = [BEM.n; B.n{B.innerGamma(qq)}];
        BEM.elements = [BEM.elements; B.elements{B.innerGamma(qq)}];
    end
    BEM.n = unique(BEM.n);
    BEM.nodes = nodes.all(nodes.nind(BEM.n, 2), :);
    BEM.elements = changem(BEM.elements, 1:length(BEM.n), BEM.n);
    
    %% SAVE AS .MAT-FILE
    %save([file.path 'Boundary.mat'], 'B', 'BEM', 'inner', 'N');   
end % endifexist
switch getDimension(elements.type)
    case 2
        disp(['>     Number of outer boundary elements (lines) ' num2str(size(B.elements{B.outerGamma}, 1))])
    case 3
        disp(['>     Number of outer boundary elements (surfaces) ' num2str(size(B.elements{B.outerGamma}, 1))])
    otherwise disp('! WARNING: Invalid dimension in getBoundary.')
end
end % endfunction

function surfaces = getUniqueSurfaces(surfaces)
    % clears non unique rows without sorting
    tmpsurfaces = sort(surfaces, 2);                                 % sort in ascending order
    [~, uniqueInd, ~] = unique(tmpsurfaces, 'rows');                 % get unique row indices
    nonuniqueInd = ~ismember(1:size(surfaces,1),uniqueInd);          % get non-unique row indices
    surfaces(nonuniqueInd,:) = [];                                   % clear non-unique rows
end

function surfaces = getBoundarySurfaces(surfaces)
    % extracts boundary surfaces (those which only occur once in 'surfaces')
    tmpsurfaces = sort(surfaces, 2);                                 % sort in ascending order
    [~, uniqueInd, ~] = unique(tmpsurfaces, 'rows');                 % get unique row indices
    nonuniqueInd = ~ismember(1:size(surfaces,1),uniqueInd);          % get non-unique row indices
    % get indices of every non-unique row (-> those of inner surfaces)
    innerInd = ismember(tmpsurfaces, tmpsurfaces(nonuniqueInd,:),'rows');
    surfaces(innerInd,:) = [];                                       % clear inner surfaces
end

function plotSurfaces(surfaces,nodes)
    % plots the given 'surfaces' using the general 'nodes' cell array.
    % mainly for debugging...
    %
    q = figure();
    hold on;
    surfaceNodes = nodes.all(nodes.nind(surfaces(:,:), 2),:);
    if size(surfaces,2) == 3
        plot3(surfaceNodes(:,1), surfaceNodes(:,2), surfaceNodes(:,3),'kx');
        linecolors = parula(size(surfaces,1));
        for el = 1:size(surfaces,1)
            surfaceNodes2 = nodes.all(nodes.nind(surfaces(el,:), 2),:);
            line([surfaceNodes2(1,1),surfaceNodes2(2,1)],[surfaceNodes2(1,2),surfaceNodes2(2,2)],[surfaceNodes2(1,3),surfaceNodes2(2,3)],'Color',linecolors(el,:));
            line([surfaceNodes2(2,1),surfaceNodes2(3,1)],[surfaceNodes2(2,2),surfaceNodes2(3,2)],[surfaceNodes2(2,3),surfaceNodes2(3,3)],'Color',linecolors(el,:));
            line([surfaceNodes2(3,1),surfaceNodes2(1,1)],[surfaceNodes2(3,2),surfaceNodes2(1,2)],[surfaceNodes2(3,3),surfaceNodes2(1,3)],'Color',linecolors(el,:));
        end
    end
end