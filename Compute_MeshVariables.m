%% Mesh variables: global/local/boundary edges/vertices
function [eout, intpts, global_points, global_edges, global_edges_to_domains, no_global_edges, subdomains_list, internal_vertices, internal_edges_ind ] = Compute_MeshVariables(p, e, t)

    eout=e(:,e(6,:)==0|e(7,:)==0); %External boundary siding with 0
    bdrpts=union(eout(1,:),eout(2,:)); %Boundary points
    intpts=setdiff(1:size(p,2),bdrpts); %Returns internal points

    % no_global_edges = unique edge numbers in edge label
    % ie = indices of different labels (one per label)
    [no_global_edges,ie] = unique(e(5,:));
    % Store left and right domains of labeled edges
    global_edges_to_domains = e(6:7, ie);
    %List of subdomains excluding outer domain 0
    subdomains_list = setdiff(unique(global_edges_to_domains), 0);
    % Setting variables to store global edges
    global_edges = zeros(2, length(no_global_edges));

    % Determining global edges, boundaries and indices
    for j = 1:length(no_global_edges) % Looping over # global edges
        %Storing local edges which are on the selected global edge j
        e_loc = e(1:2, ismember(e(5,:), no_global_edges(j)));
        % Selecting edges which appear only once = boundary edges
        start_end_ind = setxor(e_loc(1,:), e_loc(2,:));
        % The global edge is given by boundary indices
        global_edges(:,j) = start_end_ind(:);
    end
    global_points = unique(global_edges)';
    %Find all global edges that are not on the boundary
    internal_vertices = intersect(global_points, intpts);
    internal_edges_ind = find(sum(global_edges_to_domains(:,1:end)~=0)-1);
end