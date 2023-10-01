function [Edge_mode, eff_no_modes,lambdas] = Compute_EDGE_Modes(a, p, e, t, System_Matrix, global_edges, global_edges_to_domains, no_global_edges, no_modes)
% Compute EDGE basis functions, see Section 3.5

% cell array for storing edge modes
Edge_mode = cell(length(no_global_edges),1);
% array for storing number of effective modes used
eff_no_modes = cell(length(no_global_edges),1);
% array for storing corresponding eigenvalues
lambdas = cell(length(no_global_edges),1);

parfor j = 1:length(no_global_edges) % loop (in parallel) over edges of domain decomposition
    e_loc = e(:, ismember(e(5,:), no_global_edges(j))); % line segments of global mesh that define j-th edge
    % Take only points coordinates and separating boundary and inner points
    e_loc_pts = union(e_loc(1,:), e_loc(2,:)); % Local edge points
    e_loc_bdr = [global_edges(1,j) global_edges(2,j)]; % Global edge boundary
    e_loc_inner = setdiff(e_loc_pts, e_loc_bdr); % Local inner edges
    [S,R] = assemr(p, e_loc, 1, 1, 1, 0); % local stiffness and mass matrix on j-th edge
    % compute domain indices of domains that touch edge
    sub_doms_for_edge = setdiff(global_edges_to_domains(:,j), 0);

    % compute local extension operator E^e, see Section 3.3
    local_basis=sparse(e_loc_inner(:)',1:length(e_loc_inner),1,size(p,2),length(e_loc_inner));
    extension_op = Extension_Eu(p, e, t, local_basis, sub_doms_for_edge, System_Matrix);

    S_loc_dir = S(e_loc_inner,e_loc_inner); % interior stiffness matrix (Dirichlet bc)
    R_loc_dir = R(e_loc_inner,e_loc_inner); % interior mass matrix (Dirichlet bc)
    S_loc_dir=(S_loc_dir'+S_loc_dir)/2; % symmetrize
    R_loc_dir=(R_loc_dir'+R_loc_dir)/2; % symmetrize
    % Solve eigenvalue problem on edges with Dirichlet bc
    %     [V, D] = eig(full(S_loc_dir), full(R_loc_dir));
    no_modes_aux = min(no_modes, length(e_loc_inner)); % if no_modes> local number of degrees of freedom, then limit
    eff_no_modes{j} = no_modes_aux;
    [V,D] = eigs(S_loc_dir, R_loc_dir,no_modes_aux,'sm');
    % post-process eigenvalues and -functions such that they are ordered properly:
    dd=diag(D);
    [dds,inds]=sort(dd,'ascend');
    V=V(:,inds);
    lambdas{j}=dds;
    % select eigenfunctions and extend them to subdomains
    Edge_mode{j} = extension_op * V(:, 1:no_modes_aux);
end