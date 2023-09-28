function [Ephi_vertex, phi] = Compute_VERTEX_Function(a, p, e, t, System_Matrix, global_points, global_edges, global_edges_to_domains,no_global_edges)
% Computation of vertex basis function, see Section 3.4

% vertex functions phi
phi = cell(size(global_points,2),1);
% extension of vertex functions
Ephi_vertex = cell(size(global_points,2),1);
% cell array for storing domain indices that connect to a vertex
subdomains_list_aux = cell(size(global_points,2),1);
% linear function with bc 1/0 on edges pe/ps
f_aux = @(x,ps,pe) sum((x-ps).*(pe-ps))/(norm(pe-ps,2))^2;

for P = 1 : size(global_points,2)
    phi{P} = zeros(size(p,2),1); %Vertex functions
end

% Compute first the edge contributions and then extending them to
% subdomains
for j = 1 : size(global_edges,2) % loop over edges of domain decomposition
    for k = 1 : 2 % loop over both boundary points of edge
        P1 = global_edges(k,j); % index within FEM mesh of vertex to be extended
        P2 = global_edges(setdiff(1:2,k),j); % index within FEM mesh of other point on bdry of edge
        e_loc = e(:, ismember(e(5,:),no_global_edges(j))); % line segments of global mesh that define j-th edge
        e_loc_pts = union(e_loc(1,:), e_loc(2,:)); % select local points on j-th edge
        e_loc_bdr = [global_edges(1,j) global_edges(2,j)]; % global edge boundary
        e_loc_inner = setdiff(e_loc_pts, e_loc_bdr); % local inner edges

        % skeleton function phi: here some function that interpolates
        % continuously from 1 to 0, (piecewise 'linear')
		ind_global_vertex=find(global_points==P1);
        phi{ind_global_vertex}(e_loc_pts) = f_aux(p(1:2,e_loc_pts), p(1:2,P2), p(1:2,P1)).';
        % assemble edge Laplace
        S=assemr(p,e_loc,a,0,0); 
        % Select interior stiffness matrix (Dirichlet bc)
        S_loc_dir = S(e_loc_inner,e_loc_inner); 
        % compute correction to phi s.t. phi-cor is harmonic on edge:
        tmp = S*phi{ind_global_vertex};
        cor = S_loc_dir\tmp(e_loc_inner);
        phi{ind_global_vertex}(e_loc_inner)=phi{ind_global_vertex}(e_loc_inner)-cor;
        % compute domain indices of domains that touch edge
        dom_loc = setdiff(global_edges_to_domains(:,j),0);
        % append domains to subdomain list for vertex
        subdomains_list_aux{ind_global_vertex} = [subdomains_list_aux{ind_global_vertex}; dom_loc(:) ];
    end
end

% compute extensions
for j = 1: size(global_points,2)
    subdomains_list_aux{j} = unique(subdomains_list_aux{j});
    Ephi_vertex{j} = Extension_Eu(p, e, t, phi{j}, subdomains_list_aux{j}, System_Matrix);
end