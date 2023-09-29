function [Eu] = Extension_Eu(p, e, t, u, subdomains_list, System_Matrix)
% compute extensions of interface functions to domains specified in subdomain_list, see Section 3.3

% number of domains to extend to
no_ext_dom = length(subdomains_list);
% array for the extension
Eu = zeros(size(p,2),size(u,2));
tsub = cell(no_ext_dom,1); %

% array of points for each subdomain
psub = cell(no_ext_dom,1);
% array of points for boundary of each subdomain
psub_bdry = cell(no_ext_dom,1);
% array of interior points of each subdomain
psub_inner = cell(no_ext_dom,1); 

for nsub_ind = 1 : no_ext_dom % loop over extension domains
    nsub = subdomains_list(nsub_ind);
   
    % boundary of subdomain nsub
    esub = e(:, ismember(e(6,:), nsub) | ismember(e(7,:), nsub)); 
    % boundary points of subdomain nsub
    psub_bdry{nsub_ind} = union(esub(1,:), esub(2,:)); 
    % elements of subdomain nsub
    tsub{nsub_ind} = t(4,:) == nsub; 
    % points of domain nsub
    psub{nsub_ind} = unique([t(1,tsub{nsub_ind}), t(2,tsub{nsub_ind}), t(3,tsub{nsub_ind})]);
    % interior points of domain nsub
    psub_inner{nsub_ind} = setdiff(psub{nsub_ind}, psub_bdry{nsub_ind});
    
    % Solve local Dirichlet problems
    % assemble prolongation matrix of interior points to global mesh
    P_dir = sparse(psub_inner{nsub_ind}, 1:length(psub_inner{nsub_ind}), 0*psub_inner{nsub_ind}+1,size(p,2), length(psub_inner{nsub_ind})); 
    % assemble local system matrix by Galerkin projection
    Bsub = P_dir' * System_Matrix * P_dir;
    % set boundary conditions of extension, using input function u
    Eu(psub_bdry{nsub_ind},:) = u(psub_bdry{nsub_ind},:);
    % set load vector for extension, see (3.5)
    RHS = -P_dir' * System_Matrix * u;
    % solve local Dirichlet problem, and define extension E^{nsub_ind}
    Eu(psub_inner{nsub_ind},:) = (Bsub\RHS) + u(psub_inner{nsub_ind},:); 
end

end