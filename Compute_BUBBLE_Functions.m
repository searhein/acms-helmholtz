function [basis_bubbles, eff_no_bubbles,lambdas] = Compute_BUBBLE_Functions(p, e, t, S, M, subdomains_list, no_bubbles)
% Computing the bubble basis functions, see Section 3.1

% number of subdomains in domain decomposition
no_ext_dom = length(subdomains_list);
% array for storing bubble functions
basis_bubbles = cell(no_ext_dom,1);
% array specifying how many bubbles are use per subdomain
eff_no_bubbles = cell(no_ext_dom,1);
% array storing the associated eigenvalues
lambdas = cell(length(no_ext_dom),1);

% array storing local eigen(bubble)functions
Vtr=cell(no_ext_dom,1);
% array for interior points of subdomains
p_inner=cell(no_ext_dom,1);

for nsub_ind=1:no_ext_dom % loop (in parallel) over subdomains of domain decomposition
    nsub = subdomains_list(nsub_ind); % index of subdomain

    % construct local mesh:
    % edges of subdomain nsub
    esub = e(:, ismember(e(6,:), nsub) | ismember(e(7,:), nsub));
    % boundary points of subdomain nsub
    psub_bdry = union(esub(1,:), esub(2,:));
    % elements on domain nsub
    tsub = t(4,:) == nsub;
    % points on domain nsub
    psub = unique([t(1,tsub), t(2,tsub), t(3,tsub)]);
    psub_inner = setdiff(psub, psub_bdry); %Internal points
    p_inner{nsub_ind}=psub_inner;

    % project matrix S do local H^1_0-conforming FEM space
    % prolongation matrix
    P_dir = sparse(psub_inner, 1:length(psub_inner), 0*psub_inner+1,size(p,2), length(psub_inner));
    % restriced system matrix
    Ssub = P_dir' * S * P_dir;
    Ssub = (Ssub+Ssub')/2; % symmetrize for eigenvalue solver // works only for self-adjoint S!
    Msub = P_dir'*M*P_dir;
    Msub = (Msub+Msub')/2; % symmetrize for eigenvalue solver
    %     if size(Msub,2)<1000
    %     [Vtr,Dtr] = eig(full(Bsub), full(P_dir'*M*P_dir));
    %     else
    no_bubbles_aux = min(no_bubbles, size(psub_inner,2)); % limit available number of bubble functions
    [Vtr_tmp,Dtr] = eigs(Ssub, Msub,no_bubbles_aux,'sm'); % solve eigenvalue problems
    % post process eigenvalues and -functions to order increasingly
    dd=diag(Dtr);
    [dds,inds]=sort(dd,'ascend');
    lambdas{nsub_ind}=dds;
    eff_no_bubbles{nsub_ind} =no_bubbles_aux;
    Vtr{nsub_ind}=Vtr_tmp(:,inds);
    % if min(norm(dd))<1e-3
    %     fprintf('Attention: In bubble computation there are close to 0 eigenvalues')
    % end
end

% set number of selected bubbles limiting to the current size
for nsub_ind=1:no_ext_dom
    nsub=subdomains_list(nsub_ind);
    basis_bubbles{nsub} = zeros(size(p,2), eff_no_bubbles{nsub_ind});
    basis_bubbles{nsub}(p_inner{nsub_ind},:) = Vtr{nsub_ind};
end
end