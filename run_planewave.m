%% Implementation of the ACMS method.
% For details see accompanying paper
% An extension of the approximate component mode synthesis method to the heterogeneous Helmholtz equation
% Elena Giammatteo, Alexander Heinlein, Matthias Schlottbom
% https://arxiv.org/abs/2303.06671

% This script is used to produce the computations for the classical Helmholtz example
% -div(a*grad(u))-kappa^2*u = f in Omega = unit disc,
% a d_n u + i k beta u = g on dOmega = unit circle,
% with planewave solution, see Section 5.1.

fprintf('Helmholtz problem with Robin bc on circle and plane wave solution\n');
current_path=pwd; path(path,[current_path '/functions']); % add function path
load data/circle.mat % load geometry as depicted in Figure 5.1 (right)
% pdemesh(p,e,t) % visualize initial mesh and corresponding domain decomposition

% parameters
a = 1; % diffusion coefficient
beta = 1; % used on the Robin boundary

% wavevector k = omega (k1, k2)
omega= 2 % frequency
k1 = 0.6;
k2 = 0.8;
c = 1; % speed of light
kappa = omega/c; % wavenumber

u_ex = @(x) exp(-1i*kappa*(k1*x(1,:)+k2*x(2,:))); % plane wave solution
Laplace_u_ex = @(x) -kappa^2 * u_ex(x); % Laplacian applied to plane wave
gR = @(x) -1i*kappa*(x(1,:)*k1+x(2,:)*k2).*u_ex(x)-1i*beta*u_ex(x); % consistent Robin boundary data
f = @(p,t,u,time) - Laplace_u_ex(p) - kappa^2* u_ex(p); % consistent volume data / f=0 for plane wave

% hyper parameters
no_p_FEM = 1e4; % number of vertices for the FEM grid to realize ACMS method
doFEM=1; % (=1) do P1-FEM computations for comparison, (=0) skip P1-FEM computations

%% FEM Computations (for comparison)
% Number of refinements for mesh convergence study of P1-FEM
if doFEM
    while size(p,2)<no_p_FEM
        [p,e,t]=refinemesh(g,p,e,t); % refine FEM mesh
        % compute mesh connectivity, eout denotes the edges at the outer
        % boundary, where Robin boundary condition is specified
        eout = Compute_MeshVariables(p, e, t);
        % assemble stiffness and mass matrix
        [K,M,~]=assema(p, t, a, 1, 0);
        % assemble load vector for internal sources, using a linear
        % interpolant
        F=M*(f(p,0,0,0).');
        % assemble boundary sesquilinear form
        [~,R,~]=assemr(p,eout,1,1,1,0); % int_{dOmega} u v ds
        % assemble system matrix
        System_Matrix = K - kappa^2 * M - 1i * beta * R;
        % assemble load vector for boundary sources, using a linear
        % interpolant
        G=R*gR(p).';
        % total load vector
        RHS = F + G; % approximates int_{Omega} fv dx + int_{Gamma} gv dGamma

        %% L2 error computation
        tic
        % solve P1-FEM solution, using a direct solver
        u_FEM=System_Matrix\RHS;
        time_FEM=toc;
        % Compute L2 and H1 error between u_FEM and plane wave u_ex
        [L2err,H1err] = Compute_L2error(1, 1, f, g, p, e, t, u_FEM, u_ex);
        fprintf('FEM dofs=%7d |u-u_ex|_{L2}=%1.3e |u-u_ex|_{H1}=%1.3e in %6.2f sec\n',size(p,2),L2err,H1err,time_FEM)
    end
end
% end FEM computations

%% ACMS method
fprintf("\n Do ACMS method\n")
no_edge_modes = 2.^(1:6); % number of edge modes per edge
% get mesh connectivity, see Section 2.4
[eout, intpts, global_points, global_edges, global_edges_to_domains, no_global_edges, subdomains_list, internal_vertices, internal_edges_ind ] = Compute_MeshVariables(p, e, t);
% use all vertices of domain decomposition for ACMS method, indices for vertex set V
select_vertex = 1:length(global_points);
% select all edges of domain decomposition, indices for edge set E
select_edges_ind = 1 : size(global_edges,2);

% construct underlying FEM matrices and load vectors, as in FEM above
[K,M,~]=assema(p, t, a, 1, 0);
F=M*(f(p,0,0,0).');
[S,R,N]=assemr(p,eout,1,1,1,0);
System_Matrix = K - kappa^2 * M - 1i * beta * R;

[Sg,Rg,Ng]=assemr(p,eout,1,1,1,0);
G=R*gR(p).';
RHS = F + G; % approximates int_{Omega} fv dx + int_{Gamma} gv dGamma

% Compute L2- and H1-errors between the linear interpolant and the exact
% solution
[L2interp,H1interp] = Compute_L2error(1, 1, f, g, p, e, t, u_ex(p).', u_ex);
fprintf('underlying FEM dofs=%7d |I_h u_ex -u_ex|_{L2}=%1.3e |I_h u_ex - u_ex|_{H1}=%1.3e\n',size(p,2),L2interp,H1interp)

% Computation of vertex based functions
fprintf('Compute Vertex modes in '),tic
[Ephi_vertex, phi] = Compute_VERTEX_Function(1, p, e, t, K - kappa^2 * M, global_points, global_edges, global_edges_to_domains,no_global_edges);
Ephi_vertex=[Ephi_vertex{select_vertex}];
fprintf('%f sec ...\n',toc)

% Computation of edge eigenfunctions on all subdomains
fprintf('Compute edge modes in '),tic
[Edge_mode, eff_no_modes,lambda_edge] = Compute_EDGE_Modes(a, p, e, t, K - kappa^2 * M, global_edges, global_edges_to_domains, no_global_edges, no_edge_modes(end));
Edge_mode=[Edge_mode{select_edges_ind}];
fprintf('%f sec ...\n',toc)

% loop over different edge modes
for ed_mo=1:length(no_edge_modes)
    % select edge modes
    ind_edges=[];
    ind=0;
    for jj=1:length(select_edges_ind)
        use=min(no_edge_modes(ed_mo),eff_no_modes{jj});
        if no_edge_modes(ed_mo)>eff_no_modes{jj}
            fprintf('Warning: Not enough edge modes available\n')
            break
        end
        ind_edges=[ind_edges ind+(1:use)];
        ind=eff_no_modes{jj}+ind;
    end

    tic
    % ACMS system assembly
    basis_matrix= [Ephi_vertex  Edge_mode(:,ind_edges)];
    ACMS_system = basis_matrix' * System_Matrix * basis_matrix;
    [i,j,k]=find(ACMS_system);
    ind=(abs(k)>1e-10);
    % sparsify ACMS system matrix
    ACMS_system=sparse(i(ind),j(ind),k(ind),size(ACMS_system,1),size(ACMS_system,2));
    % assemble ACMS load vector
    F_system = basis_matrix' * RHS;
    tic
    % solve ACMS system
    u_system = ACMS_system\F_system; % coordinates in ACMS basis
    time_ACMS=toc;
    u_ACMS = basis_matrix * u_system; % coordinates in standard P1 basis

    % L2- and H1-error of ACMS approximation
    [L2err_ACMS_aux, H1err_ACMS_aux,err,pf,tf] = Compute_L2error(1, 1, f, g, p, e, t, u_ACMS, u_ex);

    % error between ACMS and FEM solution
    err_f=u_ACMS-u_FEM;
    err_FA_L2=sqrt(real(err_f'*M*err_f));
    err_FA_H1=sqrt(real(err_f'*(K+M)*err_f));
    fprintf('ACMS dofs=%5d #edge=%4d  |u-u_ex|_{L2}=%1.3e |u-u_ex|_{H1}=%1.3e |u_FEM-u_ACMS|_{L2}=%1.3e |u_FEM-u_ACMS|_{H1}=%1.3e in %6.2f sec\n',size(u_system,1),length(ind_edges),real(L2err_ACMS_aux),real(H1err_ACMS_aux),err_FA_L2,err_FA_H1,time_ACMS)
end