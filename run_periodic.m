%% Implementation of the ACMS method.
% For details see accompanying paper
% An extension of the approximate component mode synthesis method to the heterogeneous Helmholtz equation
% Elena Giammatteo, Alexander Heinlein, Matthias Schlottbom
% https://arxiv.org/abs/2303.06671

% This script is used to produce the computations for the periodic structure example
% -div(a*grad(u))-kappa^2*u = f in Omega = unit square,
% a d_n u + i k beta u = g on dOmega,
% with planewave solution, see Section 5.4.

fprintf('Helmholtz problem for periodic structure with Robin bc on unit square\n');
% load geometry, see Figure 5.3
load periodic.mat
% refine initial mesh sufficiently
no_p_FEM = 1e5; % number of vertices for the FEM grid to realize ACMS method
while size(p)<no_p_FEM
    [p,e,t]=refinemesh(g,p,e,t);
end
[p,eee,tt]=remove_inner_structure(p,e,t,n); % remove rectangular inclusions, which where included only to resolve diffusion coefficient properly
a=set_diffusion_coeffcient(p,t,n); % define diffusion coefficient
e=eee;
t=tt;
% pdesurf(p,t,a) % visualize diffusion coefficient
% pdemesh(p,e,t) % visualize initial mesh and corresponding domain decomposition

% parameters
beta = 1; % used on the Robin boundary
% wavevector k = omega (k1, k2)
omega= 100; % frequency
k1 = 1.0;
k2 = 0.0;
c = 1; % speed of light
kappa = omega/c; % wavenumber

gR = @(x) exp(-1i*kappa*(k1*x(1,:)+k2*x(2,:))).*exp(-100*(x(1,:).^2+(x(2,:)-0.5).^2)); % Robin boundary data source
f = @(p,t,u,time) 0*p(1,:); % volume data

%% FEM Computations (for use as reference solution)
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

tic
% solve P1-FEM solution, using a direct solver
u_FEM=System_Matrix\RHS;
time_FEM=toc;
% Compute L2- and H1-norms of u_FEM
[u_FEM_L2,u_FEM_H1] = Compute_L2error(1, 1, f, g, p, e, t, u_FEM, @(p)0*p(1,:));

fprintf('FEM dofs=%7d. u_FEM computed in %6.2f sec and |u_FEM|_{L2}=%1.2e |u_FEM|_{H1}=%1.2e\n',size(p,2),time_FEM,u_FEM_L2,u_FEM_H1)
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
G=R*gR(p).';
RHS = F + G; % approximates int_{Omega} fv dx + int_{Gamma} gv dGamma

% Computation of vertex functions
fprintf('Compute Vertex modes in '),tic
[Ephi_vertex, phi] = Compute_VERTEX_Function(1, p, e, t, K - kappa^2 * M, global_points, global_edges, global_edges_to_domains,no_global_edges);
Ephi_vertex=[Ephi_vertex{select_vertex}];
fprintf('%f sec ...\n',toc)

% Computation of edge modes
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

    % error between ACMS and FEM solution
    err_f=u_ACMS-u_FEM;
    err_FA_L2=sqrt(real(err_f'*M*err_f));
    err_FA_H1=sqrt(real(err_f'*(K+M)*err_f));
    fprintf('ACMS dofs=%5d #edge=%4d |u_FEM-u_ACMS|_{L2}/|u_FEM|_{L2}=%1.3e |u_FEM-u_ACMS|_{H1}/|u_FEM|_{H1}=%1.3e in %6.2f sec\n',size(u_system,1),length(ind_edges),err_FA_L2/u_FEM_L2,err_FA_H1/u_FEM_H1,time_ACMS)
end
