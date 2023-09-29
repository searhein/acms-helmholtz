function plot_solutions(p, t, u_FEM, u_ACMS, u_ex)

%% Figures of solution and approximation
figure
pdesurf(p,t,real(u_FEM))
set(gca,'FontSize',18)
title('FEM approx - Real part')

figure
pdesurf(p,t,imag(u_FEM))
set(gca,'FontSize',18)
title('FEM approx - Imaginary part')

figure
pdesurf(p,t,real(u_ACMS))
set(gca,'FontSize',18)
title('ACMS approx - Real part')

figure
pdesurf(p,t,imag(u_ACMS))
set(gca,'FontSize',18)
title('ACMS approx - Complex part')

figure
pdesurf(p,t,real(u_ex(p(1,:),p(2,:))).')
set(gca,'FontSize',18)
title('Exact Solution - Real part')

figure
pdesurf(p,t,imag(u_ex(p(1,:),p(2,:)))')
set(gca,'FontSize',18)
title('Exact Solution - Complex part')
end
