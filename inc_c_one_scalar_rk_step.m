function [k_step_i ] = inc_c_one_scalar_rk_step( c, phi, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, pf, D)
nx=nx-1;
k_step = zeros(nx,1);

K_eff=k_rat*phi+(1-phi);
D_eff=D1*k_rat*phi+D2*(1-phi);
c_tilde=(bnd(c))./(bnd(K_eff));

r_flux=u*(D.I_r)*c-((D.I_r)*D_eff).*((D.D1_r)*c_tilde);
l_flux=u*(D.I_l)*c-((D.I_l)*D_eff).*((D.D1_l)*c_tilde);
r_flux(end)=r_flux(end)*(1-neumann);
l_flux(1)=l_flux(1)*(1-neumann);
k_step_i=-dt/dx*(r_flux-l_flux);
