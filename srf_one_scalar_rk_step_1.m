function [ k_step ] = srf_one_scalar_rk_step_1(srf_os, phi, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, pf, D, D_s, r_a1, r_a2, r_d1, r_d2, s_inf)
nx=nx-1;

Delta_s=Compute_Delta(phi,D);
K_eff=k_rat*phi+(1-phi);
D_eff=D1*k_rat*phi+D2*(1-phi);


f=r_d1*phi+r_a1*s_inf*Delta_s-srf_os*r_a1;
c1=(-f+sqrt(f.^2+4*r_a1*r_d1*phi.*srf_os))./(2*r_a1);
c1_tilde=(bnd(c1))./(bnd(phi));

r_flux=u*(D.I_r)*srf_os-D1*((D.I_r)*phi).*((D.D1_r)*c1_tilde)-(circshift(S,-1)).*((D.I_r)*c1_tilde)-((D.I_r)*(D_s*Delta_s)).*((D.D1_r)*(bnd(srf_os-c1)./bnd(Delta_s)));
l_flux=u*(D.I_l)*srf_os-D1*((D.I_l)*phi).*((D.D1_l)*c1_tilde)-S.*((D.I_l)*c1_tilde)-((D.I_l)*(D_s*Delta_s)).*((D.D1_l)*(bnd(srf_os-c1)./bnd(Delta_s)));

r_flux(end)=r_flux(end)*(1-neumann);
l_flux(1)=l_flux(1)*(1-neumann);
k_step=-dt/dx*(r_flux-l_flux);