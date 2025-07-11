function [ k_step ] = srf_one_scalar_rk_step_0(srf_os, phi, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, pf, D, D_s, r_a1, r_a2, r_d1, r_d2, s_inf)
nx=nx-1;

Delta_s=Compute_Delta(phi,D);
K_eff=k_rat*phi+(1-phi);
D_eff=D1*k_rat*phi+D2*(1-phi);

f=r_d1*K_eff+r_a1*k_rat*s_inf*Delta_s-srf_os*r_a1*k_rat;
c_os=(-f+sqrt(f.^2+4*r_a1*k_rat*r_d1*K_eff.*srf_os))./(2*r_a1*k_rat);
c_tilde=(bnd(c_os))./(bnd(K_eff));

r_flux=u*(D.I_r)*srf_os-((D.I_r)*D_eff).*((D.D1_r)*c_tilde)-(k_rat-1)*(circshift(S,-1)).*((D.I_r)*c_tilde)-((D.I_r)*(D_s*Delta_s)).*((D.D1_r)*(bnd(srf_os-c_os)./bnd(Delta_s)));
l_flux=u*(D.I_l)*srf_os-((D.I_l)*D_eff).*((D.D1_l)*c_tilde)-(k_rat-1)*S.*((D.I_l)*c_tilde)-((D.I_l)*(D_s*Delta_s)).*((D.D1_l)*(bnd(srf_os-c_os)./bnd(Delta_s)));

r_flux(end)=r_flux(end)*(1-neumann);
l_flux(1)=l_flux(1)*(1-neumann);
k_step=-dt/dx*(r_flux-l_flux);