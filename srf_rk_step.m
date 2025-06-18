function [ k_step ] = srf_rk_step( srf, phi, S, c1, c2, nx, u, dt, dx, gam, eps, dirichlet, neumann, pf, D, D_s, r_a1, r_a2, r_d1, r_d2, s_inf)
nx=nx-1;
%S=zeros(nx,1);
k_step = zeros(nx,1);
c1_tilde=bnd(c1)./bnd(phi);
c2_tilde=bnd(c2)./bnd(1-phi);
Delta_s=Compute_Delta(phi,D);
srf_tilde=bnd(srf)./bnd(Delta_s);

r_flux=u*(D.I_r)*srf-D_s*((D.I_r)*Delta_s).*((D.D1_r)*srf_tilde);
l_flux=u*(D.I_l)*srf-D_s*((D.I_l)*Delta_s).*((D.D1_l)*srf_tilde);
J_1s=Delta_s.*(r_a1*c1_tilde.*(s_inf-srf_tilde)-r_d1*srf_tilde);
J_2s=Delta_s.*(r_a2*c2_tilde.*(s_inf-srf_tilde)-r_d2*srf_tilde);

r_flux(end)=r_flux(end)*(1-neumann);
l_flux(1)=l_flux(1)*(1-neumann);
k_step=-dt/dx*(r_flux-l_flux)+dt*J_1s+dt*J_2s;
%k_step=k_step+dt/dx*6*(1-2*phi).*(circshift(S,-1)-S);
