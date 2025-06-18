function [ k_step ] = c2_rk_step( c1, c2, phi, S, srf, nx, u, dt, dx, D1, D2, my_eps, A, k_rat, dirichlet, neumann, pf, D, r_a1, r_a2, r_d1, r_d2, s_inf)
nx=nx-1;
%S=zeros(nx,1);
k_step = zeros(nx,1);
c2_tilde=bnd(c2)./bnd(1-phi);
Delta_s=Compute_Delta(phi,D);
srf_tilde=bnd(srf)./bnd(Delta_s);

r_flux=u*(D.I_r)*c2+(circshift(S,-1)).*((D.I_r)*c2_tilde)-D2*(1-(D.I_r)*phi).*((D.D1_r)*c2_tilde);
l_flux=u*(D.I_l)*c2+S.*((D.I_l)*c2_tilde)-D2*(1-(D.I_l)*phi).*((D.D1_l)*c2_tilde);
D_m=D1*D2./(k_rat*D1*(1-phi)+D2*phi+eps);
J_12=-A*D_m.*(k_rat*c2.*phi-c1.*(1-phi))+D_m.*((D.D1)*phi).*((D.D1)*(c1+k_rat*c2));
J_s2=Delta_s.*-(r_a2*c2_tilde.*(s_inf-srf_tilde)-r_d2*srf_tilde);
r_flux(end)=r_flux(end)*(1-neumann);
l_flux(1)=l_flux(1)*(1-neumann);
k_step=-dt/dx*(r_flux-l_flux)+dt*J_12+dt*J_s2;

