function [ k_step, S] = rk_step( phi, nx, u, dt, dx, gam, eps, dirichlet, neumann, pf, D)
nx=nx-1;
k_step = zeros(nx,1);
mu = zeros(nx,1);

% compute mobility
if (pf.mobility == "nonDegenerate")
    if (pf.model=="conservativeDI")
        Mobility=gam*ones(nx,1);
    elseif (pf.model=="CahnHilliard")
        Mobility=pf.coeffMobility*ones(nx,1);
    end
elseif (pf.mobility == "Degenerate_1")
    if (pf.model=="conservativeDI")
        Mobility=gam*(1-(2*pull_back(phi)-1).^2);
    elseif (pf.model=="CahnHilliard")
        Mobility=pf.coeffMobility*(1-(2*pull_back(phi)-1).^2);
    end
elseif (pf.mobility == "Degenerate_2")
    if (pf.model=="conservativeDI")
        Mobility=gam*(1-(2*pull_back(phi)-1).^2).^2;
    elseif (pf.model=="CahnHilliard")
        Mobility=pf.coeffMobility*(1-(2*pull_back(phi)-1).^2).^2;
    end
elseif (pf.mobility == "Degenerate_3")
    if (pf.model=="conservativeDI")
        Mobility=gam*(2*pull_back(phi)-1).^2;
    elseif (pf.model=="CahnHilliard")
        Mobility=pf.coeffMobility*(2*pull_back(phi)-1).^2;
    end
end

% compute right and left fluxes
if (pf.model=="conservativeDI")
    r_flux=-((D.I_r)*Mobility).*(eps*(D.D1_r)*phi-(D.I_r)*(phi.*(1-phi).*(my_sign(circshift(phi,-1)-circshift(phi,1)))));
    l_flux=-((D.I_l)*Mobility).*(eps*(D.D1_l)*phi-(D.I_l)*(phi.*(1-phi).*(my_sign(circshift(phi,-1)-circshift(phi,1)))));
end

if (pf.model=="CahnHilliard")
    free_energy_prime = 4*phi.*(2*phi-1).*(phi-1);
    mu=(pf.sigma/eps)*free_energy_prime-pf.sigma*eps*(D.D2)*phi;
    r_flux=-((D.I_r)*Mobility*eps).*((D.D1_r)*mu);
    l_flux=-((D.I_l)*Mobility*eps).*((D.D1_l)*mu);
end

S=-l_flux;

l_flux=l_flux+u*(D.I_l)*phi;
r_flux=r_flux+u*(D.I_r)*phi;

k_step=-dt/dx*(r_flux-l_flux);
%% Precalculate mu
% theta=10;
% 
% for j = 1:nx
%     jm1 = j-1;
%     jp1 = j+1;
%     jm2 = j-2;
%     jp2 = j+2;
%     if j == 1
%         jm1 = nx;
%         jm2 = nx-1;
%     elseif j == 2
%         jm2 = nx;
%     elseif j == nx-1
%         jp2 = 1;
%     elseif j == nx
%         jp1 = 1;
%         jp2 = 2;
%     end
%     if (pf.model == "CahnHilliard")
%         free_energy_prime = 4* phi(j) * (2*phi(j) - 1) * (phi(j) - 1);
%         mu(j) = (pf.sigma/eps)*free_energy_prime - pf.sigma*(eps/dx^2)*(phi(jm1)-2*phi(j)+phi(jp1));
%     elseif (pf.model == "FloryHuggins")
%         %free_energy_prime = log(pull_back(phi(j)))-log(pull_back(1-phi(j)))-theta*(phi(j)-0.5);
%         mu(j) = - pf.sigma*(eps/dx^2)*(phi(jm1)-2*phi(j)+phi(jp1));
%     end
% end
% 
% 
% 
% %% Compute the RK step 
% for j = 1:nx
%     jm1 = j-1; 
%         jp1 = j+1;
%     jm2 = j-2;
%     jp2 = j+2;
%     if j == 1
%         jm1 = nx;
%         jm2 = nx-1;
%     elseif j == 2
%         jm2 = nx;
%     elseif j == nx-1
%         jp2 = 1;
%     elseif j == nx
%         jp1 = 1;
%         jp2 = 2;
%     end
% 
%     if (pf.model == "conservativeDI")
%         k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1))...
%             +gam*eps*dt/dx^2*(phi(jm1)-2*phi(j)+phi(jp1))...
%             +gam*dt/(2*dx)*(phi(jp1)*(phi(jp1)-1)*my_sign(phi(jp2)-phi(j))...
%             -(phi(jm1)*(phi(jm1)-1)*my_sign(phi(j)-phi(jm2))));
%     elseif (pf.model == "CahnHilliard")
%         if (pf.mobility == "nonDegenerate")
%             k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1)) + pf.coeffMobility*eps*(dt/dx^2)*(mu(jm1)-2*mu(j)+mu(jp1));
%         elseif (pf.mobility == "Degenerate_1")
%             Flux_j = pf.coeffMobility*eps*((2*pull_back(0.5*(phi(j)+phi(jm1)))-1)^2 - 1)^2*(mu(jm1)-mu(j))/dx;
%             Flux_jp1 = pf.coeffMobility*eps*((2*pull_back(0.5*(phi(j)+phi(jp1)))-1)^2 - 1)^2*(mu(j)-mu(jp1))/dx;
%             %Flux_jp1 = pf.coeffMobility*eps*( (((2*phi(j)-1)^2 - 1)^2*mu(j))-(((2*phi(jp1)-1)^2 - 1)^2*mu(jp1)) )/dx;
%             k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1)) + (dt/dx)*(Flux_j-Flux_jp1);
%         elseif (pf.mobility == "Degenerate_2")
%             Flux_j = pf.coeffMobility*eps*(-((2*pull_back(0.5*(phi(j)+phi(jm1)))-1)^2 - 1))*(mu(jm1)-mu(j))/dx;
%             Flux_jp1 = pf.coeffMobility*eps*(-((2*pull_back(0.5*(phi(j)+phi(jp1)))-1)^2 - 1))*(mu(j)-mu(jp1))/dx;
%             %Flux_jp1 = pf.coeffMobility*eps*( (((2*phi(j)-1)^2 - 1)^2*mu(j))-(((2*phi(jp1)-1)^2 - 1)^2*mu(jp1)) )/dx;
%             k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1)) + (dt/dx)*(Flux_j-Flux_jp1);
%         end
%     elseif(pf.model == "FloryHuggins")
%         Flux_j = pf.coeffMobility*eps*(1*(pf.sigma/eps)*(1-theta*0.5*(phi(j)+phi(jm1))*(1-0.5*(phi(j)+phi(jm1))))*(phi(jm1)-phi(j))/dx+...
%             (-((2*pull_back(0.5*(phi(j)+phi(jm1)))-1)^2 - 1))*(mu(jm1)-mu(j))/dx);
%         Flux_jp1 = pf.coeffMobility*eps*(1*(pf.sigma/eps)*(1-theta*0.5*(phi(jp1)+phi(j))*(1-0.5*(phi(jp1)+phi(j))))*(phi(j)-phi(jp1))/dx+...
%             (-((2*pull_back(0.5*(phi(jp1)+phi(j)))-1)^2 - 1))*(mu(j)-mu(jp1))/dx);
%         %Flux_jp1 = pf.coeffMobility*eps*( (((2*phi(j)-1)^2 - 1)^2*mu(j))-(((2*phi(jp1)-1)^2 - 1)^2*mu(jp1)) )/dx;
%         k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1)) + (dt/dx)*(Flux_j-Flux_jp1);
%     end
% end
% 
% 
%     if(dirichlet)    
%         if((j==1)||(j==nx))
%             k_step(j)=0;
%         elseif((j==2)||(j==nx-1))
%             k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1))...
%         +gam*eps*dt/dx^2*(phi(jm1)-2*phi(j)+phi(jp1))...
%         +gam*dt/(2*dx)*(phi(jp1)*(phi(jp1)-1)...
%                       -(phi(jm1)*(phi(jm1)-1)));
%         end
%     end
%     if(neumann)
%               if(j==1)
%                  jm1=j; 
%                k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1))...
%         +gam*eps*dt/dx^2*(phi(jm1)-2*phi(j)+phi(jp1))...
%         +gam*dt/(2*dx)*(phi(jp1)*(phi(jp1)-1)*my_sign(phi(jp2)-phi(j))...
%                       -(phi(jm1)*(phi(jm1)-1)*my_sign(phi(jp1)-phi(j))));
%               elseif(j==nx)
%                   jp1=j;
%                      k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1))...
%         +gam*eps*dt/dx^2*(phi(jm1)-2*phi(j)+phi(jp1))...
%         +gam*dt/(2*dx)*(phi(jp1)*(phi(jp1)-1)*my_sign(phi(j)-phi(jm1))...
%                       -(phi(jm1)*(phi(jm1)-1)*my_sign(phi(j)-phi(jm2))));
%         elseif((j==2)||(j==nx-1))
%             k_step(j) = u*dt/(2*dx)*(phi(jm1)-phi(jp1))...
%         +gam*eps*dt/dx^2*(phi(jm1)-2*phi(j)+phi(jp1))...
%         +gam*dt/(2*dx)*(phi(jp1)*(phi(jp1)-1)*my_sign(phi(jp1)-phi(jm1))...
%                       -(phi(jm1)*(phi(jm1)-1)*my_sign(phi(jp1)-phi(jm1))));
%         end 
%     end
end

