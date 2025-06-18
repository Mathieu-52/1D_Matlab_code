function [ S ] = compute_S( phi, nx, u, dt, dx, gam, eps, dirichlet, neumann, pf)
nx=nx-1;
mu = zeros(nx,1);
S = zeros(nx,1);
%% Precalculate mu
if (pf.model == "CahnHilliard")
    for j = 1:nx
        jm1 = j-1; 
        jp1 = j+1;
        jm2 = j-2;
        jp2 = j+2;
        if j == 1
            jm1 = nx;
            jm2 = nx-1;
        elseif j == 2
            jm2 = nx;
        elseif j == nx-1
            jp2 = 1;
        elseif j == nx
            jp1 = 1;
            jp2 = 2;
        end
        free_energy_prime = 4* phi(j) * (2*phi(j) - 1) * (phi(j) - 1);
        mu(j) = (pf.sigma/eps)*free_energy_prime - pf.sigma*(eps/dx^2)*(phi(jm1)-2*phi(j)+phi(jp1));
    end
end

%%
for j = 1:nx
    jm1 = j-1; 
        jp1 = j+1;
    jm2 = j-2;
    jp2 = j+2;
    if j == 1
        jm1 = nx;
        jm2 = nx-1;
    elseif j == 2
        jm2 = nx;
    elseif j == nx-1
        jp2 = 1;
    elseif j == nx
        jp1 = 1;
        jp2 = 2;
    end
    %S(j) =gam*eps*(phi(j)-phi(jm1))/dx+gam*0.5*(phi(j)*(phi(j)-1)*my_sign(phi(jp1)-phi(jm1))...
    %                  +(phi(jm1)*(phi(jm1)-1)*my_sign(phi(j)-phi(jm2))));
    if (pf.model == "conservativeDI")
        S(j) =gam*eps*(phi(j)-phi(jm1))/dx+gam*0.5*(phi(j)*(phi(j)-1)*my_sign(phi(jp1)-phi(jm1))...
                      +(phi(jm1)*(phi(jm1)-1)*my_sign(phi(j)-phi(jm2))));
    elseif (pf.model == "CahnHilliard")
        if (pf.mobility == "nonDegenerate")
            S(j) = pf.coeffMobility*eps*(mu(j)-mu(jm1))/dx;
        elseif (pf.mobility == "Degenerate_1")
            S(j) = pf.coeffMobility*eps*((2*pull_back(0.5*(phi(j)+phi(jm1)))-1)^2 - 1)^2*(mu(j)-mu(jm1))/dx;
        elseif (pf.mobility == "Degenerate_2")
            S(j) = pf.coeffMobility*eps*( -((2*pull_back(0.5*(phi(j)+phi(jm1)))-1)^2 - 1))*(mu(j)-mu(jm1))/dx;
           elseif (pf.mobility == "Degenerate_3")
            S(j) = pf.coeffMobility*eps*(0.5*(phi(j)+phi(jm1))-0.5).^2*(mu(j)-mu(jm1))/dx;
        end
    end
    %S(j)=0;% would use to test inconsistent formulations
end

