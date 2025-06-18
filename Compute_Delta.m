function [ Delta ] = Compute_Delta(phi, D)
    Delta=abs(6*phi.*(1-phi));
    %Delta=abs((D.D1)*phi);
end