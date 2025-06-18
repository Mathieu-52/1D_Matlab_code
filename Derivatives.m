function D = Derivatives(n,dx,dirichlet)
%DERIVATIVES Input desired size and mesh size, output is first derivative
%or interpolator function. Boundary condition is periodic for now

%% Derivatives

%1st derivative, central
D.D1 = diag(-1*ones(1,n-1),-1)+diag(ones(1,n-1),1);
D.D1(1,n) = -1;
D.D1(n,1) = 1;
D.D1=D.D1/(2*dx);

%Left sided first derivative
D.D1_l = diag(ones(1,n),0)+diag(-ones(1,n-1),-1);
D.D1_l(1,n) = -1;
D.D1_l=D.D1_l/dx;

%Right sided first derivative
D.D1_r = diag(-ones(1,n),0)+diag(ones(1,n-1),1);
D.D1_r(n,1) = 1;
D.D1_r=D.D1_r/dx; 

% Central difference second derivative
D.D2 = diag(ones(1,n-1),-1)+diag(-2*ones(1,n),0)+diag(ones(1,n-1),1);
D.D2(1,n) = 1;
D.D2(n,1) = 1;
D.D2 = D.D2./(dx^2);

%Left sided interpolation
D.I_l = (1/2)*diag(ones(1,n),0)+(1/2)*diag(ones(1,n-1),-1);
D.I_l(1,n) = 1/2;

%Right sided interpolation
D.I_r = (1/2)*diag(ones(1,n),0)+(1/2)*diag(ones(1,n-1),1);
D.I_r(n,1) = 1/2;

if(dirichlet)
    Dir_Mat=eye(n,n);
    Dir_Mat(1,:)=zeros(1,n);
    Dir_Mat(n,:)=zeros(1,n);

    D.D1=Dir_Mat*D.D1;
    D.D1_l=Dir_Mat*D.D1_l;
    D.D1_r=Dir_Mat*D.D1_r;
    D.D2=Dir_Mat*D.D2;
    D.I_l=Dir_Mat*D.I_l;
    D.I_r=Dir_Mat*D.I_r;    
end

% 
% 
% if isequal(BC,'D')
%     D.D1(n,n-1) = 0;
%     D.D1(n,n) = 0;
%     D.D1_l(n,n-1) = 0;
%     D.D1_l(n,n) = 0;
%     D.D1_r(n,n-1) = 0;
%     D.D1_r(n,n) = 0;
% end
% 
% if isequal(BC,'N')
%     D.D1(n,n-1) = 1;
%     D.D1(n,n) = -1;
%     D.D1_l(n,n-1) = 1;
%     D.D1_l(n,n) = -1;
%     D.D1_r(n,n-1) = 0;
%     D.D1_r(n,n) = 0;
% 
%     D.D1(1,1) = -1;
%     D.D1(1,2) = 1;
%     D.D1_l(1,1) = 0;
%     D.D1_l(1,2) = 0;
%     D.D1_r(1,1) = -1;
%     D.D1_r(1,2) = 1;
% end

end

