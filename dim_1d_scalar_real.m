% 1D two-phase solver + scalar and surfactant transport modeling
% Based on Mirjalili et al. JCP,2020, Mirjalili et al. IJHMT 2022,
% Mirjalili et al. CTR ARB 2022, and extending Teigen et al. CMS 2009 
% Second order finite differences in space
% Various explicit time stepping schemes
% Various options for mobility and phase field

clear all; close all;format long;


% Set desired plots to true, you can also change font_size, marker_size, and line_width here
my_plot=Plot_Setting(phi=true,srf=true,srf_os=false,srf_total=false,c1=true,c2=true,c_one_scalar=false,line_width=3);


% time stepping scheme
my_time_stepper=Set_Time_Stepping(use_euler=true);
[rk_order,pre_coeff,post_coeff]=my_time_stepper.params();

% Run settings: Only change if you want to do a series of runs while varying epsilon, gamma, mesh size, ...
my_run=Run_Setting();
my_run.set_loops();


% Simulation settings, including physical parameters
nx_base = 129;

% Euler dt (unnecessarily small if not using Euler time-stepping)
use_euler_dt = true;

% constant velocity
u = 0;

% free parameter determining speed of thermal/concentration equilibrium
Amp=1000;

centered = true; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% TEST CASE %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model options include CahnHilliard, conservativeDI. 
% Mobility options are Degenerate_1, Degenerate_2, Degenerate_3 
% and nonDegenerate
% BUT WE ONLY USE NONDEGENERATE MOBILITY 

my_pf=PF_Setting(model="conservativeDI", mobility="nonDegenerate");

%choice (1 = soluble in both and finite s_inf ;' ...
% 2 = realistic solubility in one phase (large diff ratio);' ...
% 3 = confined to one phase and interface;' ...
% 4 = interface confined; ' ...
% 5 = limit of s_inf=0 (not on interface), including fake interface ) 

choice = input("choice = " ); 

switch choice
    case 1 

            k_rat = 1;
            D1 = 1.0; D2 = 1.0; D_s = 0.2;
            c1_t0 = 2.5; c2_t0 = 0.8; srf_0 = 0.0;
            t_final = 3;
        
            r_a1 = 3; r_a2 = 1; r_d1 = 1;
            r_d2 = (r_d1 * r_a2) / (k_rat * r_a1);  
            s_inf = 1.0;



    case 2 
            k_rat = 1;
            D1 = 1.0; D2 = 0.000001; D_s = 0.2;
            c1_t0 = 2.0; c2_t0 = 0.0; srf_0 = 0.3;
            t_final = 3.;
        
            r_a1 = 3; r_a2 = 0.1; r_d1 = 1;
            r_d2 = (r_d1 * r_a2) / (k_rat * r_a1);  
            s_inf = 1.0;

    case 3 
        % choice phase 1 + interface 
            k_rat = 1;
            D1 = 0.5; D2 = 0.0; D_s = 0.5;
            c1_t0 = 2.0; c2_t0 = 0.0; srf_0 = 1.5;
            t_final = 3;
        
            r_a1 = 3; r_a2 = 0; r_d1 = 2;
            r_d2 = 0; 
            s_inf = 1.0;
        
    case 4 
            k_rat = 1;
            D1 = 0.0; D2 = 0.0; D_s = 0.1;
            c1_t0 = 0.0; c2_t0 = 0.0; srf_0 = 0.5;
            t_final = 0.4;
        
            r_a1 = 0; r_a2 = 0; r_d1 = 0;
            r_d2 = 0;
            s_inf = 1.0;
        
    case 5 
            k_rat = 1;
            D1 = 1.0; D2 = 1.0; D_s = 0.0;
            c1_t0 = 2.; c2_t0 = 1.0; srf_0 = 0.0;
            t_final = 2;   % enought to get equilibrium state 
        
            r_a1 = 2; r_a2 = 2; r_d1 = 1;
            r_d2 = (r_d1 * r_a2) / (k_rat * r_a1);  
            s_inf = 0.0;
        
    otherwise
            error('Invalid choice. Please enter 1, 2, 3, 4 or 5.');
end

          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Henry's constant (c1/c2 across the interface at equilibrium)
HC1=1;
HC2=HC1/k_rat;

equilibrium_reached = false; 

c1_all = [];
c2_all = [];
srf_all = [];
time_vec = [];
c_sum = [];
c_all_os = [];
srf_sum = [];



if  ismember(choice, [2 3])    % differencier selon le cas confined = 1 ou 0 
   
    srf_all_os_0 = [];      % confined = 0 -> naive model (assumong that all is at equilibrium)
    error_t_srf_0 = zeros(length(time_vec)); 

    srf_all_os_1 = [];     % confined = 1 -> good model
    error_t_srf_1 = zeros(length(time_vec)); % contains error at each time


else 
    error_t_srf = zeros(length(time_vec)); % contains error at each time
    srf_all_os = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% boundary conditions, default is periodic
dirichlet=0;
neumann=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    % loop over resolutions
    for idx = 1:my_run.ndx
        
        % number of mesh cells
        nx = (nx_base-1)*2^(idx-1)+1;
        my_nx(idx) = nx;
        

        % mesh spacing
        dx = 2/(nx-1);
    
        % edge locations
        x_e = -1:dx:1;
        x_e=x_e';
        
        % cell center locations
        x_c=(x_e(2:end)+x_e(1:end-1))*0.5;
        
        gams = zeros(1,my_run.ngam+1); % all gamma values
        epss = zeros(my_run.neps+1,1); % all epsilon values


        
        % loop over gamma values
        for r = 1:my_run.ngam+1
            gam = u*((my_run.max_gam_u-my_run.min_gam_u)*(r-1)/my_run.ngam+my_run.min_gam_u);
            if my_run.single_run || my_run.eps_run
                gam =1.0; % default value for gamma
            end
            gams(r) = gam;
            
            % loop over epsilon values
            for q = 1:my_run.neps+1
                
                my_eps = dx*((my_run.max_eps_dx-my_run.min_eps_dx)*(q-1)/my_run.neps+my_run.min_eps_dx);
                if my_run.single_run || my_run.gam_run
                    if (my_pf.model == "CahnHilliard"||my_pf.model == "FloryHuggins")
                        my_eps = 2*dx;
                        my_pf.coeffMobility = (3*my_eps^2);%original
                        %coeffMobility = (3*my_eps);
                        %coeffMobility = 3;
                    else
                        my_eps = dx; % default value for epsilon (my_eps is used becaues eps is a constant in matlab)
                    end
                end
                epss(q) = my_eps;
                          
                D=Derivatives(nx-1,dx,dirichlet);


                
                % initial profile
                % phi = ((1/2*(1+tanh((x_c+0.25)/2/my_eps)) - 1/2*(1+tanh((x_c-0.25)/2/my_eps)))).*(abs(x_c)<0.25)+...
                %     ((1/2*(1+tanh((x_c+0.25)/2/(15*my_eps))) - 1/2*(1+tanh((x_c-0.25)/2/(5*my_eps))))).*(abs(x_c)>0.25);
                %phi=1/2*(1+tanh((x_c+.5)/2/my_eps));
                %phi=0.999;
                  %phi = ((1/2*(1+tanh((x_c+2.25)/2/my_eps)) - 1/2*(1+tanh((x_c-2.25)/2/my_eps))));
                phi = ((1/2*(1+tanh((x_c+.25)/2/my_eps)) - 1/2*(1+tanh((x_c-.25)/2/my_eps))));
                  % scalar field concentration in phase 1 per total volume
                c1 = c1_t0.*phi;
                c1_init = c1;
                
                % scalar field concentration in phase 2 per total volume
                c2 = c2_t0.*(1-phi);
                c2_init = c2;
                
                % total concentration
                c_tot =  c1+c2;
                c_tot_init = c_tot;
                c_one_scalar = c_tot; % we start the one scalar model with the c_tot
                c_one_scalar_i = c_tot;
    
                srf_init=srf_0.*Compute_Delta(phi,D);
                srf=srf_init;
                if  ismember(choice, [2 3])  
                    srf_os_0=srf_init+c1_init+c2_init; 
                    srf_os_1=srf_init+c1_init+c2_init; 
                else 
                    srf_os=srf_init+c1_init+c2_init; 
                end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                c1_all(:,end+1) = c1_init;
                c2_all(:,end+1) = c2_init;
                srf_all(:,end+1) = srf_init;
                time_vec(end+1) = 0;
                c_sum(:,end+1) = c1_init + k_rat*c2_init;
                c_all_os(:,end+1) = c_tot;
                srf_sum(:,end+1) = c1_init + c2_init + srf_init;

                if  ismember(choice, [2 3])  
                    srf_all_os_0(:,end+1) = srf_os_0;   
                    srf_all_os_1(:,end+1) = srf_os_1;
                else 
                    srf_all_os(:,end+1) = srf_os;
                end                

                c_total = dx*sum(c1_init + c2_init + srf_init);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

                % dt calculation       
                if ~use_euler_dt
                    if (my_time_stepper.use_rk2)
                        dt = my_run.cfactor*min(min(dx*dx/(2*gam*my_eps),sqrt(3)*dx/(u+gam)), dx*dx/(2*max(D1,D2))); %RK2
                    end
                    if (my_time_stepper.use_rk3)
                        dt = my_run.cfactor*min(min(2.51*dx*dx/(4*gam*my_eps),sqrt(3)*dx/(u+gam)),2.51*dx*dx/(4*max(D1,D2))); %RK3
                    end
                    if (my_time_stepper.use_rk4)
                        dt = my_run.cfactor*min(min(2.79*dx*dx/(4*gam*my_eps),2.83*dx/(u+gam)),2.79*dx*dx/4/max(D1,D2)); %RK4
                    end
                end
                
                if (my_pf.model == "CahnHilliard"||my_pf.model == "FloryHuggins")
                    dt = my_run.cfactor*(dx*dx);
                else 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
                    if choice == 4 
                        dt = 0.0001;  % I have to impose a dt in this case "interface confined"
                    else 
                        dt=my_run.cfactor*min(min(dx*dx/(2*gam*my_eps),1/(gam*my_eps/dx^2+(u+gam)^2/4/gam/my_eps)),dx*dx/2/max(D1,D2)); % euler
                    end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
                end
                if (my_pf.model == "CahnHilliard")
                    dt=dt/5.0;%ad hoc penalty in time
                end
                
                if (my_run.one_step)
                    t_final = dt;
                end
                
                
                my_dt=dt/4
    
    
                time = 0.0;
                tic
                counter=1;
                is_final_time=false;
    
                it = 0; % compteur d'iterations, instant initial deja fait avant 
                
                while (time < t_final)
                    it = it + 1; 
                    % special treatment for the final time step
                    if (time+dt) > t_final
                        dt = t_final-time;
                        is_final_time=true;
                    end
                    
                    time = time + dt;
     
                     phi_int=phi;
                     c1_int=c1;
                     c2_int=c2;
                     c_one_scalar_int=c_one_scalar;
                     c_one_scalar_i_int=c_one_scalar_i;
                     srf_int=srf;
                     if ismember (choice, [2,3])   
                         srf_os_int_0 = srf_os_0;
                         srf_os_int_1 = srf_os_1;
                     else 
                         srf_os_int=srf_os;
                     end                     


                     phi_np1=phi;
                     c1_np1=c1;
                     c2_np1=c2;
                     c_one_scalar_np1=c_one_scalar;
                     c_one_scalar_i_np1=c_one_scalar_i;
                     srf_np1=srf;
                     if ismember (choice, [2,3])   
                         srf_os_np1_0 = srf_os_0;
                         srf_os_np1_1 = srf_os_1;
                     else 
                         srf_os_np1=srf_os;
                     end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if mod(it, 1000) == 0
                            c1_all(:, end+1) = c1_np1;
                            c2_all(:, end+1) = c2_np1;
                            srf_all(:,end+1) = srf_np1;
                            time_vec(end+1) = time;
                            c_sum(:,end+1) = c1_np1 + k_rat*c2_np1;
                            c_all_os(:,end+1) = c_one_scalar_np1;
                            srf_sum(:,end+1) = c1_np1 + c2_np1 + srf_np1;
                            if  ismember(choice, [2 3])  
                                srf_all_os_0(:,end+1) = srf_os_np1_0;   
                                srf_all_os_1(:,end+1) = srf_os_np1_1;
                            else 
                                srf_all_os(:,end+1) = srf_os_np1;
                            end 
                        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for rk_count=1:1:rk_order
                        [d_phi,S]= rk_step( phi_int, nx, u, dt, dx, gam, my_eps, dirichlet, neumann, my_pf, D);
                        d_c1=c1_rk_step( c1_int, c2_int, phi_int, S, srf_int, nx, u, dt, dx, D1, D2, my_eps, Amp, k_rat, dirichlet, neumann, my_pf, D, r_a1, r_a2, r_d1, r_d2, s_inf);
                        d_c2=c2_rk_step( c1_int, c2_int, phi_int, S, srf_int, nx, u, dt, dx, D1, D2, my_eps, Amp, k_rat, dirichlet, neumann, my_pf, D, r_a1, r_a2, r_d1, r_d2, s_inf);
                        d_c_one_scalar=c_one_scalar_rk_step( c_one_scalar_int, phi_int, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, my_pf, D);
                        d_c_one_scalar_i=inc_c_one_scalar_rk_step( c_one_scalar_i_int, phi_int, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, my_pf, D);
                        d_srf=srf_rk_step(srf_int, phi_int, S, c1_int, c2_int, nx, u, dt, dx, gam, my_eps, dirichlet, neumann, my_pf, D, D_s, r_a1, r_a2, r_d1, r_d2, s_inf);
                        if ismember(choice, [2 3])  
                            d_srf_os_0=srf_one_scalar_rk_step_0(srf_os_int_0, phi_int, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, my_pf, D, D_s, r_a1, r_a2, r_d1, r_d2, s_inf);
                            
                            d_srf_os_1=srf_one_scalar_rk_step_1(srf_os_int_1, phi_int, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, my_pf, D, D_s, r_a1, r_a2, r_d1, r_d2, s_inf);
                        else 
                            d_srf_os=srf_one_scalar_rk_step(srf_os_int, phi_int, S, nx, u, dt, dx, D1, D2, my_eps, k_rat, dirichlet, neumann, my_pf, D, D_s, r_a1, r_a2, r_d1, r_d2, s_inf);
                       
                        end 

                        phi_int= phi + pre_coeff(rk_count)*d_phi;
                        c1_int= c1 + pre_coeff(rk_count)*d_c1;
                        c2_int= c2 + pre_coeff(rk_count)*d_c2;
                        c_one_scalar_int= c_one_scalar + pre_coeff(rk_count)*d_c_one_scalar;
                        c_one_scalar_i_int= c_one_scalar_i + pre_coeff(rk_count)*d_c_one_scalar_i;
                        srf_int= srf + pre_coeff(rk_count)*d_srf;
                        if ismember(choice,[2,3])
                            srf_os_int_0 =srf_os_0 + pre_coeff(rk_count)*d_srf_os_0;
                            srf_os_int_1 =srf_os_1 + pre_coeff(rk_count)*d_srf_os_1;
                        else 
                            srf_os_int=srf_os + pre_coeff(rk_count)*d_srf_os;
                        end 
    
                        phi_np1= phi_np1 + post_coeff(rk_count)*d_phi;
                        c1_np1= c1_np1 + post_coeff(rk_count)*d_c1;
                        c2_np1= c2_np1 + post_coeff(rk_count)*d_c2;
                        c_one_scalar_np1= c_one_scalar_np1 + post_coeff(rk_count)*d_c_one_scalar;
                        c_one_scalar_i_np1= c_one_scalar_i_np1 + post_coeff(rk_count)*d_c_one_scalar_i;   
                        srf_np1= srf_np1 + post_coeff(rk_count)*d_srf;
                        if ismember(choice,[2,3])
                            srf_os_np1_0= srf_os_np1_0 + post_coeff(rk_count)*d_srf_os_0;
                            srf_os_np1_1= srf_os_np1_1 + post_coeff(rk_count)*d_srf_os_1;
                        else 
                            srf_os_np1= srf_os_np1 + post_coeff(rk_count)*d_srf_os;
                        end 
               
                   end
                    
                     phi=phi_np1;
                     c1=c1_np1;
                     c2=c2_np1;
                     c_one_scalar=c_one_scalar_np1;
                     c_one_scalar_i=c_one_scalar_i_np1;
                     srf=srf_np1;
                     if ismember(choice,[2,3])
                          srf_os_0 = srf_os_np1_0;
                          srf_os_1 = srf_os_np1_1;
                     else 
                          srf_os = srf_os_np1;
                     end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
                    % tol = 1e-6;
                    % 
                    % mask = (phi > 1e-2) & (phi < 1 - 1e-2);  % avoid divisions by 0
                    % 
                    % 
                    % eq_diff = abs(c1_np1(mask)./phi(mask) - k_rat * c2_np1(mask)./(1 - phi(mask)));
                    % 
                    % max_err = max(eq_diff);
                    % 
                    % if max_err < tol
                    %     disp("Equilibrium satisfied --------------- YESSSSSSSSSSS")
                    %     disp("time  = "+ time)   % to get a precise final time ( = equilibrium state )
                    % else
                    %     disp("Equilibrium not satisfied")
                    %     continue          
                    % 
                    % end



                    % get c2_tilde_numerical 
                    c2_tilde_num = c2_np1./(1-phi);
                    
                    % get c1_tilde_numerical 
                    %c1_tilde_num = bnd(c1_np1)./bnd(phi);
                    epsilon = 1e-6;
                    phi_threshold = 1e-3;  % seuil minimal pour dire qu'on est dans la phase 1 ou l'interface
                    vec_c1_tilde_num = zeros(size(c1_np1));  % initialisation
                    
                    % Calcul uniquement là où phi est assez grand
                    mask = phi > phi_threshold;
                    vec_c1_tilde_num(mask) = c1_np1(mask) ./ (phi(mask) + epsilon);

                    c1_tilde_num = max(vec_c1_tilde_num); 



                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    T=c_one_scalar./(HC1*phi+HC2*(1-phi));
                    T_i=c_one_scalar_i./(HC1*phi+HC2*(1-phi));
                    T_1=(c1+c2)./(HC1*phi+HC2*(1-phi));%(HC1*phi);
                    T_2=c2./(HC2*(1-phi));
        
                    % plotting evolution of data
              
                    % if((mod(counter,my_plot.plotting_frequency)==1)||(is_final_time))
                    %     my_plotter(x_c,phi,c1,c2,c_one_scalar,c_one_scalar_i,srf,srf_os,T,T_i,T_1,T_2,my_plot);
                    %     pause(0.1)
                    %     hold on
                    %     %plot(x_c,1*(c1_t0*phi-c1),'-.k','linewidth',2);
                    %     hold off
                    %     %               saveas(gcf,sprintf('image_vor_%d',counter),'fig')
                    %     %               set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 12]);
                    %     %     print(gcf,'-depsc','-r100',sprintf('image_%d',counter));
                    % 
                    %     xlabel('$x$','interpreter','latex','fontsize',my_plot.font_size);
                    %     set(gca,'Fontsize',my_plot.font_size)
                    %     pause(0.01)
                    %     hold off
                    %     time_percent=100*time/t_final;
                    %     total=sum(c1+c2+srf);
                    %     total_os=sum(srf_os);
                    % end
                    counter=counter+1;




                end
                
                
                 toc
                
            end
        end

    end


% to save the final c1 and c2 
c1_all(:, end+1) = c1_np1;
c2_all(:, end+1) = c2_np1;
srf_all(:,end+1) = srf_np1;
time_vec(end+1) = time;
c_sum(:,end+1) = c1_np1 + k_rat*c2_np1;
c_all_os(:,end+1) = c_one_scalar_np1;
srf_sum(:,end+1) = c1_np1 + c2_np1 + srf_np1;
if  ismember(choice, [2 3])  
    srf_all_os_0(:,end+1) = srf_os_np1_0;   
    srf_all_os_1(:,end+1) = srf_os_np1_1;
else 
    srf_all_os(:,end+1) = srf_os_np1;
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Validation Test  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MASS CONSERVATION (validate for 1 & 3 scalars model)
mass_init  =  dx*sum(c1_init + c2_init + srf_init);
mass_final = dx*sum(c1_np1 + c2_np1 + srf_np1);


            
if abs(mass_init - mass_final) < 1e-6 % machin accuracy  
    disp("Test " + choice + " validated -> total mass conserved globally");
else
    disp("!!!!!!! Test "  + choice + " Invalid -> mass not conserved !!!!!!!!");
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SEMI-ANALYTICAL SOLUTION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get S_tilde_numerical 
%S_tilde_num = srf_np1./Compute_Delta(phi,D);
% disp(['t = 1 correspond à ', num2str(time_vec(2)), ' secondes.']);
% disp(['t = 2 correspond à ', num2str(time_vec(3)), ' secondes.']);
% disp(['t = 20 correspond à ', num2str(time_vec(21)), ' secondes.']);


if (choice == 1)
    S_ana = ((r_a2*s_inf*mean(c2_tilde_num))/(r_d2+r_a2* ...
        mean(c2_tilde_num)))*Compute_Delta(phi,D);
    S_error = sum(abs(srf_np1-S_ana))/sum(S_ana);
    %S_error = norm((srf_np1-S_ana),inf);
    disp("error test S = " + S_error)
    disp("                                    ")
    figure;
    plot(x_c, srf_np1  , 'r-x', 'LineWidth', 2);
    hold on;
    plot(x_c, S_ana, 'b--', 'LineWidth', 2);
    hold off;
    
    xlabel('x');
    ylabel('');
    legend('S num','S ana');
    title('');
    grid off;


    c1_ana = mean(c2_tilde_num)*k_rat.*phi;
    c1_error = sum(abs(c1_np1-c1_ana))/sum(c1_ana);
    disp(" error test c1  = " + c1_error)
    disp("                                    ") 
    figure;
    plot(x_c, c1_np1, 'r-x', 'LineWidth', 2);
    hold on ;
    plot(x_c, c1_ana  , 'b--', 'LineWidth', 2);
    hold off;
    xlabel('x');
    ylabel('');
    legend('c1 num','c1 ana');
    title('');
    grid off;


    c2_ana = c1_tilde_num.*(1-phi)/k_rat;
    c2_error = sum(abs(c2_np1-c2_ana))/sum(c2_ana);
    disp(" error test c2  = " + c2_error)
    disp("                                    ") 
    figure;
    plot(x_c, c2_np1  , 'r-x', 'LineWidth', 2);
    hold on;
    plot(x_c, c2_ana, 'b--', 'LineWidth', 2);
    hold off;
    xlabel('x');
    ylabel('');
    legend('c2 num','c2 ana');
    title('');
    grid off;


    error_t_srf = zeros(length(time_vec)); % contains error at each time 

    figure;
    hold on;
    for k = 1:length(time_vec)
        plot(x_c, srf_sum(:, k), 'r--'); 
        plot(x_c, srf_all_os(:,k), 'b--');
        error_t_srf(k) = sum(abs(srf_sum(:,k)- ...
            srf_all_os(:,k)))/sum(srf_sum(:,k));
    end
    xlabel('x'); ylabel('');
    title('');
    legend('C1 + C2 + S',' Srf one scalar');
    hold off;


    figure;
    semilogy(time_vec, error_t_srf, 'r-o');
    xlabel('Time (s)');
    ylabel('Error');
    title('');
    grid on;

    
    error_x_srf = zeros(nx_base - 1, length(time_vec));
    for k = 1:length(time_vec)
        for i = 1:nx_base - 1
            error_x_srf(i, k) = abs(srf_sum(i, k) - srf_all_os(i, k)) / abs(srf_sum(i, k));
        end
    end
    figure;
    hold on;
    
    % semilogy(x_c, error_x_srf(:, 1),  'g--'); % t = 0
    % semilogy(x_c, error_x_srf(:, 2), 'r--'); 
    % semilogy(x_c, error_x_srf(:, 3), 'k--');
    % semilogy(x_c, error_x_srf(:, 5), 'm--');
    semilogy(x_c, error_x_srf(:, end), 'b--'); 
    
    xlabel('x');
    ylabel('');
    %legend('t = 0 s', 't = 0.06 s', 't = 0.12 s'  , 't = 0.24 s','Final time = 3 s');
    legend('Final time = 3 s ');
    title('');
    grid off;
    hold off;


    % c1,c2,S,phi at each time on the same graph
    figure;
    hold on; 
    for i = 1 : length(time_vec)
        plot(x_c, c1_all(:,i) , 'b--');
        plot(x_c, c2_all(:,i) , 'r--');
        plot(x_c, srf_all(:,i) ,'m--');  
        plot(x_c, phi , 'k-');
    end 
    xlabel('x');
    ylabel('');
    legend('c1', 'c2', 'S','Phi');
    grid off;
    hold off;


elseif ismember(choice, [2 3])


    % S_ana_tilde = ((r_a1*s_inf*c1_tilde_num)/(r_d1+r_a1* ...
    %     c1_tilde_num));
    % S_ana = S_ana_tilde*Compute_Delta(phi,D);
    % S_error = sum(abs(srf_np1-S_ana))/sum(S_ana);
    % disp("error test S = " + S_error)
    % disp("                                    ")
    % figure;
    % plot(x_c, srf_np1  , 'r-x', 'LineWidth', 2);
    % hold on;
    % plot(x_c, S_ana, 'b--', 'LineWidth', 2);
    % hold off;
    % 
    % xlabel('x');
    % ylabel('');
    % legend('S num','S ana');
    % title('');
    % grid off;

   % I draw c1_ana thanks to the relation between C1 and S at equilibrium
    % c1_ana = (r_d1/r_a1)*(S_ana_tilde/(s_inf-S_ana_tilde)).*phi;
    % c1_error = sum(abs(c1_np1-c1_ana))/sum(c1_ana);
    % disp(" error test c1  = " + c1_error)
    % disp("                                    ") 
    % figure;
    % plot(x_c, c1_np1  , 'r-x', 'LineWidth', 2);
    % hold on;
    % plot(x_c, c1_ana, 'b--', 'LineWidth', 2);
    % hold off;
    % 
    % xlabel('x');
    % ylabel('');
    % legend('c1 num','c1 ana');
    % title('');
    % grid off;




    for k = 1:length(time_vec)
        error_t_srf_1(k) = sum(abs(srf_sum(:,k)- ...
            srf_all_os_1(:,k)))/sum(srf_sum(:,k));
        error_t_srf_0(k) = sum(abs(srf_sum(:,k)- ...
            srf_all_os_0(:,k)))/sum(srf_sum(:,k));
    end 

   


% plot all models

    figure;
    hold on;
    %for k = 1:length(time_vec)
        plot(x_c, srf_sum(:, end), 'r--o'); 
        plot(x_c, srf_all_os_0(:,end), 'b-x');
        plot(x_c, srf_all_os_1(:,end), 'g-x');
    %end
    xlabel('x'); ylabel('');
    title('');
    legend('3-scalars model',' 1-confined-scalar model', '1-naive-scalar model' ...
        , 'Location', 'northeast', 'FontSize', 6);
    hold off;


    figure;
    hold on; 
    semilogy(time_vec,error_t_srf_0, 'r-x');
    semilogy(time_vec,error_t_srf_1, 'b-o');
    hold off;
    xlabel('Time (s)');
    ylabel('');
    title('');
    legend('naive','confined');
    grid on;

    error_x_srf = zeros(nx_base - 1, length(time_vec));
    for k = 1:length(time_vec)
        for i = 1:nx_base - 1
            error_x_srf(i, k) = abs(srf_sum(i, k) - srf_all_os_1(i, k)) / abs(srf_sum(i, k));
        end
    end

   % at each time 
    figure;
    hold on;
    semilogy(x_c, error_x_srf(:, 1),  'g--'); % t = 0
    semilogy(x_c, error_x_srf(:, 2), 'r--'); 
    semilogy(x_c, error_x_srf(:, 5), 'm--');
    semilogy(x_c, error_x_srf(:, 21), 'k--');
    semilogy(x_c, error_x_srf(:, end), 'b--'); 
    xlabel('x');
    ylabel('');
    legend('t = 0 s', 't = 0.06 s' , 't = 0.24 s','t = 1.22 s' ,'Final time = 3 s');
    title('');
    grid off;
    hold off;

   % at equilibrium 
    figure;
    hold on;
    semilogy(x_c, error_x_srf(:, end), 'b--'); 
    xlabel('x');
    ylabel('');
    legend('Final time = 3 s ');
    title('');
    grid off;
    hold off;



    % c1,c2,S,phi at each time on the same graph
    figure;
    hold on; 
    for i = 1 : length(time_vec)
        plot(x_c, c1_all(:,i) , 'b--');
        plot(x_c, c2_all(:,i) , 'r--');
        plot(x_c, srf_all(:,i) ,'m--') ;  
        plot(x_c, phi , 'k-');
    end 
    xlabel('x');
    ylabel('');
    legend('c1', 'c2', 'S','Phi');
    grid off;
    hold off;





elseif (choice == 4)
    S_ana = srf_0.*Compute_Delta(phi,D);
    S_error = sum(abs(srf_np1-S_ana))/sum(S_ana);

    disp("error test S = " + S_error)
    disp("                                    ")

    figure;
    plot(x_c, srf_np1  , 'r-x', 'LineWidth', 2);
    hold on;
    plot(x_c, S_ana, 'b--', 'LineWidth', 2);
    hold off;
    
    xlabel('x');
    ylabel('');
    legend('S num','S ana');
    title('');
    grid off;

elseif (choice == 5) 

    c1_ana = mean(c2_tilde_num)*k_rat.*phi;
    c1_error = sum(abs(c1_np1-c1_ana))/sum(c1_ana);
    disp(" error test c1  = " + c1_error)
    disp("                                    ") 
    figure;
    plot(x_c, c1_np1  , 'r-x', 'LineWidth', 2);
    hold on;
    plot(x_c, c1_ana, 'b--', 'LineWidth', 2);
    hold off;
    
    xlabel('x');
    ylabel('');
    legend('c1 num','c1 ana');
    title('');
    grid off;


    c2_ana = c1_tilde_num.*(1-phi)/k_rat;
    c2_error = sum(abs(c2_np1-c2_ana))/sum(c2_ana);
    disp(" error test c2  = " + c2_error)
    disp("                                    ") 
    figure;
    plot(x_c, c2_np1  , 'r-x', 'LineWidth', 2);
    hold on;
    plot(x_c, c2_ana, 'b--', 'LineWidth', 2);
    hold off;
    
    xlabel('x');
    ylabel('');
    legend('c2 num','c2 ana');
    title('');
    grid off;    


    figure;
    hold on;
    for k = 1:length(time_vec)
        plot(x_c, c1_all(:, k), 'b--');  
        plot(x_c, c2_all(:, k), 'r--'); 
        plot(x_c, srf_all(:,k), 'm-');
        plot(x_c, phi, 'k-');
    end
    xlabel('x'); ylabel('');
    title('');
    legend('c1', 'c2' , 'S' , 'Phi');
    hold off;

    figure;
    hold on;
    for k = 1:length(time_vec)
        plot(x_c, c_sum(:, k), 'r-x'); 
        plot(x_c, c_all_os(:,k), 'b-')
    end
    xlabel('x'); ylabel('');
    title('');
    legend('2-scalars model',' 1-scalar model');
    hold off;


% ERROR     
    error_t_c = zeros(length(time_vec)); % contains error at each time 
    error_x_c = zeros(nx_base - 1, length(time_vec));
    for k = 1:length(time_vec)
        error_t_c(k) = sum(abs(c_sum(:,k)- ...
            c_all_os(:,k)))/sum(srf_sum(:,k));
        for i = 1:nx_base - 1
            error_x_c(i, k) = abs(c_sum(i, k) - c_all_os(i, k)) / abs(c_sum(i, k));
        end
    end

    % error time
    figure;
    semilogy(time_vec,error_t_c, 'r-o');
    xlabel('Time (s)');
    ylabel('Error');
    title('');
    grid on;


    %error space
    figure;
    hold on;
    
    % semilogy(x_c, error_x_c(:, 1),  'g--'); % t = 0
    % semilogy(x_c, error_x_c(:, 2), 'r--'); 
    % semilogy(x_c, error_x_c(:, 3), 'k--');
    % semilogy(x_c, error_x_c(:, 11), 'm--');
    semilogy(x_c, error_x_c(:, end), 'b--'); 
    
    xlabel('x');
    ylabel('');
    %legend('t = 0 s', 't = 0.06 s', 't = 0.12 s'  , 't = 0.62 s','Final time = 2 s');
    legend('Final time = 2 s ');
    title('');
    grid off;
    hold off;
    

end 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % construction of the function phi from data
% xq = linspace(min(x_c), max(x_c), 1000);
% f_phi = @(xq) interp1(x_c, phi, xq, 'linear' , 'extrap');  
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Verification : if the curves of f_phi and phi are stacked in the x domain
% 
% %figure;
% %set(gcf, 'Color', 'w');  % Fond blanc de la figure
% %plot(x_c, f_phi(x_c), 'r-x', 'LineWidth', 2);
% %hold on;
% %plot(x_c, phi, 'b-x', 'LineWidth', 2);
% %hold o1
% % ff;
% 
% %xlabel('x');
% %ylabel('phi');
% %legend('Phi function', 'Phi vector');
% %title('Curves of phi function and vector');
% %grid on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%  Integral of f_phi in the x domain [-1;1]   %%%%%%%%%%%%%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% I1 = sum(phi)*dx;
% I2 = sum(Compute_Delta(phi,D))*dx;
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%    Calcul of c2 analytical     %%%%%%%%%%%%%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if (choice ~= 4 )
%     c2_ana_tilde = c2_analytical(r_a2,r_d2,c_total,s_inf,k_rat,I1,I2); 
% 
%     c2_ana = c2_ana_tilde*(1-phi); 
% 
%     figure;
%     plot(x_c, c2_np1, 'r-x', 'LineWidth', 2);
%     hold on;
%     plot(x_c, c2_ana, 'b-o', 'LineWidth', 2);
%     hold off;
% 
%     xlabel('x');
%     ylabel('c');
%     legend('c2','c2 ana');
%     title('Numerical and analytical c2 results at equilibrium');
%     grid on;
% else 
%     c2_ana_tilde = 0;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%    Calcul of S analytical     %%%%%%%%%%%%%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % using delta = 6phi(1-phi)
% if ~ismember(choice, [4 5])
%     S_tilde_ana = ((r_a2*s_inf*c2_ana_tilde)/(r_d2+r_a2*c2_ana_tilde));
%     S_ana = S_tilde_ana*Compute_Delta(phi,D); 
% 
%     figure;
%     plot(x_c, srf_np1  , 'r-x', 'LineWidth', 2);
%     hold on;
%     plot(x_c, S_ana, 'b-o', 'LineWidth', 2);
%     hold off;
% 
% 
%     xlabel('x');
%     ylabel('S');
%     legend('S','S ana');
%     title('Numerical and analytical S results at equilibrium');
%     grid on;
% elseif (choice == 4)
%     S_ana = srf_0.*Compute_Delta(phi,D);
% 
%     figure;
%     plot(x_c, srf_np1  , 'r-x', 'LineWidth', 2);
%     hold on;
%     plot(x_c, S_ana, 'b-o', 'LineWidth', 2);
%     hold off;
% 
% 
%     xlabel('x');
%     ylabel('S');
%     legend('S','S ana');
%     title('Numerical and analytical S results at equilibrium');
%     grid on;
% end 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%    Calcul of c1 analytical     %%%%%%%%%%%%%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if ~ismember(choice, [3 4])
%     c1_ana = c2_ana_tilde*k_rat*phi;
% 
%     figure;
%     plot(x_c, c1_np1  , 'r-x', 'LineWidth', 2);
%     hold on;
%     plot(x_c, c1_ana, 'b-o', 'LineWidth', 2);
%     hold off;
% 
%     xlabel('x');
%     ylabel('c');
%     legend('c1','c1 ana');
%     title('Numerical and analytical c1 results at equilibrium');
%     grid on;
% end 