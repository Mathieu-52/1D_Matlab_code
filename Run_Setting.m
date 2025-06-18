classdef Run_Setting 
    properties 
        % one gamma and epsilon value
        single_run = true;
 
        % multiple runs with varying epsilon or gamma values
        eps_run = false;
        gam_run = false;

        % max iterations for intitialization
        max_iter = 1000;

        % single time step
        one_step = false;

        % use log scale
        logit = true;

        % initialization for phi
        initialize = false;

        % reinitialization in order to enforce interface equiblibrium
        reinitialize=false;

        % constant in front of dt
        cfactor = 0.5;

        % lils: small values that are used for 0/0 prevention...
        lil_bnd = eps;
        lil_den = 1E-8;
        lil_cfl = 1E-8;

        % parameter interval for when a parameter sweep is intended
        min_gam_u = 0.02;
        max_gam_u = 2.0;
        min_eps_dx = 0.02;
        max_eps_dx = 8.0;

        % total number of gamma, epsilon, dx values
        ngam=0
        neps=0
        ndx=1    % I have to change it to make spatial convergence 
    end
    methods
        function obj = PF_Setting(opts)
            arguments
                opts.?PF_Setting
            end
           
            for prop = string(fieldnames(opts))'
                obj.(prop) = opts.(prop);
            end
        end
        function set_loops(obj)
            if obj.gam_run
                obj.min_gam_u = 1.0;
            elseif obj.eps_run
                obj.min_eps_dx = 1.0;
            end

            if obj.single_run
                obj.ngam = 0;
                obj.neps = 0;
            elseif obj.gam_run
                obj.neps = 0;
            elseif obj.eps_run
                obj.ngam = 0;
            end

            %obj.ndx = 1; % this is the number of mesh resolutions that are looped over
            if ~obj.single_run
                obj.ndx = 1;   % I have to change it to make spatial convergence 
            end
        end
    end
end


