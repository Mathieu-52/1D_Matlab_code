classdef Set_Time_Stepping
    properties 
        use_euler=false 
        use_rk2=false
        use_rk3=false
        use_rk4=false
    end
    methods
        function obj = Set_Time_Stepping(opts)
            arguments
                opts.?Set_Time_Stepping
            end
           
            for prop = string(fieldnames(opts))'
                obj.(prop) = opts.(prop);
            end
        end
        function [rk_order,pre_coeff,post_coeff]=params(obj)
            if(obj.use_euler)
               rk_order=1;
               pre_coeff=[0];
               post_coeff=[1];
            elseif(obj.use_rk2)
               rk_order=2;
               pre_coeff=[0.5 0.5];
               post_coeff=[0 1];
            elseif(obj.use_rk4)
               rk_order=2;
               pre_coeff=[0.5 0.5 1 1.0/6.0];
               post_coeff=[1.0/6.0 1.0/3.0 1.0/3.0 1.0/6.0];
            end
        end
    end
end