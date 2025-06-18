classdef Plot_Setting 
    properties 
        phi=false 
        c1=false
        c2=false
        c1pc2=false 
        c_one_scalar=false 
        c_one_scalar_i=false 
        srf=false
        srf_os=false
        srf_total=false
        T=false 
        T_i=false
        T_1=false 
        T_2=false
        font_size = 15;
        marker_size = 10;
        line_width = 2.0;
        plotting_frequency=200;
    end
    methods
        function obj = Plot_Setting(opts)
            arguments
                opts.?Plot_Setting
            end
           
            for prop = string(fieldnames(opts))'
                obj.(prop) = opts.(prop);
            end
        end
    end
end