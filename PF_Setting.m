classdef PF_Setting 
    properties 
        model="conservativeDI" 
        mobility="nonDegenerate"
        coeffMobility=1.0 
        sigma=1.0
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
    end
end