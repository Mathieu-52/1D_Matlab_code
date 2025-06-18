function [output] = pull_back(input)
    output=min(max(input,0),1);
end