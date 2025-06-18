
% this is a quadratic function verified by c2 analytical 
% Here we calculate the roots of it 

function c2_ana_tilde = c2_analytical( r_a2, r_d2, c_total, s_inf, k_rat,integral_,integral_2)


    a = r_a2 * (integral_ * (k_rat - 1) + 2);
    b = (k_rat - 1) * integral_ * r_d2 + ...
        6 * r_a2 * (integral_ - integral_2) * s_inf + ...
        2*r_d2 - 2 * r_a2 * c_total;
    c = -2 * r_d2 * c_total;

    delta = b^2 - 4*a*c;

    if delta < 0
        c2_ana_tilde = NaN;  % No real root
    elseif (delta == 0)
       disp('The equation has only one real roots');
       disp(-b./(2*a)); 
    else 
        root1 = (-b + sqrt(delta)) / (2*a);
        root2 = (-b - sqrt(delta)) / (2*a);

        % You can also use logic here to choose based on physical relevance
        %c2_ana_tilde = min(root1, root2);  
        c2_ana_tilde = (root1+root2)/2; 

    end
end
