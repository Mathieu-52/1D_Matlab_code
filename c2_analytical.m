
% this is a quadratic function verified by c2 analytical 
% Here we calculate the roots of it 

function c2_ana_tilde = c2_analytical( r_a2, r_d2, c_total, s_inf, k_rat,I1,I2)


    a = r_a2 * (I1 * (k_rat - 1) + 2);
    b = (r_d2*(I1*(k_rat-1)+2) + r_a2*(s_inf*I2-2*c_total));
    c = -2* r_d2 * c_total;

    delta = b^2 - 4*a*c;

    if delta < 0
        c2_ana_tilde = NaN;  % No real root
       
    elseif (delta == 0)
       disp('The equation has only one real roots');
       disp(-b./(2*a)); 
    else 
        root1 = (-b + sqrt(delta)) / (2*a);
        root2 = (-b - sqrt(delta)) / (2*a);
    

        c2_ana_tilde = max(root1, root2);  

    end
end
