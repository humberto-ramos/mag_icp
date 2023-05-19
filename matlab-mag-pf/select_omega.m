function omega = select_omega(k)
    if k<=50
        omega = 0.1;
    elseif k>50 && k<= 80
        omega = -0.2;
    elseif k>80 && k<= 100
        omega = -0.1;
    else
        omega = -0.1;
    end
end