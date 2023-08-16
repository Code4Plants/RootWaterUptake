function kk = piecewk(t,k,z,dx)
    % conductance profiles from piecewise function values
    % dx < 0  return k  else d k/dt 
    for i = 1:length(t)
        if t(i)<=z(1) 
            if (dx < 0)
                kk(i) = (k(1) + (k(2)-k(1))/z(1)*t(i)); 
            else
                kk(i)  = ((k(2)-k(1))/z(1)); 
            end
            continue ;
        end
        if t(i)<=z(2)
            if (dx < 0)
                kk(i)  = (k(2) + (k(3)-k(2))/(z(2)-z(1))*(t(i)-z(1))) ; 
            else
                kk(i)  = (k(3)-k(2))/(z(2)-z(1)); 
            end
            continue;
        end
        if t(i)<=z(3)
            if (dx < 0)
                kk(i) = (k(3) + (k(4)-k(3))/(z(3)-z(2))*(t(i)-z(2))); 
            else
                kk(i)  =  (k(4)-k(3))/(z(3)-z(2)); 
            end
            continue;
        end
        if (dx < 0)
            kk(i) = k(4);
        else
            kk(i) = 0.0; 
        end
    end

