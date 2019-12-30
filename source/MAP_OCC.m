function [v1, v2] = MAP_OCC(fr1, fr2, v1_init, v2_init, num_iters, sigma)
  [h, w] = size(fr1);
  v1 = v1_init;
  v2 = v2_init;
  To = 500;
  lambda_d = 20;
  lambda_o = 10;
  tau = 0.2;
  
  for k = 1:num_iters
    T = tau/log(k+1);
    if T < 0.001
      break;
    end
    
    for j = 2:h-1
      for i = 2:w-1
        % Likelihood energy
        Ug = (fr2(j,i)-INTP(fr1,i-v1(j,i),j-v2(j,i)))^2/(2*sigma^2);
        
        % Occlusion field using singleton clique
        if Ug > To
          occ = 1;
        else
          occ = 0;
        end
        Uo = occ * To;
        
        % Motion field using 2-pixel clique
        Ud = 0;
        Ud = Ud + (v1(j,i)-v1(j,i+1))^2 + (v2(j,i)-v2(j,i+1))^2 + (v1(j,i)-v1(j,i-1))^2 + (v2(j,i)-v2(j,i-1))^2;
        Ud = Ud + (v1(j,i)-v1(j+1,i))^2 + (v2(j,i)-v2(j+1,i))^2 + (v1(j,i)-v1(j-1,i))^2 + (v2(j,i)-v2(j-1,i))^2;
        
        % Potential function
        p = Ug*(1-occ) + lambda_d*Ud + lambda_o*Uo;
        
        while(1)
          % Perturb displacement vector
          d1 = v1(j,i) + 2*(rand-0.5)*v1(j,i);
          d2 = v2(j,i) + 2*(rand-0.5)*v2(j,i);
          
          % Update the energy
          Ug = (fr2(j,i)-INTP(fr1,i-d1,j-d2))^2/(2*sigma^2);
          Ud = 0;
          Ud = Ud + (d1-v1(j,i+1))^2 + (d2-v2(j,i+1))^2 + (d1-v1(j,i-1))^2 + (d2-v2(j,i-1))^2;
          Ud = Ud + (d1-v1(j+1,i))^2 + (d2-v2(j+1,i))^2 + (d1-v1(j-1,i))^2 + (d2-v2(j-1,i))^2;
          
          % Occlusion field using 1 clique
          if Ug > To
            occ = 1;
          else
            occ = 0;
          end
          Uo = occ * To;
          
          % Update the potential function
          p_new = Ug*(1-occ) + lambda_d*Ud + lambda_o*Uo;
          %p_new = Ug + lambda_d*Ud + lambda_o*Uo;
          
          % Accept or not
          dp = p_new - p;
          if (dp<=0) || ((exp(-dp/T)>rand) && (exp(-dp/T)<1))
            v1(j,i) = d1;
            v2(j,i) = d2;
            break;
          end
        end
      end
    end
  end
end

