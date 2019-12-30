function [v1, v2] = MAP_LINE(fr1, fr2, v1_init, v2_init, num_iters, sigma)
  [h, w] = size(fr1);
  v1 = v1_init;
  v2 = v2_init;
  lambda_d = 20;
  lambda_l = 100;
  tau = 0.2;
  Tl = 3;
  alpha = 0.01;
  eps = 1e-6;
  hx = [-1, 1];
  hy = [-1; 1];
  
  % Calculate spatial gradient of image for 1-clique line field
  dx = imfilter(fr2, hx);
  dy = imfilter(fr2, hy);
  
  % The Metropolis algorithm
  for k = 1:num_iters
    T = tau/log(k+1);
    if T < 0.001
      break;
    end
    
    for j = 2:h-1
      for i = 2:w-1
        % Likelihood engergy
        Ug = (fr2(j,i)-INTP(fr1,i-v1(j,i),j-v2(j,i)))^2/(2*sigma^2);
        
        % Calculate line field by checking 4 directions
        lf = zeros(4,1);
        if (abs(v1(j,i)-v1(j,i+1))>Tl) || (abs(v2(j,i)-v2(j,i+1))>Tl)
          lf(1) = 1;
        end
        if (abs(v1(j,i)-v1(j,i-1))>Tl) || (abs(v2(j,i)-v2(j,i-1))>Tl)
          lf(2) = 1;
        end
        if (abs(v1(j,i)-v1(j+1,i))>Tl) || (abs(v2(j,i)-v2(j+1,i))>Tl)
          lf(3) = 1;
        end
        if (abs(v1(j,i)-v1(j-1,i))>Tl) || (abs(v2(j,i)-v2(j-1,i))>Tl)
          lf(4) = 1;
        end
        
        % Calculate line field potential using 1-clique
        % For horizontal cliques
        if (lf(1)==1) || (lf(2)==1)
          lp_x1 = alpha/(dx(j,i).^2+eps);
        else
          lp_x1 = 0;
        end
        % For vertical cliques
        if (lf(3)==1) || (lf(4)==1)
          lp_x2 = alpha/(dy(j,i).^2+eps);
        else
          lp_x2 = 0;
        end
        lp_c1 = lp_x1 + lp_x2;
        
        % Calculate line field potential using 4-cliques
        if sum(lf)==0
          lp_c4 = 0;
        elseif (sum(lf)==1) || (sum(lf)==4)
          lp_c4 = 2.7;
        elseif (sum(lf)==2) && ((lf(1)+lf(2)==2) || (lf(3)+lf(4)==2))
          lp_c4 = 0.9;
        elseif (sum(lf)==2) || (sum(lf)==3)
          lp_c4 = 1.8;
        end
        
        % Calculate line field energy
        Ul = lp_c1 + lp_c4;
        
        % Motion field using 2-pixel cliques
        Ud = 0;
        Ud = Ud + ((v1(j,i)-v1(j,i+1))^2 + (v2(j,i)-v2(j,i+1))^2)*(1-lf(1));
        Ud = Ud + ((v1(j,i)-v1(j,i-1))^2 + (v2(j,i)-v2(j,i-1))^2)*(1-lf(2));
        Ud = Ud + ((v1(j,i)-v1(j+1,i))^2 + (v2(j,i)-v2(j+1,i))^2)*(1-lf(3));
        Ud = Ud + ((v1(j,i)-v1(j-1,i))^2 + (v2(j,i)-v2(j-1,i))^2)*(1-lf(4));
        
        % Calculate potential function
        p = Ug + lambda_d*Ud + lambda_l*Ul;
        
        while(1)
          % Perturb displacement vector randomly1
          d1 = v1(j,i) + 2 * (rand-0.5)*v1(j,i);
          d2 = v2(j,i) + 2 * (rand-0.5)*v2(j,i);
          
          % Update likelihood energy
          Ug = (fr2(j,i)-INTP(fr1,i-d1,j-d2))^2/(2*sigma^2);
          
          % Check line field
          lf_new = zeros(4,1);
          if (abs(d1-v1(j,i+1))>Tl) || (abs(d2-v2(j,i+1))>Tl)
            lf_new(1) = 1;
          end
          if (abs(d1-v1(j,i-1))>Tl) || (abs(d2-v2(j,i-1))>Tl)
            lf_new(2) = 1;
          end
          if (abs(d1-v1(j+1,i))>Tl) || (abs(d2-v2(j+1,i))>Tl)
            lf_new(3) = 1;
          end
          if (abs(d1-v1(j-1,i))>Tl) || (abs(d2-v2(j-1,i))>Tl)
            lf_new(4) = 1;
          end
          
          % Calculate line field potential using 1-clique
          if (lf_new(1)==1) || (lf_new(2)==1)
            lp_x1 = alpha/(dx(j,i).^2+eps);
          else
            lp_x1 = 0;
          end
          if (lf_new(3)==1) || (lf_new(4)==1)
            lp_x2 = alpha/(dy(j,i).^2+eps);
          else
            lp_x2 = 0;
          end
          lp_c1 = lp_x1 + lp_x2;
          
          % Calculate line field potential using 4-cliques
          if sum(lf_new)==0
            lp_c4 = 0;
          elseif (sum(lf_new)==1) || (sum(lf_new)==4)
            lp_c4 = 2.7;
          elseif (sum(lf_new)==2) && ((lf_new(1)+lf_new(2)==2) || (lf_new(3)+lf_new(4)==2))
            lp_c4 = 0.9;
          elseif (sum(lf_new)==2) || (sum(lf_new)==3)
            lp_c4 = 1.8;
          end
          
          % Update line field energy
          Ul = lp_c1 + lp_c4;
          
          % Update motion field energy
          Ud = 0;
          Ud = Ud + ((d1-v1(j,i+1))^2 + (d2-v2(j,i+1))^2)*(1-lf_new(1));
          Ud = Ud + ((d1-v1(j,i-1))^2 + (d2-v2(j,i-1))^2)*(1-lf_new(2));
          Ud = Ud + ((d1-v1(j+1,i))^2 + (d2-v2(j+1,i))^2)*(1-lf_new(3));
          Ud = Ud + ((d1-v1(j-1,i))^2 + (d2-v2(j-1,i))^2)*(1-lf_new(4));
          
          % Update potential function
          p_new = Ug + lambda_d*Ud + lambda_l*Ul;
          
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

