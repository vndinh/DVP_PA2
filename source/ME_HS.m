function [v1, v2] = ME_HS(fr1, fr2, block_size)
  %% Motion Estimation using Horn-Schunck algorithm
  %% Arguments
  % fr1         : the frame at time t
  % fr2         : the frame at time t+1
  % block_size  : the size of pixel block
  % alpha       : constant in the Horn-Schunck iterative equation
  % v1, v2      : the motion vectors
  
  %% Implementation
  [height, width] = size(fr1);
  alpha = 10;
  
  nb_x1 = floor(width/block_size);  % Number of blocks in the horizontal dimension
  nb_x2 = floor(height/block_size); % Number of blocks in the vertical dimension
  
  % Motion vectors
  v1 = zeros(nb_x2, nb_x1);
  v2 = zeros(nb_x2, nb_x1);
  
  % Partial derivatives w.r.t x1, x2 and t
  Kx1 = 0.25*[-1 1; -1 1];
  dx1 = conv2(fr1, Kx1, 'same') + conv2(fr2, Kx1, 'same');
  Kx2 = 0.25 * [-1 -1; 1 1];
  dx2 = conv2(fr1, Kx2, 'same') + conv2(fr2, Kx2, 'same');
  Kt = 0.25 * ones(2);
  dt = conv2(fr1, Kt, 'same') + conv2(fr2, -Kt, 'same');
  
  % Iteration
  % Laplacian filter
  L = [1/12 1/6 1/12; 1/6 0 1/6; 1/12 1/6 1/12];
  num_iter = 50;
  for i = 1:num_iter
    v1_bar = conv2(v1, L, 'same');
    v2_bar = conv2(v2, L, 'same');
    for x2 = 1:block_size:(height-block_size+1)
      for x1 = 1:block_size:(width-block_size+1)
        x2_blk = floor((x2-1)/block_size) + 1;
        x1_blk = floor((x1-1)/block_size) + 1;
        x2_dest = x2 + block_size - 1;
        x1_dest = x1 + block_size - 1;
        dx1_blk = dx1(x2:x2_dest, x1:x1_dest);
        dx2_blk = dx2(x2:x2_dest, x1:x1_dest);
        dt_blk = dt(x2:x2_dest, x1:x1_dest);
        
        D = alpha^2 + dx1_blk.^2 + dx2_blk.^2;
        N1 = dx1_blk .* v1_bar(x2_blk, x1_blk) + dx2_blk .* v2_bar(x2_blk, x1_blk) + dt_blk;
        
        v1(x2_blk, x1_blk) = v1_bar(x2_blk, x1_blk) - mean(mean(dx1_blk.*N1./D));
        v2(x2_blk, x1_blk) = v2_bar(x2_blk, x1_blk) - mean(mean(dx2_blk.*N1./D));
      end
    end
  end
end

