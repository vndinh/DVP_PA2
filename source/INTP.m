function value = INTP(img, v1, v2)
  [h, w] = size(img);
  
  % Boundary condition
  v1(v1 < 1) = 1;
  v2(v2 < 1) = 1;
  v1(v1 > w) = w;
  v2(v2 > h) = h;
  
  % Get values of 4 pixels at 4 corners of the image
  v11 = floor(v2);
  v12 = v11 + 1;
  v12(v11 == v2) = v11;
  v21 = floor(v1);
  v22 = v21 + 1;
  v22(v21 == v1) = v21;
  
  % Coefficient for bilinear interpolation
  p1 = v1 - floor(v1);
  p2 = v2 - floor(v2);
  
  % Bilinear interpolation
  value = (img(v11,v21)*(1-p1)+img(v11,v22)*p1)*(1-p2) + (img(v12,v21)*(1-p1)+img(v12,v22)*p1)*p2;
end
