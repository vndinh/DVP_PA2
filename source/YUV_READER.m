function frame = YUV_READER(file,width,height,YUV_type,fr_num,channel)

fr_size = width*height;
YUVfr_size = (YUV_type(1).^2 + YUV_type(2).^2 + YUV_type(3).^2) * fr_size;
fseek(file,(fr_num-1)*YUVfr_size,'bof');
f_temp = fread(file,YUVfr_size,'uchar'); 

Y = reshape(f_temp(1:fr_size),width,height);
frame(:,:,1) = uint8(Y)';

if(channel==3)
    U = reshape(f_temp(1+fr_size:(YUV_type(1).^2+YUV_type(2).^2)*fr_size),YUV_type(2)*width,YUV_type(2)*height);
    V = reshape(f_temp(1+(YUV_type(1).^2+YUV_type(2).^2)*fr_size:(YUV_type(1).^2+YUV_type(2).^2+YUV_type(2).^2)*fr_size),YUV_type(3)*width,YUV_type(3)*height);
    frame(:,:,2) = uint8(imresize(U,2)'); 
    frame(:,:,3) = uint8(imresize(V,2)');
end
frame = double(frame);
end