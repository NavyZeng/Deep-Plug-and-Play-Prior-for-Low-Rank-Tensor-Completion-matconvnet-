function hsi = hsv2hsi(hsv)

[h,w,c,t] = size(hsv);
hsi=zeros(h,w,c*t);
for i=1:t
   s = (i-1)*c+1;
   e = i*c;
   hsi(:,:,s:e) = hsv(:,:,:,i); 
end

end