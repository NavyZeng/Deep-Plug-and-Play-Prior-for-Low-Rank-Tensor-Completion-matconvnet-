function hsv = hsi2hsv(hsi,c,t)

[h,w,tt] = size(hsi);
hsv=zeros(h,w,c,t);
for i=1:c:tt
   hsv(:,:,:,fix(i/c)+1)=hsi(:,:,i:(i+c-1)); 
end

end