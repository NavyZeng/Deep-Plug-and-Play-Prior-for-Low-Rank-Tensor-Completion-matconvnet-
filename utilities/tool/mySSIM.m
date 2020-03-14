function [s, svec] = mySSIM(y,x)
sz = size(x);
N = prod(sz(3:end));
svec = zeros(1,N);
for k=1:1:N
    svec(k)=ssim(x(:,:,k)*256,y(:,:,k)*256);
end
s = mean(svec);
end

