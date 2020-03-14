function [Xm, Omega] = rgb2mosaic(X)

% To rggb
[h, w, c] = size(X);
Xm = zeros(h,w);
Omega = [];

for tub = 1:c
    for hor = 1:w
        for ver = 1:h   
            if tub==1
                tmp = (hor-1)*h+ver;
                if mod(ver,2)==1 && mod(hor,2)==1
                    Omega = [Omega,tmp];
                    Xm(ver,hor) = X(ver,hor,tub);
                end
            elseif tub==2
                tmp = w*h+(hor-1)*h+ver;
                if mod(hor,2)==1
                   if mod(ver,2)==0
                       Omega = [Omega,tmp];
                       Xm(ver,hor) = X(ver,hor,tub);
                   end
                else
                    if mod(ver,2)==1
                        Omega = [Omega,tmp];
                        Xm(ver,hor) = X(ver,hor,tub);
                    end
                end
            else
                tmp = 2*w*h+(hor-1)*h+ver;
                if mod(hor,2)==0 && mod(ver,2)==0
                    Omega = [Omega,tmp];
                    Xm(ver,hor) = X(ver,hor,tub);
                end
            end
        end
    end
end

end