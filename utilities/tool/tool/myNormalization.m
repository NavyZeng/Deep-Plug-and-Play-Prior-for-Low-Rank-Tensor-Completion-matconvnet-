function [X] = myNormalization(label)

sz=size(label);
N = prod(sz(3:end));
X = zeros(sz);
if max(label(:))>1
    for i = 1:N
        Xtemp = label(:,:,i);
        maxX = max(Xtemp(:));
        minX = min(Xtemp(:));
        X(:,:,i) =  (Xtemp-minX)/(maxX-minX);
    end
else
    X = label;
end
end

