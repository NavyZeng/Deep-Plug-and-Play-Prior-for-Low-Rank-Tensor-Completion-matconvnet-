function [X] = Unfold(X,s,i)
d=[i 1:(i-1) (i+1):ndims(X)];
X=reshape(permute(X,d),s(i),[]);
end