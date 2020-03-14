function [X] = Fold(X,s,i)
%s是想要生成的张量的规模
d=[i 1:(i-1) (i+1):length(s)];
X=reshape(X,s(d));
X=ipermute(X,d);
end