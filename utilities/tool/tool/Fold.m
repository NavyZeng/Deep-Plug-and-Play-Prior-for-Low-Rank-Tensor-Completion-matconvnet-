function [X] = Fold(X,s,i)
%s����Ҫ���ɵ������Ĺ�ģ
d=[i 1:(i-1) (i+1):length(s)];
X=reshape(X,s(d));
X=ipermute(X,d);
end