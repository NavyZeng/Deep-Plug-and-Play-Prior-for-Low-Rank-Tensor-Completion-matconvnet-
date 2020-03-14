%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADM algorithm: tensor completion
% paper: Tensor completion for estimating missing values in visual data
% date: 05-22-2011
% min_X: \sum_i \alpha_i \|X_{i(i)}\|_*
% s.t.:  X_\Omega = T_\Omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, his] = myHaLRTC(T, Omega, alpha, beta, maxIter, epsilon, X)

Xtrue = X;
    X = T;%��ʼ����
    X(logical(1-Omega)) = mean(T(Omega));%���ʼ�����еĶ�ʧ������������ֵƽ����
    % X(logical(1-Omega)) = 10;



res = zeros(maxIter, 1);
dim = size(T);
Y = cell(ndims(T), 1);%����������T��ά��ͬ��ģ��cell����, �����滻�ı���
M = Y;%����

% normT = norm(T(:));
for i = 1:ndims(T)
    Y{i} = X;%Y��ÿ��cell��ʼ��ΪX
    M{i} = zeros(dim);%M��ÿ��cell��ʼ��Ϊ��ģͬT��������
end

Msum = zeros(dim);
Ysum = zeros(dim);
for k = 1: maxIter
   
    beta = beta * 1.05;
    
    % update Y
    Msum = 0*Msum;
    Ysum = 0*Ysum;
    for i = 1:ndims(T)
        Y{i} = Fold(Pro2TraceNorm(Unfold(X-M{i}/beta, dim, i), alpha(i)/beta), dim, i);
        Msum = Msum + M{i};
        Ysum = Ysum + Y{i};
    end
    
    % update X
    %X(logical(1-Omega)) = ( Msum(logical(1-Omega)) + beta*Ysum(logical(1-Omega)) ) / (ndims(T)*beta);
    lastX = X;
    X = (Msum + beta*Ysum) / (ndims(T)*beta);
    X(Omega) = T(Omega);
    
    % update M
    for i = 1:ndims(T)
        M{i} = M{i} + beta*(Y{i} - X);
    end
    
    % compute the error
    real(k) = norm(X(:) - Xtrue(:)) / norm(Xtrue(:));
    res(k) = norm(X(:)-lastX(:)) / norm(lastX(:));
%     if mod(k, 20) == 0
%         fprintf('HaLRTC: iterations = %d   difference=%f\n', k, res(k));%ÿ��20�ε�������м���
%     end
   
    if res(k) < epsilon
        break;
    end
end

his.res = res(1:k);
his.real = real;
% fprintf('HaLRTC ends: total iterations = %d   difference=%f\n\n', k, res(k));

