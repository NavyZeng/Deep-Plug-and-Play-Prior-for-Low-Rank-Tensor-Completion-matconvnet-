function [X, his] = LRTC_TNN(B, Omega, opts)

tol = 1e-5;
maxit = 200;
rho = 1.1;
beta = 1e-4;
max_beta = 1e10;

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'maxit');    maxit = opts.maxit;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'beta');          beta = opts.beta;                end
if isfield(opts, 'max_beta');      max_beta = opts.max_beta;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
if isfield(opts, 'Xtrue');       Xtrue = opts.Xtrue;         end

szB = size(B);
X = rand(szB);
X(Omega) = B(Omega);
X(logical(1-Omega)) = mean(B(Omega)); %%令初始估计中的丢失像素区域像素值平均化

Y = zeros(szB); %% 辅助变量
M = zeros(szB); %% 拉格朗日乘子

his = [];

for k = 1:maxit
    Xold = X;
    % solve Y-subproblem
    [Y] = prox_TNN(Xold - M/beta,1/beta);
    
    % solve X-subproblem
    X = Y + M / beta;
    
    X(Omega) = B(Omega);
    
    % check the convergence
    if isfield(opts, 'Xtrue')
        real(k) = norm(X(:) - Xtrue(:)) / norm(Xtrue(:));
    end
    res(k) = norm(X(:) - Xold(:)) / norm(Xold(:));
%     if mod(k, 20) == 0
%         %fprintf('TNN: iter = %d   diff=%f\n', k, res(k));
%     end
    
    if res(k) < tol && k>50
        break;
    end
    
    % update Lagrange multiplier
    M = M + beta * (Y - X); 

    beta = min(rho * beta, max_beta);
end


%fprintf('TNN ends: total iterations = %d   difference=%f\n', k, res(k));

end