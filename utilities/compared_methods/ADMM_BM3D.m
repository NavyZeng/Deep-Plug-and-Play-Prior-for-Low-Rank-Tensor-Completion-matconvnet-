function [Y, Out] = ADMM_BM3D(Xmiss,Omega,opts)
if isfield(opts,'sigma');       sigma = opts.sigma;                  else  sigma1 = 0.05;                end
if isfield(opts,'beta1');       beta1 = opts.beta1;                  else  beta1 = 1e-2;                end
if isfield(opts,'beta2');       beta2 = opts.beta2;                  else  beta2 = 1e-2;                end
if isfield(opts,'maxit');      maxit = opts.maxit;      else   maxit=150;     end
if isfield(opts,'tol');       tol= opts.tol;                  else  tol = 1e-4;                end
if isfield(opts,'debug');       debug= opts.debug;                  else  debug = 0;                end
Xtrue = opts.Xtrue;
%% initialization
[w, h, ~] = size(Xmiss);

Y = rand(size(Xmiss));
Y(Omega) = Xmiss(Omega);

Z = rand(size(Xmiss));
Z(Omega) = Xmiss(Omega);

Lambda1 = zeros(size(Y));
Lambda2 = zeros(size(Y)); %乘子

t0=tic;
for l = 1:maxit
    Xlast = Y;
    %% update X
    X = (beta1*Y-Lambda1+beta2*Z-Lambda2)/(beta1+beta2);
    X(Omega) = Xmiss(Omega);
    %% update Y
    [Y] = prox_tnn(X+Lambda1/beta1,1/beta1);   
    %% update Z   
    input  = X+Lambda2/beta2;   
    max_in = max(input(:));min_in = min(input(:)); 
    input  = (input-min_in)/(max_in-min_in); %归一化
    sigmas = sigma/(max_in-min_in); %参数缩放同倍
    
    [~, output] = CBM3D(1,input,sigmas);
    
    output(output<0)=0;output(output>1)=1; %归一化
    Z = output*(max_in-min_in)+min_in; %复位数据       
    %% stopping criterion
    psnrrec(l) = myPSNR(X,Xtrue);
    relerr(l) = abs(norm(X(:)-Xlast(:)) / norm(Xlast(:)));
    if  mod(l-1,20) ==0 && debug ==1
        fprintf('%d: RSE:   %f \n',l-1,relerr(l));
        fprintf('    PSNR:  %f \n',psnrrec(l));
    end
    real(l) = abs(norm(X(:)-Xlast(:)) / norm(Xlast(:)));
    if l > maxit || relerr(l) < tol
        break
    end    
    %% update Lambda
    Lambda1 = Lambda1 + beta1*(X-Y);
    Lambda2 = Lambda2 + beta2*(X-Z);   
end
    Out.time = toc(t0);
    Out.res = relerr;
    Out.real = real;
    Out.psnr = psnrrec;
    fprintf('total iterations = %d;  time = %4.2fs \n',l-1,Out.time);
end  