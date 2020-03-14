function [X, Out_tSVD_TNN] = tSVD_TNN(Xmiss,Omega,opts)


if isfield(opts,'beta');       beta = opts.beta;                  else  beta = 1e-2;                end
if isfield(opts,'maxit');      maxit = opts.maxit;      else   maxit=500;     end
if isfield(opts,'tol');       tol= opts.tol;                  else  tol = 1e-4;                end
Xtrue = opts.Xtrue;


%% initialization

X = Xmiss;
Y = Xmiss;
Lambda = zeros(size(Y));

t0=tic;
for r = 1:maxit
    Xlast = X;
    %% update X
    [Y] = prox_tnn_dft(Xlast - Lambda/beta,1/beta);
    
    
    %% update Y
    X = Y + Lambda/beta;
    X(Omega) = Xmiss(Omega);
    
    psnrrec(r) = myPSNR(X,Xtrue);
    relerr(r) = abs(norm(X(:)-Xlast(:)) / norm(Xlast(:)));
    %% stopping criterion
     if  mod(r-1,20) ==0
        fprintf('%d: RSE:   %f \n',r-1,relerr(r));
        fprintf('    PSNR:  %f \n',psnrrec(r));
    end
     
    real(r) = abs(norm(X(:)-Xtrue(:)) / norm(Xtrue(:)));
    if r > maxit || relerr(r) < tol
        break
    end
    
    %% update Lambda
    Lambda = Lambda + beta*(Y-X);

    %      disp(['step ' num2str(r) '  ||M-A_k||_F^2 = ' num2str(error)]);

end

Out_tSVD_TNN.time = toc(t0);
Out_tSVD_TNN.res = relerr;
Out_tSVD_TNN.real = real;
Out_tSVD_TNN.psnr = psnrrec;

fprintf('total iterations = %d;  time = %4.2fs \n',r-1,Out_tSVD_TNN.time);

end