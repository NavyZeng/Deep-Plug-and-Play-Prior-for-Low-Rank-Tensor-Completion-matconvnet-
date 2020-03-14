function [X, Out_tSVD_FFDnet] = admm_DP3LR(Xmiss,Omega,opts)
global sigmas
minsigma = 0.01;
if isfield(opts,'sigma');       sigma = opts.sigma;                  else  sigma = 0.1;                end
% if isfield(opts,'lambda');       lambda = opts.lambda;                  else  lambda = 0.05;                end
if isfield(opts,'beta1');       beta1 = opts.beta1;                  else  beta1 = 1e-2;                end
if isfield(opts,'beta2');       beta2 = opts.beta2;                  else  beta2 = 1e-2;                end

if isfield(opts,'maxit');      maxit = opts.maxit;      else   maxit=200;     end
if isfield(opts,'tol');       tol= opts.tol;                  else  tol = 1e-4;                end

if isfield(opts,'rho');       rho= opts.rho;                  else  rho = 1;                end
if isfield(opts,'maxbeta');       maxbeta= opts.maxbeta;                  else  maxbeta = 10000;                end
if isfield(opts,'debug');       debug= opts.debug;                  else  debug = 0;                end
Xtrue = opts.Xtrue;
useGPU      = 1;
% maxsigma = 0.01;
%% initialization
[w, h, c] = size(Xmiss);
% eleNum = prod(size(Xmiss));
X = rand(size(Xmiss));
X(Omega) = Xmiss(Omega);
% Y = X;
Lambda1 = zeros(size(X));
Lambda2 = zeros(size(X));

%% FFDnet parameter

if c == 3
    load(fullfile('FFDNet_Clip_color.mat'));
else
    load(fullfile('FFDNet_Clip_gray.mat'));
end

net = vl_simplenn_tidy(net);
if useGPU
    net = vl_simplenn_move(net, 'gpu') ;
end
% Out_tSVD_FFDnet.z = {};
process = {};
for r = 1:maxit
    Xlast = X;
    %% update Y
%     tempY = permute(X+Lambda1/beta1,[3 2 1]);
        [Y] = prox_TNN(X+Lambda1/beta1,1/beta1);
%     Y = ipermute(Y,[3 2 1]);
    
    %% update Z
    if c==3
        input = X + Lambda2/beta2;
        sigX = Xmiss;
    else
        input = unorigami(X + Lambda2/beta2,[w h c]);
        sigX = unorigami(Xmiss,[w h c]);
    end
    
    %     sigma = sqrt(lambda/beta2^sqrt(1/r));
    %     sigma = max(sigma^sqrt(r),0.02)
    %       sigma = 0.1*sqrt(abs(opts.sigma^2-(norm(input(:)-sigX(:))^2)/eleNum));
    
    
    input = single(input); %
    if c==3
        if mod(w,2)==1
            input = cat(1,input, input(end,:,:)) ;
        end
        if mod(h,2)==1
            input = cat(2,input, input(:,end,:)) ;
        end
    else
        if mod(w,2)==1
            input = cat(1,input, input(end,:)) ;
        end
        if mod(h,2)==1
            input = cat(2,input, input(:,end)) ;
        end
    end
    
    if useGPU
        input = gpuArray(input);
    end
    max_in = max(input(:));min_in = min(input(:));
    input = (input-min_in)/(max_in-min_in);

    sigmas = sigma/(max_in-min_in);
    
    res    = vl_simplenn(net,input,[],[],'conserveMemory',true,'mode','test');
    output = res(end).x;
    
    output(output<0)=0;output(output>1)=1;
    output = output*(max_in-min_in)+min_in;
    
    if c==3
        if mod(w,2)==1
            output = output(1:end-1,:,:);
        end
        if mod(h,2)==1
            output = output(:,1:end-1,:);
        end
    else
        if mod(w,2)==1
            output = output(1:end-1,:);
        end
        if mod(h,2)==1
            output = output(:,1:end-1);
        end
    end
    
    if useGPU
        output = gather(output);
    end
    if c==3
        Z = double(output);
    else
        Z = origami(double(output),[w h c]);
    end
    
    %% update X
    X = (beta1*Y+beta2*Z-Lambda1-Lambda2)/(beta1+beta2);
    X(Omega) = Xmiss(Omega);
    X(X>1) = 1;
    X(X<0) = 0;
    
    %% stopping criterion
    psnrrec(r) = myPSNR(X,Xtrue);
    relerr(r) = abs(norm(X(:)-Xlast(:)) / norm(Xlast(:)));
    if  mod(r-1,20) ==0 && debug ==1
        process = [process,X-Xtrue];
        fprintf('%d: RSE:   %f \n',r-1,relerr(r));
        fprintf('    PSNR:  %f \n',psnrrec(r));
    end
    real(r) = abs(norm(X(:)-Xtrue(:)) / norm(Xtrue(:)));
    
    if r > maxit || relerr(r) < tol
        break
    end
    

    
    %% update Lambda
    Lambda1 = Lambda1 + beta1*(X-Y);
    Lambda2 = Lambda2 + beta2*(X-Z);
    
    if rho ~= 1  &&  r>20
        beta1 = min(maxbeta,beta1*rho);
        beta2 = min(maxbeta,beta2*rho);
        sigma = max(minsigma, sigma/rho);
        if debug ==2
            fprintf('beta=%.4f\n',beta)
        end
    end
    
end
process = [process,X];
Out_tSVD_FFDnet.res = relerr;
Out_tSVD_FFDnet.real = real;
Out_tSVD_FFDnet.psnr = psnrrec;
Out_tSVD_FFDnet.process = process;
fprintf('total iterations = %d.',r-1);
end


