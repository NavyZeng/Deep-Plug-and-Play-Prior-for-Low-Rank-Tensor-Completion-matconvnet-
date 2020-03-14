function [X, Out_snn_FFDnet] = admm_snn_DP3LR(Xmiss,Omega,opts)
global sigmas
if isfield(opts,'sigma');       sigma = opts.sigma;                  else  sigma = 0.1;      end
if isfield(opts,'alpha');       alpha = opts.alpha;                  else  alpha = [1,1,1]/3;         end
if isfield(opts,'beta1');       beta1 = opts.beta1;                  else  beta1 = 1;         end
if isfield(opts,'beta2');       beta2 = opts.beta2;                  else  beta2 = 1;         end
if isfield(opts,'maxit');       maxit = opts.maxit;                  else   maxit=200;        end
if isfield(opts,'tol');         tol= opts.tol;                       else  tol = 1e-4;        end
if isfield(opts,'debug');       debug= opts.debug;                   else  debug = 0;         end

Xtrue = opts.Xtrue;
useGPU      = 1;

%% initialization
[h, w, c] = size(Xmiss);     %Œ¨∂»¥Û–°
dim = size(Xmiss);
X = rand(dim);
X(Omega) = Xmiss(Omega);

N = ndims(Xmiss);
Phi = cell(N, 1);
Y = Phi;
for i = 1:N
    Y{i} = X;
    Phi{i} = zeros(dim);
end
Lambda = zeros(dim);

Ysum = zeros(dim);
Psum = zeros(dim);

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

for r = 1:maxit
    Xlast = X;
    %% update Y (SNN)
    Psum = 0*Psum;
    Ysum = 0*Ysum;
    for i = 1:N
        Y{i} = Fold(Pro2TraceNorm(Unfold(X-Phi{i}/beta1, dim, i), alpha(i)/beta1), dim, i);
        Psum = Psum + Phi{i};
        Ysum = Ysum + Y{i};
    end
    
    
    %% update Z £®FFDnet£©
    if c==3
        input = X - Lambda/beta2;
    else
        input = unorigami(X + Lambda/beta2,dim);
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
    X = (Ysum*beta1+Psum+beta2*Z+Lambda)/(N*beta1+beta2);
    X(Omega) = Xmiss(Omega); 
    
    %% stopping criterion
    psnrrec(r) = psnr(X,Xtrue);
    relerr(r) = abs(norm(X(:)-Xlast(:)) / norm(Xlast(:)));
    if  mod(r-1,20) ==0 && debug ==1
        fprintf('%d: RSE:   %f \n',r-1,relerr(r));
        fprintf('    PSNR:  %f \n',psnrrec(r));
    end
    real(r) = abs(norm(X(:)-Xtrue(:)) / norm(Xtrue(:)));
    
    if r > maxit || relerr(r) < tol
        break
    end
    
    %% update Multiplier
    
    for i = 1:N
       Phi{i} = Phi{i} +beta1*(Y{i}-X); 
    end
    Lambda = Lambda + beta2*(Z-X);   
end

Out_snn_FFDnet.res = relerr;
Out_snn_FFDnet.real = real;
Out_snn_FFDnet.psnr = psnrrec;

fprintf('total iterations = %d.',r-1);
end


