function [X, Out] = LRTC_tnn2DTV_ADMM(B, Omega, opts)
%%% Solve the 3DTNN-based LRTC problem by ADMM
%
% min_X ||X||_3DTNN, s.t. P_Omega(X) = P_Omega(B)
% where ||X||_3DTNN=sum_{k=1}^{3} alpha_k||X||_TNN_k
% --------------------------------------------
% Input:
%       B       -    the observed tensor with missing elements, 
%                    please make sure B is size of d1*d2*d3 and in range [0, 1].
%       Omega   -    index of the known elements.
%       opts    -    structure value in Matlab. The fields are
%           opts.alpha      -   weights,
%           opts.tol        -   termination tolerance,
%           opts.maxit      -   maximum number of iterations,
%           opts.beta       -   stepsize for dual variable updating in ADMM,
%           opts.max_beta   -   maximum stepsize,
%           opts.rho        -   rho>=1, ratio used to increase beta.
%         
% Output:
%       X       -    the recovery tensor.
%       Out     -    structure value in Matlab. The fields are
%           Out.PSNR        -   PSNR of each step,
%           Out.Res         -   Res of each step, 
%               Res: the relative square error of two successive recovered tensors,
%    stepsize for dual variable updating in ADMM       Out.ResT        -   ResT of each step, 
%               ResT: the relative square error of the recovered tensor and the ground truth tensor.

% Date 20/06/2018
% Written by Yu-Bang Zheng (zhengyubang@163.com)
%%

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'maxit');       maxit = opts.maxit;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'beta');        beta = opts.beta;                end
if isfield(opts, 'lambda');      lambda = opts.lambda;                end
if isfield(opts, 'mu');          mu = opts.mu;              end
if isfield(opts, 'beta');        beta = opts.beta;                end
if isfield(opts, 'max_beta');    max_beta = opts.max_beta;        end
if isfield(opts, 'max_lambda');  max_lambda = opts.max_lambda;        end
if isfield(opts, 'alpha');       alpha = opts.alpha;        end

X = B;
Nway = size(B);

Y = zeros(Nway); %% the auxiliary tensor
M = Y; %% the Lagrange multiplier

Z1=Y;
Lambda1=Z1;
Z2=Z1;
Lambda2=Z2;

D1=zeros(Nway(1),Nway(1),Nway(3));%差分张量D1的大小n1*n1*n3
a=diag(-ones(Nway(1),1));
b=diag(ones(Nway(1)-1,1),1);
tempD=a+b;
tempD(end,1)=1;
D1(:,:,1) = tempD;

D2=zeros(Nway(2),Nway(2),Nway(3));
a=diag(-ones(Nway(2),1));
b=diag(ones(Nway(2)-1,1),-1);
tempD=a+b;
tempD(1,end)=1;
D2(:,:,1) = tempD;

I=zeros(Nway(1),Nway(1),Nway(3));%单位张量
I(:,:,1)=eye(Nway(1));%the first frontal slice 是单位阵

Fa = (1/sqrt(Nway(1)))*fft(eye(Nway(1)));
Fb = (1/sqrt(Nway(2)))*fft(eye(Nway(2)));


Out.Res=[];Out.ResT=[]; Out.PSNR=[];

for k = 1:maxit
    Xold = X;
    
    %% solve Y-subproblem
    tau = 1/beta;
    Y = prox_tnn(X + M/beta,tau);
    
    %% solve Z-subproblem
    tau2 = mu./lambda;
    Z1= prox_l1(tprod(D1,X)+Lambda1/lambda(1),tau2(1));
    Z2= prox_l1(tprod(X,D2)+Lambda2/lambda(2),tau2(2));
   
    
    %% solve X-subproblem
    temp = Y - M/beta;
    
    tempm = beta*temp;
    tempt1 = lambda(1)*tprod(tran(D1),Z1-Lambda1/lambda(1));
    tempt2 = lambda(2)*tprod(Z2-Lambda2/lambda(2),tran(D2));
    
    tempC=tempm+tempt1+tempt2;
    
    tempA = sum(beta)*I+lambda(1)*tprod(tran(D1),D1);
    tempB = lambda(2)*tprod(D2,tran(D2));
    
    tempAf=fft(tempA,[],3);
    tempBf=fft(tempB,[],3);
    tempCf=fft(tempC,[],3);
    
    Xf = zeros(Nway);
    for i=1:Nway(3)
        Ai=tempAf(:,:,i);
        Bi=tempBf(:,:,i);
        Ci=tempCf(:,:,i);
        da=Ai(:,1); deigA=fft(da);
        db=Bi(:,1); deigB=fft(db);        
        Sig=repmat(deigA,1,Nway(2))+repmat(deigB',Nway(1),1);       
        Sig=1./Sig;       
        temp=Sig.*(Fa*Ci*Fb');
        Xf(:,:,i)=Fa'*temp*Fb;
    end 
        
    X = real(ifft(Xf,[],3));
    
    X(Omega) = B(Omega);
    
    %% check the convergence
    if isfield(opts, 'Xtrue')
        XT=opts.Xtrue;
        resT=norm(X(:)-XT(:))/norm(XT(:));
        psnr=myPSNR(X,XT);
        Out.ResT = [Out.ResT,resT];
        Out.PSNR = [Out.PSNR,psnr];
    end
    res=norm(X(:)-Xold(:))/norm(Xold(:));
    Out.Res = [Out.Res,res];
    
    if k==1 || mod(k, 10) == 0
        if isfield(opts, 'Xtrue')
            fprintf('TNN_2DTV: iter = %d   PSNR=%f   res=%f   real-res=%f\n', k, psnr, res, resT);
        else
            fprintf('TNN_2DTV: iter = %d   res=%f   \n', k, res);
        end
    end
    if res < tol 
        break;
    end
    
    %% update Lagrange multiplier
    M = M + beta * (X-Y);
    Lambda1 =  Lambda1 + lambda(1) * (tprod(D1,X)-Z1);
    Lambda2 =  Lambda2 + lambda(2) * (tprod(X,D2)-Z2);
  
    beta = min(rho * beta, max_beta);
    lambda = min(rho * lambda, max_lambda);
    
end
end