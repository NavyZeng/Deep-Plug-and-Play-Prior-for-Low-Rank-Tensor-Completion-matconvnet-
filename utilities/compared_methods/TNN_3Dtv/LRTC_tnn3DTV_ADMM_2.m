function [X, Out] = LRTC_tnn3DTV_ADMM(B,Omega,opts)

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'maxit');       maxit = opts.maxit;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'beta');        beta = opts.beta;                end
if isfield(opts, 'lambda');      lambda = opts.lambda;                end
if isfield(opts, 'lambda4');     lambda4 = opts.lambda4;                end
if isfield(opts, 'mu');          mu = opts.mu;              end
if isfield(opts, 'max_beta');    max_beta = opts.max_beta;        end
if isfield(opts, 'max_lambda');  max_lambda = opts.max_lambda;        end
if isfield(opts, 'max_lambda4'); max_lambda4 = opts.max_lambda4;        end

X = B;
Nway = size(B);

Y = zeros(Nway); %% the auxiliary tensor
M = zeros(Nway); %% the Lagrange multiplier

V = zeros(Nway); %% the auxiliary tensor
N = zeros(Nway); %% the Lagrange multiplier

Z1=Y;
Lambda1=Z1;
Z2=Z1;
Lambda2=Z2;
Z3=Z1;
Lambda3=Z3; 

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


D3=zeros(Nway(1),Nway(1),Nway(3));
D3(:,:,1) = -eye(Nway(1));
D3(:,:,end) = eye(Nway(1));

I=zeros(Nway(1),Nway(1),Nway(3));%单位张量
I(:,:,1)=eye(Nway(1));%the first frontal slice 是单位阵

Fa = (1/sqrt(Nway(1)))*fft(eye(Nway(1)));
Fb = (1/sqrt(Nway(2)))*fft(eye(Nway(2)));

Out.Res=[];Out.ResT=[]; Out.PSNR=[];

for k = 1:maxit
    Xold = X;

%% solve Y-subproblem
    Y = prox_tnn(X + M/beta,1/beta);
    
%% solve Z-subproblem
    tau2 = mu./lambda;
    Z1= prox_l1(tprod(D1,X)+Lambda1/lambda(1),tau2(1));
    Z2= prox_l1(tprod(X,D2)+Lambda2/lambda(2),tau2(2));
    Z3= prox_l1(tprod(D3,X)+Lambda3/lambda(3),tau2(3));
%% solve V-subproblem
    V = X + N/lambda4;
    V(Omega) = B(Omega);
        
%% solve X-subproblem
    temp = Y - M/beta;
    tempm = beta*temp;
    
    tempt1 = lambda(1)*tprod(tran(D1),Z1-Lambda1/lambda(1));
    tempt2 = lambda(2)*tprod(Z2-Lambda2/lambda(2),tran(D2));
    tempt3 = lambda(3)*tprod(tran(D3),Z3-Lambda3/lambda(3));
    
    tempt4 = lambda4*(V-N/lambda4);
    
    tempC = tempm+tempt1+tempt2+tempt3+tempt4;
    
    tempA = beta*I+lambda(1)*tprod(tran(D1),D1)+lambda(3)*tprod(tran(D3),D3)+lambda4*I;
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

%% check the convergence
   if isfield(opts, 'Xtrue')
        XT=opts.Xtrue;
        resT=norm(X(:)-XT(:))/norm(XT(:));
        
        Out.ResT = [Out.ResT,resT];
        Out.PSNR = [Out.PSNR,psnr(X,XT)];
    end
    res=norm(X(:)-Xold(:))/norm(Xold(:));
    Out.Res = [Out.Res,res];

%     if k==1 || mod(k, 10) == 0
%         if isfield(opts, 'Xtrue')
%             fprintf('TNN_3DTV: iter = %d   PSNR=%f   res=%f   real-res=%f\n', k, psnr, res, resT);
%         else
%             fprintf('TNN_3DTV: iter = %d   res=%f   \n', k, res);
%         end
%     end
    if res < tol
        break;
    end
%% update Lagrange multiplier
    M = M + beta * (X-Y);

    Lambda1 =  Lambda1 + lambda(1) * (tprod(D1,X)-Z1);
    Lambda2 =  Lambda2 + lambda(2) * (tprod(X,D2)-Z2);
    Lambda3 =  Lambda3 + lambda(3) * (tprod(D3,X)-Z3);

    N = N + lambda4 * (X-V);
     
    beta = min(rho * beta, max_beta);
    lambda = min(rho * lambda, max_lambda);
    lambda4 = min(rho * lambda4, max_lambda4);
end