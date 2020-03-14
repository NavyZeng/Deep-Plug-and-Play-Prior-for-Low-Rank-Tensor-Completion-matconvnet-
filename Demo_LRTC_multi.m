clc
clear; close all;
addpath(genpath(cd));
rng(213412);

methodname  = {'Observed','FaLRTC', 'TNN','TNN_3DTV','DP3LRTC'};
EN_FaLRTC   = 1;
EN_TNN     = 1;
EN_TNN_3DTV  = 1;
EN_DP3LRTC    = 1;
Mnum = length(methodname);

folderTest  = 'testsets';
folderResult= 'results';
imageSets   = {'img','msi','video'};         % testing datasets
setTestCur  = imageSets{2};      % current testing dataset

printFig  = 1;
figSlice = [1 2 31];
figTimes = 2;
pauseTime   = 0;
SR_case = 0;
for SR = [0.05 0.1 0.2]  % 0.2 0.3 0.4
    fno = 32;
    fprintf('SR=%.2f, Fno=%d\n',SR,fno);
    SR_case = SR_case+1;
    folderResultCur       =  fullfile(folderResult, [setTestCur,'_',num2str(SR)]);
    if ~isfolder(folderResultCur)
        mkdir(folderResultCur)
    end
    
    % read images
    ext         =  {'*.mat'};
    filePaths   =  [];
    
    for jj = 1 : length(ext) %length(ext)
        filePaths = cat(1,filePaths, dir(fullfile(folderTest,setTestCur,ext{jj})));
    end
    Re_tensor	=  cell(Mnum,length(filePaths));
    MPSNRALL	=  zeros(Mnum,length(filePaths));
    SSIMALL     =  zeros(Mnum,length(filePaths));
    FSIMALL     =  zeros(Mnum,length(filePaths));
    Time_all    =  zeros(Mnum,length(filePaths));
    Image_name_all  = cell(length(filePaths),1);
    for jj = [1 2 3 4 6 7 8] %length(filePaths) 7 8 10
        % read images
        disp([filePaths(jj).name])
        load([filePaths(jj).name])
        sz=size(X);
        fno = min(fno,sz(3));
        label = X(:,:,1:fno);
        sz=size(label);
        w = sz(2)*figTimes;
        h = sz(1)*figTimes;
        [~,Image_Name,extCur] = fileparts(filePaths(jj).name);
        disp(['Image #' num2str(jj) ': ' Image_Name ])
        Image_name_all{jj} = Image_Name;
        X_GT = myNormalization(label);
        
        % Generate known data
        P = round(SR*prod(sz));      % prod·µ»Ø³Ë»ý
        Omega = randsample(prod(sz),P);
        [Omega,~] = sort(Omega);   % Known = Omega
        % Missing data
        Xkn          = X_GT(Omega);
        X_Observed        = zeros(sz);
        X_Observed(Omega) = Xkn;
        
        Omega2     = zeros(sz);
        Omega2(Omega) = 1;
        Omega2     = (Omega2 > 0);
        
 %% Observed data
        ii = 1;
        disp([methodname{ii}, ':']);
        Re_tensor{ii,jj} = X_Observed;
        % quantitative metrics
%         MPSNRALL(ii,jj) = psnr(Re_tensor{ii,jj},X_GT);
%         SSIMALL(ii,jj)  = ssim(uint8(Re_tensor{ii,jj}*255),uint8(X_GT*255));
        PPPSNRV = [];SSSIMV = [];
            for kkk = 1:size(X_Observed,3)
                PPPSNRV(kkk) = psnr(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk),1);
            end
            MPSNRALL(ii,jj) = mean(PPPSNRV);%psnr(Re_tensor{ii,jj},X_GT);
            for kkk = 1:size(X_Observed,3)
                SSSIMV(kkk) = ssim(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk));
            end
            SSIMALL(ii,jj) = mean(SSSIMV);
        Time_all(ii,jj) = 0;
        fprintf(' %8.8s    %5.4s    %5.4s    %5.4s \n','method','PSNR', 'SSIM', 'time');
        fprintf(' %8.8s    %5.3f    %5.3f    %.1f \n',methodname{ii},MPSNRALL(ii,jj) , SSIMALL(ii,jj) , Time_all(ii,jj));
        
%                 if printFig ==1
%                     filename =[ 'F:\matlab\tSVD FFDnet\figs\',setTestCur,'\',newname,'-',num2str(100*SR), '.jpg'];
%                     figure,
%                     imshow(Xmiss(:,:,figSlice),'border','tight','initialmagnification','fit');
%                     set (gcf,'Position',[0,0,w,h]);
%                     axis normal;
%                     print(gcf,'-djpeg',filename);
%         
%                     filename2 =[ 'F:\matlab\tSVD FFDnet\figs\',setTestCur,'\',newname, '.jpg'];
%                     figure,
%                     imshow(Xtrue(:,:,figSlice),'border','tight','initialmagnification','fit');
%                     set (gcf,'Position',[0,0,w,h]);
%                     axis normal;
%                     print(gcf,'-djpeg',filename2);
%                 end
        
        %% use SNN
        ii = 2;
        if EN_FaLRTC
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 200;
            alpha = [1, 1, 1];
            alpha = alpha/sum(alpha);
            rho = 1e-2;
            
            fprintf('\n');
            disp(['performing ',methodname{ii}, ' ... ']);
            
             tStart = tic;
            [X_FaLRTC, Out_SNN] = myHaLRTC(...
                X_Observed, ...                       % a tensor whose elements in Omega are used for estimating missing value
                Omega,...               % the index set indicating the obeserved elements
                alpha,...                  % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_*
                rho,...                      % the initial value of the parameter; it should be small enough
                opts.maxit,...               % the maximum iterations
                opts.tol,...                 % the tolerance of the relative difference of outputs of two neighbor iterations
                X_GT...
                );
            % calculate PSNR, SSIM and save results
            Time_all(ii,jj) = toc(tStart);
            Re_tensor{ii,jj} = X_FaLRTC;
            % quantitative metrics
%             MPSNRALL(ii,jj) = psnr(Re_tensor{ii,jj},X_GT);
%             SSIMALL(ii,jj) = ssim(uint8(Re_tensor{ii,jj}*255),uint8(X_GT*255));
            PPPSNRV = [];SSSIMV = [];
            for kkk = 1:size(X_Observed,3)
                PPPSNRV(kkk) = psnr(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk),1);
            end
            MPSNRALL(ii,jj) = mean(PPPSNRV);%psnr(Re_tensor{ii,jj},X_GT);
            for kkk = 1:size(X_Observed,3)
                SSSIMV(kkk) = ssim(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk));
            end
            SSIMALL(ii,jj) = mean(SSSIMV);
            fprintf(' %8.8s    %5.4s    %5.4s    %5.4s \n','method','PSNR', 'SSIM', 'time');
            fprintf(' %8.8s    %5.3f    %5.3f    %.1f \n',methodname{ii},MPSNRALL(ii,jj) , SSIMALL(ii,jj) , Time_all(ii,jj));
        end
        %% use TNN
        ii = 3;
        if EN_TNN
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 200;
            opts.Xtrue = X_GT;
            opts.beta = 1e-4;
            
            fprintf('\n');
            disp(['performing ',methodname{ii}, ' ... ']);
            
            
            tStart = tic;
            [X_TNN, Out_TNN] = LRTC_TNN(X_Observed, Omega, opts);
            Time_all(ii,jj) = toc(tStart);
            Re_tensor{ii,jj} = X_TNN;
            % quantitative metrics
%             MPSNRALL(ii,jj) = psnr(Re_tensor{ii,jj},X_GT);
%             SSIMALL(ii,jj) = ssim(uint8(Re_tensor{ii,jj}*255),uint8(X_GT*255));
            PPPSNRV = [];SSSIMV = [];
            for kkk = 1:size(X_Observed,3)
                PPPSNRV(kkk) = psnr(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk),1);
            end
            MPSNRALL(ii,jj) = mean(PPPSNRV);%psnr(Re_tensor{ii,jj},X_GT);
            for kkk = 1:size(X_Observed,3)
                SSSIMV(kkk) = ssim(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk));
            end
            SSIMALL(ii,jj) = mean(SSSIMV);
            fprintf(' %8.8s    %5.4s    %5.4s    %5.4s \n','method','PSNR', 'SSIM', 'time');
            fprintf(' %8.8s    %5.3f    %5.3f    %.1f \n',methodname{ii},MPSNRALL(ii,jj) , SSIMALL(ii,jj) , Time_all(ii,jj));
        end
        %% use EN_TNN_3DTV
        ii = 4;
        if EN_TNN_3DTV
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 100;
            opts.max_beta = 1e10;
            opts.rho = 1.05;
            for mu12 =[0.5]
                for mu3 = [0.1]
                    opts.beta =1e-2;
                    opts.lambda = 1e-4*[1 1 1];%paii(iii)*[1 1 1];
                    opts.lambda4 = 1e-2;
                    opts.mu= [mu12, mu12, mu3];%[pajj(jjj) pajj(jjj) papp(ppp)]
                    opts.max_lambda =1e10;
                    opts.max_lambda4 =1e10;
                    opts.Xtrue=X_GT;
                    
                    
                    fprintf('\n');
                    disp(['performing ',methodname{ii}, ' ... ']);
                    
                    tStart = tic;
                    [X_TNN_3DTV, Out_TNN_3Dtv] = LRTC_tnn3DTV_ADMM_2(X_Observed, Omega, opts);
                    Time_all(ii,jj) = toc(tStart);
                    Re_tensor{ii,jj} = X_TNN_3DTV;
                    % quantitative metrics
%                     MPSNRALL(ii,jj) = psnr(Re_tensor{ii,jj},X_GT);
%                     SSIMALL(ii,jj) = ssim(uint8(Re_tensor{ii,jj}*255),uint8(X_GT*255));
            PPPSNRV = [];SSSIMV = [];
            for kkk = 1:size(X_Observed,3)
                PPPSNRV(kkk) = psnr(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk),1);
            end
            MPSNRALL(ii,jj) = mean(PPPSNRV);%psnr(Re_tensor{ii,jj},X_GT);
            for kkk = 1:size(X_Observed,3)
                SSSIMV(kkk) = ssim(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk));
            end
            SSIMALL(ii,jj) = mean(SSSIMV);
                    fprintf(' %8.8s    %5.4s    %5.4s    %5.4s \n','method','PSNR', 'SSIM', 'time');
                    fprintf(' %8.8s    %5.3f    %5.3f    %.1f \n',methodname{ii},MPSNRALL(ii,jj) , SSIMALL(ii,jj) , Time_all(ii,jj));
                end
            end
        end
        %% use tSVD-FFDnet
        ii = 5;
        if EN_DP3LRTC
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 200;
            opts.Xtrue = X_GT;
            opts.debug = 1;
            
            opts.maxbeta = 1000;
            
            %%%%%
            disp(['performing ',methodname{ii}, ' ... ']);
            
            rho = 1;
            %                 for lambda = [1e-3]
            sigma = 0.1;
            opts.rho = rho;
            opts.sigma=sigma;
            beta = 1;
            opts.beta1=beta;
            opts.beta2=beta;
            
            tStart = tic;
            [X_DP3LRTC, Out_TSVD_FFDnet] = admm_DP3LR( X_Observed, Omega, opts );
            Time_all(ii,jj) = toc(tStart);
            Re_tensor{ii,jj} = X_DP3LRTC;
            % quantitative metrics
%             MPSNRALL(ii,jj) = psnr(Re_tensor{ii,jj},X_GT);
%             SSIMALL(ii,jj) = ssim(uint8(Re_tensor{ii,jj}*255),uint8(X_GT*255));
            PPPSNRV = [];SSSIMV = [];
            for kkk = 1:size(X_Observed,3)
                PPPSNRV(kkk) = psnr(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk),1);
            end
            MPSNRALL(ii,jj) = mean(PPPSNRV);%psnr(Re_tensor{ii,jj},X_GT);
            for kkk = 1:size(X_Observed,3)
                SSSIMV(kkk) = ssim(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk));
            end
            SSIMALL(ii,jj) = mean(SSSIMV);
            fprintf(' %8.8s    %5.4s    %5.4s    %5.4s \n','method','PSNR', 'SSIM', 'time');
            fprintf(' %8.8s    %5.3f    %5.3f    %.1f \n',methodname{ii},MPSNRALL(ii,jj) , SSIMALL(ii,jj) , Time_all(ii,jj));
        end
    end
    a1 = mean(MPSNRALL,2);
    a2 = mean(MPSNRALL,2);
    a3 = mean(Time_all,2);
    disp([' averaged metrics (Images, SR=' num2str(SR)]);
    for ii = 1:5
            fprintf(' %8.8s    %5.3f    %5.3f    %.1f \n',methodname{ii},a1(ii) , a2(ii) , a3(ii));
    end
    Curpath = pwd;
    result_name = [ Curpath '\results_new\MSI_result_sr_' num2str(100*SR) '.mat'];
    save(result_name,'Re_tensor','MPSNRALL','SSIMALL','Time_all','X_GT','X_Observed','Omega');
end