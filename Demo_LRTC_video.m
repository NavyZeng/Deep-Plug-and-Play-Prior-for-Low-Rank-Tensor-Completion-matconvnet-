clc
clear; close all;
addpath(genpath(cd));
rng('default');


rng(213412);


methodname  = {'Observed','FaLRTC', 'TNN','TNN_3DTV','DP3LRTC'};
EN_FaLRTC   = 1;
EN_TNN     = 1;
EN_TNN_3DTV  = 1;
EN_DP3LRTC    = 1;
Mnum = length(methodname);


folderTest      = 'testsets';
folderResult    = 'results';
imageSets       = {'img','msi','video'};         % testing datasets
setTestCur      = imageSets{3};      % current testing dataset

Save_Images     = 1;
Show_Figure     = 1;
Enlarge_Factor  = 1;
%%
SR_case = 0;
for SR = [0.05 0.1 0.2]  % Smapling Rate
    fprintf('Sampling Rate = %d%%\n',SR*100);
    SR_case = SR_case+1;
    folderResultCur       =  fullfile(folderResult, [setTestCur,'_',num2str(SR)]);  % new folder for storing results
    if ~isfolder(folderResultCur)
        mkdir(folderResultCur)
    end
    ext         =  {'*.jpg','*.png','*.bmp','*.tif','*.mat'};
    filePaths   =  [];
    
    for i = 1 : length(ext) %length(ext)
        filePaths = cat(1,filePaths, dir(fullfile(folderTest,setTestCur,ext{i})));
    end
	Re_tensor	=  cell(Mnum,length(filePaths));
    MPSNRALL	=  zeros(Mnum,length(filePaths));
    SSIMALL     =  zeros(Mnum,length(filePaths));
    FSIMALL     =  zeros(Mnum,length(filePaths));
    Time_all    =  zeros(Mnum,length(filePaths));
    Image_name_all  = cell(length(filePaths),1);
    for jj = [1 7 4 5 3]%1:length(filePaths)
        
        % read images
        [~,Image_Name,extCur] = fileparts(filePaths(jj).name);
        disp(['Image #' num2str(jj) ': ' Image_Name ])
        load(filePaths(jj).name);
        %Input_Image = double(imread(fullfile(folderTest,setTestCur,filePaths(Image_Num).name)));
        Input_Image = X(:,:,1:30);
        Img_Size 	= size(Input_Image);
        Figure_Width       = Img_Size(2)*Enlarge_Factor;
        Figure_Height      = Img_Size(1)*Enlarge_Factor;
        
        X_GT        = zeros(Img_Size);
        if max(Input_Image(:))>2
            X_GT    = Input_Image/255;
        else
            X_GT   = Input_Image;
        end
        
        
        % Generate known data
        P           = round(SR*prod(Img_Size));      % prod·µ»Ø³Ë»ý
        Omega       = randsample(prod(Img_Size),P);
        [Omega,~]   = sort(Omega);   % Known = Omega
        % Missing data
        Known_Entries       = X_GT(Omega);
        
       %% Observed data
        ii = 1;
        disp([methodname{ii}, ':']);
        X_Observed          = zeros(Img_Size);
        X_Observed(Omega)   = Known_Entries;
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
        % save results
        
        
        %% FaLRTC
        ii = 2;
        if EN_FaLRTC
            %parameter setting
            opts = [ ];
            opts.tol   = 1e-5;
            opts.maxit = 200;
            alpha = [1, 1, 1];
            alpha = alpha/sum(alpha);
            rho = 1e-2;
            
            fprintf('\n');
            disp([methodname{ii}, ':']);
            
            tStart = tic;
            [X_FaLRTC, Out_FaLRTC] = myHaLRTC(...
                X_Observed, ...                       % a tensor whose elements in Omega are used for estimating missing value
                Omega,...               % the index set indicating the obeserved elements
                alpha,...                  % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_*
                rho,...                      % the initial value of the parameter; it should be small enough
                opts.maxit,...               % the maximum iterations
                opts.tol,...                 % the tolerance of the relative difference of outputs of two neighbor iterations
                X_GT...
                );
            % quantitative metrics
            Time_all(ii,jj) = toc(tStart);
            Re_tensor{ii,jj} = X_FaLRTC;
            % quantitative metrics
            PPPSNRV = [];SSSIMV = [];
            for kkk = 1:size(X_Observed,3)
                PPPSNRV(kkk) = psnr(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk),1);
            end
            MPSNRALL(ii,jj) = mean(PPPSNRV);%psnr(Re_tensor{ii,jj},X_GT);
            for kkk = 1:size(X_Observed,3)
                SSSIMV(kkk) = ssim(Re_tensor{ii,jj}(:,:,kkk),X_GT(:,:,kkk));
            end
            SSIMALL(ii,jj) = mean(SSSIMV);%ssim(uint8(Re_tensor{ii,jj}*255),uint8(X_GT*255));
            fprintf(' %8.8s    %5.4s    %5.4s    %5.4s \n','method','PSNR', 'SSIM', 'time');
            fprintf(' %8.8s    %5.3f    %5.3f    %.1f \n',methodname{ii},MPSNRALL(ii,jj) , SSIMALL(ii,jj) , Time_all(ii,jj));
        end
        %% TNN
        ii = 3;
        if EN_TNN
            %parameter setting
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 200;
            opts.Xtrue = X_GT;
            opts.beta = 1e-4;
%             opts.rho = 1.2;
            
            disp([methodname{ii}, ':']);
            tStart = tic;
            [X_TNN, Out_TNN] = LRTC_TNN(X_Observed, Omega, opts);
            % quantitative metrics
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
            %parameter setting
            opts = [ ];
            opts.tol   = 1e-5;
            opts.maxit = 100;
            opts.max_beta = 1e10;
            opts.rho = 1.05;
            for mu12 = [0.01]
                for mu3 = [0]
                    opts.beta =1e-2;
                    opts.lambda = 1e-4*[1 1 1];
                    opts.lambda4 = 1e-2;
                    opts.mu= [mu12, mu12, mu3];
                    opts.max_lambda =1e10;
                    opts.max_lambda4 =1e10;
                    opts.Xtrue=X_GT;
                    
                    disp([methodname{ii}, ':']);
                    tStart = tic;
                    [X_TNN_3DTV, Out_TNN_3DTV] = LRTC_tnn3DTV_ADMM_2(X_Observed, Omega, opts);
                    % quantitative metrics
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
            %parameter setting
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 200;
            opts.Xtrue = X_GT;
            opts.debug = 1;
            opts.sigma = 0.25;
            opts.beta1 = 2;
            opts.beta2 = 2;
            %opts.rho   = 1.01;
            %             opts.maxbeta = 5;
            
            %%%%%
            disp([methodname{ii}, ':']);
            tStart = tic;
            
            [X_DP3LRTC, Out_DP3LRTC] = admm_DP3LR( X_Observed, Omega, opts );
            % quantitative metrics
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
    result_name = [ Curpath '\results\video_result_sr_' num2str(100*SR) '.mat'];
    save(result_name,'Re_tensor','MPSNRALL','SSIMALL','Time_all','X_GT','X_Observed','Omega');
end
