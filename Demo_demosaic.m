clc
clear; close all;
addpath(genpath(cd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EN_Demosaic  = 0;
EN_SNN   = 1;
EN_TNN     = 1;
EN_TNN_3DTV  = 1;
EN_TSVDFFDnet = 1;
methodname    = { 'Demosaic','SNN', 'TNN','TNN_3DTV','DP3LRTC'};
Mnum = length(methodname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderTest  = 'testsets';
folderResult= 'results';
imageSets   = {'img'};              %,'msi'}; % testing datasets
setTestCur  = imageSets{1};         % current testing dataset
printFig    = 1;
pauseTime   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderResultCur  =  fullfile(folderResult, 'demosaic');
if ~isfolder(folderResultCur)
    mkdir(folderResultCur)
end
%% read images
ext         =  {'*.jpg','*.png','*.bmp','*.tif'};
filePaths   =  [];
for i = 1 : length(ext)                       % length(ext)
    filePaths = cat(1,filePaths, dir(fullfile(folderTest,setTestCur,ext{i})));
end
for i = 3:8             % read images
    disp([filePaths(i).name])
    label = double(imread(fullfile(folderTest,setTestCur,filePaths(i).name)));
    [h,w,c]=size(label);
    sz = size(label);
    [~,nameCur,extCur] = fileparts(filePaths(i).name);
    Xtrue = zeros(h,w,c);
    
    if max(max(max(label)))>1
        for s = 1:c
            maxX = max(max(label(:,:,s)));
            Xtrue(:,:,s) =  label(:,:,s)/maxX; % πÈ“ªªØ
        end
    else
        Xtrue = label;
    end
    
    %% Missing data
    Xmiss = zeros(sz);
    [Xmos, Omega] = rgb2mosaic(Xtrue);
    Xmiss(Omega) = Xtrue(Omega);
    
    newname = filePaths(i).name(1:end-4);
    mkdir(fullfile(folderResultCur,newname))
    lastname = [newname,'.mat'];
    imname=fullfile(folderResultCur,newname,lastname);
    save(imname,'Xtrue','Xmos','Omega');
    
    filename =[ 'figs\',newname,'-mosaic', '.jpg'];
    figure,
    imshow(Xmos,'border','tight','initialmagnification','fit');
    set (gcf,'Position',[0,0,w,h]);
    axis normal;
    if printFig ==1
        print(gcf,'-djpeg',filename);
    end
    
    %% use Demosaic
    j = 1;
    if EN_Demosaic
        fprintf('\n');
        disp(['performing ',methodname{j}, ' ... ']);
        Xtemp = uint8(Xmos*255);
        
        X_Demosaic = double(demosaic(Xtemp,'rggb'))/255;
        
        
        % calculate PSNR, SSIM and save results
        psnr_Demosaic = myPSNR(Xtrue,X_Demosaic);
        %             Rse_SNN_FFDnet  = RSE(X_SNN(:),Xtrue(:));
        SSIM_Demosaic = mySSIM(Xtrue,X_Demosaic);
        
        fprintf('psnr=%.2f, ssim=%.4f \n',psnr_Demosaic, SSIM_Demosaic);
        fprintf('==================================\n');
        
        lastname = [newname,'_Demosaic_',num2str(psnr_Demosaic,'%.2f'),'_ssim_',num2str(SSIM_Demosaic,'%.4f'),'.mat'];
        imname=fullfile(folderResultCur,newname,lastname);
        save(imname,'X_Demosaic');
        if printFig==1
            filename =[ 'figs\',newname,'-Demosaic', '.jpg'];
            figure,
            imshow(X_Demosaic,'border','tight','initialmagnification','fit');
            set (gcf,'Position',[0,0,w,h]);
            axis normal;
            print(gcf,'-djpeg',filename);
        end
    end
    
     %% use SNN
        j = j+1;
        if EN_SNN
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 200;
            alpha = [1, 1, 1];
            alpha = alpha/sum(alpha);
            rho = 1;
            
            fprintf('\n');
            disp(['performing ',methodname{j}, ' ... ']);
            
            time_SNN = tic;
            [X_SNN, Out_SNN] = myHaLRTC(...
                Xmiss, ...                       % a tensor whose elements in Omega are used for estimating missing value
                Omega,...               % the index set indicating the obeserved elements
                alpha,...                  % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_*
                rho,...                      % the initial value of the parameter; it should be small enough
                opts.maxit,...               % the maximum iterations
                opts.tol,...                 % the tolerance of the relative difference of outputs of two neighbor iterations
                Xtrue...
                );
            Out_SNN.time =toc(time_SNN);
            
            % calculate PSNR, SSIM and save results
            psnr_SNN = myPSNR(Xtrue,X_SNN);
            SSIM_SNN = mySSIM(Xtrue,X_SNN);
            
            fprintf('psnr=%.2f, ssim=%.4f, time=%.2f\n',psnr_SNN, SSIM_SNN, Out_SNN.time)
            fprintf('==================================\n')
            
            lastname = [newname,'_SNN_psnr_',num2str(psnr_SNN,'%.2f'),'_ssim_',num2str(SSIM_SNN,'%.4f'),'.mat'];
            imname=fullfile(folderResultCur,newname,lastname);
            save(imname,'X_SNN','opts','Out_SNN');
            
            if printFig==1
                filename =['figs\',newname,'-snn-demos.jpg'];
                figure,
                imshow(X_SNN,'border','tight','initialmagnification','fit');
                set (gcf,'Position',[0,0,w,h]);
                axis normal;
                print(gcf,'-djpeg',filename);
            end
        end
        %% use TNN
        j = j+1;
        if EN_TNN
            opts = [ ];
            opts.tol   = 1e-5;
            opts.maxit = 200;
            opts.Xtrue = Xtrue;
            opts.beta = 1e-4;
            
            fprintf('\n');
            disp(['performing ',methodname{j}, ' ... ']);
            
            
            time_TNN = tic;
            [X_TNN, Out_TNN] = LRTC_TNN(Xmiss, Omega, opts);
            Out_TNN.time =toc(time_TNN);
            
            % calculate PSNR, SSIM and save results
            psnr_TNN = myPSNR(Xtrue,X_TNN);
            SSIM_TNN = mySSIM(Xtrue,X_TNN);
            
            fprintf('psnr=%.2f, ssim=%.4f, time=%.2f\n',psnr_TNN, SSIM_TNN, Out_TNN.time)
            fprintf('==================================')
            
            lastname = [newname,'_TNN_psnr_',num2str(psnr_TNN,'%.2f'),'_ssim_',num2str(SSIM_TNN,'%.4f'),'.mat'];
            imname=fullfile(folderResultCur,newname,lastname);
            save(imname,'X_TNN','opts','Out_TNN');
            
            if printFig==1
                filename =['figs\',newname,'-tnn-demos.jpg'];
                figure,
                imshow(X_TNN,'border','tight','initialmagnification','fit');
                set (gcf,'Position',[0,0,w,h]);
                axis normal;
                print(gcf,'-djpeg',filename);
            end
        end
        %% use EN_TNN_3DTV
        j = j+1;
        if EN_TNN_3DTV
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 400;
            opts.max_beta = 1e10;
            opts.rho = 1.05;
            for mu12 =[0.01]
                for mu3 = [0]
                    opts.beta =1e-2;
                    opts.lambda = 1e-4*[1 1 1];%paii(iii)*[1 1 1];
                    opts.lambda4 = 1e-2;
                    opts.mu= [mu12, mu12, mu3];%[pajj(jjj) pajj(jjj) papp(ppp)]
                    opts.max_lambda =1e10;
                    opts.max_lambda4 =1e10;
                    opts.Xtrue=Xtrue;
                    
                    
                    fprintf('\n');
                    disp(['performing ',methodname{j}, ' ... ']);
                    
                    time_3Dtv = tic;
                    [X_TNN_3Dtv, Out_TNN_3Dtv] = LRTC_tnn3DTV_ADMM_2(Xmiss, Omega, opts);
                    Out_TNN_3Dtv.time =toc(time_3Dtv);
                    
                    % calculate PSNR, SSIM and save results
                    psnr_TNN_3Dtv = myPSNR(Xtrue,X_TNN_3Dtv);
                    SSIM_TNN_3Dtv = mySSIM(Xtrue,X_TNN_3Dtv);
                    
                    fprintf('psnr=%.2f, ssim=%.4f, time=%.2f',psnr_TNN_3Dtv, SSIM_TNN_3Dtv, Out_TNN_3Dtv.time)
                    fprintf('mu12=%.5f, mu3=%.4f\n',opts.mu(1),opts.mu(3))
                    fprintf('==================================\n')
                    
                    lastname = [newname,'_TNN_3Dtv_psnr_',num2str(psnr_TNN_3Dtv,'%.2f'),'_ssim_',num2str(SSIM_TNN_3Dtv,'%.4f'),'.mat'];
                    imname=fullfile(folderResultCur,newname,lastname);
                    save(imname,'X_TNN_3Dtv','opts','Out_TNN_3Dtv');
                    
                    
                    if printFig ==1
                        filename =[ 'figs\',newname,'-3d-demos.jpg'];
                        figure,
                        imshow(X_TNN_3Dtv,'border','tight','initialmagnification','fit');
                        set (gcf,'Position',[0,0,w,h]);
                        axis normal;
                        print(gcf,'-djpeg',filename);
                    end
                end
            end
        end

    %% use tSVD-FFDnet
    j = j+1;
    if EN_TSVDFFDnet
        opts = [ ];
        opts.tol   = 1e-4;
        opts.maxit = 300;
        opts.Xtrue = Xtrue;
        opts.debug = 0;
        %             opts.rho = 1;
        %             opts.maxbeta = 1;
        for beta = [2]
            opts.beta1 = beta;
            opts.beta2 = beta;
            for sigma = [0.1]
                opts.sigma =  sigma;
                opts.max_beta = 1e10;
                disp(['performing ',methodname{j}, ' ... ']);
                time_FFDnet = tic;
                [X_TSVD_FFDnet, Out_TSVD_FFDnet] = admm_DP3LR( Xmiss, Omega, opts );
                Out_TSVD_FFDnet.time =toc(time_FFDnet);
                % calculate PSNR, SSIM and save results
                psnr_TSVD_FFDnet = myPSNR(Xtrue,X_TSVD_FFDnet);
                %             Rse_TSVD_FFDnet  = RSE(X_TSVD_FFDnet(:),Xtrue(:));
                SSIM_TSVD_FFDnet = mySSIM(Xtrue,X_TSVD_FFDnet);
                
                fprintf('psnr=%.2f, ssim=%.4f, time=%.2f\n',psnr_TSVD_FFDnet, SSIM_TSVD_FFDnet, Out_TSVD_FFDnet.time)
                fprintf('beta=%.3f, sigma=%.3f\n',beta,opts.sigma)
                fprintf('==================================\n')
                lastname = [newname,'_demos_','TSVD_FFDnet_psnr_',num2str(psnr_TSVD_FFDnet,'%.2f'),'_ssim_',num2str(SSIM_TSVD_FFDnet,'%.4f'),'.mat'];
                imname=fullfile(folderResultCur,newname,lastname);
                save(imname,'X_TSVD_FFDnet','opts','Out_TSVD_FFDnet');
                figure,
                imshow(X_TSVD_FFDnet,'border','tight','initialmagnification','fit');
                if printFig==1
                    %Save fig
                    filename =[ 'figs\',newname,'-dp3-demos.jpg'];
                    
                    set (gcf,'Position',[0,0,w,h]);
                    axis normal;
                    print(gcf,'-djpeg',filename);
                end
            end
        end
    end
end