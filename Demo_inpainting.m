clc;
clear; close all;
addpath(genpath(cd));

EN_SNN   = 1;
EN_TNN     = 1;
EN_TNN_3DTV  = 1;
EN_TSVDFFDnet    = 1;
methodname  = {'SNN', 'TNN','TNN_3DTV','DP3LRTC'};
Mnum = length(methodname);

folderTest  = 'testsets';
folderResult= 'results';
imageSets   = {'img','msi','inpainting'};         % testing datasets
setTestCur  = imageSets{3};      % current testing dataset


folderResultCur       =  fullfile(folderResult, setTestCur);
if ~isfolder(folderResultCur)
    mkdir(folderResultCur)
end

% read images
ext         =  {'*.jpg','*.png','*.bmp','*.tif'};
filePaths   =  [];
for i = 1 : length(ext)
    filePaths = cat(1,filePaths, dir(fullfile(folderTest,setTestCur,ext{i})));
end

for i = 1:length(filePaths)
    
    % load images
    disp([filePaths(i).name])
    label = double(imread(fullfile(folderTest,setTestCur,filePaths(i).name)));
    sz=size(label);
    
    [~,nameCur,extCur] = fileparts(filePaths(i).name);
    Xtrue = zeros(sz);
    
    % normalize
    if max(label)>1
        for s = 1:sz(3)
            maxX = max(max(label(:,:,s)));
            Xtrue(:,:,s) =  label(:,:,s)/maxX;
        end
    else
        Xtrue = label;
    end
    
    mask   = double(imread('mask3.png'));       % loading mask
    mask(  mask < 255 ) = 0;
    mask=mask/255;
    
    Xtemp = mask;
    
    Omega = find( Xtemp == 1 );                 % The support of the observed entries
    Xkn  = Xtrue(Omega);                        % Values of the observed entries
    Xmiss  = zeros(sz);
    Xmiss(Omega) = Xkn;                         % Observed data
    
    % save observed image
    newname = filePaths(i).name(1:end-4);
    mkdir(fullfile(folderResultCur,newname))
    lastname = [newname,'.mat'];
    imname=fullfile(folderResultCur,newname,lastname);
    save(imname,'Xtrue','Xmiss','sz');
    
    figName =[  'results\inpainting\',newname,'-mask5', '.jpg'];
    figure,
    imshow(Xmiss,'border','tight','initialmagnification','fit');
    set (gcf,'Position',[0,0,500,500]);
    axis normal;
    print(gcf,'-djpeg',figName);
    
    
    %% use SNN
    j = 1;
    if EN_SNN
        opts = [ ];
        opts.tol   = 1e-5;
        opts.maxit = 200;
        alpha = [1, 1, 1];
        alpha = alpha/sum(alpha);
        rho = 1e-2;
        
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
        
        fprintf('psnr=%.2f, ssim=%.4f,time=%.2f\n',psnr_SNN, SSIM_SNN,Out_SNN.time)
        fprintf('==================================\n')
        
        lastname = [newname,'_inpainting_','_SNN_psnr_',num2str(psnr_SNN,'%.2f'),'_ssim_',num2str(SSIM_SNN,'%.4f'),'.mat'];
        imname=fullfile(folderResultCur,newname,lastname);
        save(imname,'X_SNN','opts','Out_SNN');
        
        filename =[ 'results\inpainting\type5-',newname,'-SNN', '.jpg'];
        figure,
        imshow(X_SNN,'border','tight','initialmagnification','fit');
        set (gcf,'Position',[0,0,500,500]);
        axis normal;
        print(gcf,'-djpeg',filename);
        
    end
    %% use TNN
    j = j+1;
    if EN_TNN
        opts = [ ];
        opts.tol   = 1e-5;
        opts.maxit = 200;
        opts.Xtrue = Xtrue;
        
        
        fprintf('\n');
        disp(['performing ',methodname{j}, ' ... ']);
        for beta = [1e-4]
            opts.beta = beta;
            
            time_TNN = tic;
            [X_TNN, Out_TNN] = LRTC_TNN(Xmiss, Omega, opts);
            Out_TNN.time =toc(time_TNN);
            
            % calculate PSNR, SSIM and save results
            psnr_TNN = myPSNR(Xtrue,X_TNN);
            SSIM_TNN = mySSIM(Xtrue,X_TNN);
            
            fprintf('psnr=%.2f, ssim=%.4f, time=%.2f\n',psnr_TNN, SSIM_TNN, Out_TNN.time)
            fprintf('==================================\n')
            
            lastname = [newname,'_inpainting_','_TNN_psnr_',num2str(psnr_TNN,'%.2f'),'_ssim_',num2str(SSIM_TNN,'%.4f'),'.mat'];
            imname=fullfile(folderResultCur,newname,lastname);
            save(imname,'X_TNN','opts','Out_TNN');
            
            filename =[ 'results\inpainting\type5-',newname,'-TNN', '.jpg'];
            figure,
            imshow(X_TNN,'border','tight','initialmagnification','fit');
            set (gcf,'Position',[0,0,500,500]);
            axis normal;
            print(gcf,'-djpeg',filename);
        end
    end
    
    
    %% use TNN_3Dtv
    j = j+1;
    if EN_TNN_3DTV
        opts.tol   = 1e-5;
        opts.maxit = 100;
        opts.max_beta = 1e10;
        opts.rho = 1.05;
        for mu12 =[0.01]
            for mu3 = [0]
                opts.beta =1e-2;
                opts.lambda = 1e-4*[1 1 1];
                opts.lambda4 = 1e-2;
                opts.mu= [mu12, mu12, mu3];
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
                
                fprintf('psnr=%.2f, ssim=%.4f, time=%.2f\n',psnr_TNN_3Dtv, SSIM_TNN_3Dtv, Out_TNN_3Dtv.time)
                fprintf('==================================\n')
                
                lastname = [newname,'_inpainting_','_TNN_3Dtv_psnr_',num2str(psnr_TNN_3Dtv,'%.2f'),'_ssim_',num2str(SSIM_TNN_3Dtv,'%.4f'),'.mat'];
                imname=fullfile(folderResultCur,newname,lastname);
                save(imname,'X_TNN_3Dtv','opts','Out_TNN_3Dtv');
                
                filename =[ 'results\inpainting\type5-',newname,'-3d', '.jpg'];
                figure,
                imshow(X_TNN_3Dtv,'border','tight','initialmagnification','fit');
                set (gcf,'Position',[0,0,500,500]);
                axis normal;
                print(gcf,'-djpeg',filename);
            end
        end
    end
    
    %% use tSVD-FFDnet
    j = j+1;
    if EN_TSVDFFDnet
        opts = [ ];
        opts.tol   = 1e-4;
        opts.maxit = 50;
        opts.Xtrue = Xtrue;
        opts.debug = 1;
        
        beta = 0.5;
        sigma = 1;
        opts.beta1 = beta;
        opts.beta2 = beta;
        opts.sigma =  sigma; %  std2(Xtrue - Xmiss); %  0.1*std2(Xtrue - Xmiss)
        
        %%%%%
        disp(['performing ',methodname{j}, ' ... ']);
        
        time_FFDnet = tic;
        [X_TSVD_FFDnet, Out_TSVD_FFDnet] = admm_DP3LR( Xmiss, Omega, opts );
        Out_TSVD_FFDnet.time =toc(time_FFDnet);
        
        % calculate PSNR, SSIM and save results
        psnr_TSVD_FFDnet = myPSNR(Xtrue,X_TSVD_FFDnet);
        SSIM_TSVD_FFDnet = mySSIM(Xtrue,X_TSVD_FFDnet);
        
        fprintf('psnr=%.2f, ssim=%.4f, time=%.2f\n',psnr_TSVD_FFDnet, SSIM_TSVD_FFDnet, Out_TSVD_FFDnet.time)
        fprintf('beta=%.4f, sigma=%.3f\n',beta, sigma)
        fprintf('==================================\n')
        
        lastname = [newname,'_inpainting_','_TSVD_FFDnet_psnr_',num2str(psnr_TSVD_FFDnet,'%.2f'),'_ssim_',num2str(SSIM_TSVD_FFDnet,'%.4f'),'.mat'];
        imname=fullfile(folderResultCur,newname,lastname);
        save(imname,'X_TSVD_FFDnet','opts','Out_TSVD_FFDnet');
        
        %Save fig
        figName =[  'results\inpainting\type5-',newname,'-propose', '.jpg'];
        figure,
        imshow(X_TSVD_FFDnet,'border','tight','initialmagnification','fit');
        set (gcf,'Position',[0,0,500,500]);
        axis normal;
        print(gcf,'-djpeg',figName);
    end
end
