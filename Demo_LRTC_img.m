clc
clear; close all;
addpath(genpath(cd));
rng('default');


Method_name  = {'Observed','FaLRTC', 'TNN','TNN_3DTV','DP3LRTC'};
EN_FaLRTC   = 1;
EN_TNN     = 1;
EN_TNN_3DTV  = 1;
EN_DP3LRTC    = 1;

Mnum = length(Method_name);

folderTest      = 'testsets';
folderResult    = 'results';
imageSets       = {'img','msi'};         % testing datasets
setTestCur      = imageSets{1};      % current testing dataset

Save_Images     = 1;
Show_Figure     = 1;
Enlarge_Factor  = 1;
%%
SR_case = 0;
for SR = [0.1 0.2 0.3]  % Smapling Rate
    fprintf('Sampling Rate = %d%%\n',SR*100);
    SR_case = SR_case+1;
    folderResultCur       =  fullfile(folderResult, [setTestCur,'_',num2str(SR)]);  % new folder for storing results
    if ~isfolder(folderResultCur)
        mkdir(folderResultCur)
    end
    ext         =  {'*.jpg','*.png','*.bmp','*.tif'};
    filePaths   =  [];
    
    for i = 1 : length(ext) %length(ext)
        filePaths = cat(1,filePaths, dir(fullfile(folderTest,setTestCur,ext{i})));
    end
    
    for Image_Num = 1:length(filePaths)
        
        % read images
        [~,Image_Name,extCur] = fileparts(filePaths(Image_Num).name);
        disp(['Image #' num2str(Image_Num) ': ' Image_Name ])
        Input_Image = double(imread(fullfile(folderTest,setTestCur,filePaths(Image_Num).name)));
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
        Method_num = 1;
        disp([Method_name{Method_num}, ':']);
        X_Observed          = zeros(Img_Size);
        X_Observed(Omega)   = Known_Entries;
        % quantitative metrics
        PSNR_Observed = psnr(X_Observed,X_GT);PSNR_all(Method_num) = PSNR_Observed;
        SSIM_Observed = ssim(X_Observed*255,X_GT*255);SSIM_all(Method_num) = SSIM_Observed;
        Time_Observed = 0;Time_all(Method_num) = Time_Observed;
        fprintf('PSNR = %.2f, SSIM =%.4f, Time=%.2f\n',PSNR_all(Method_num), SSIM_all(Method_num), Time_all(Method_num))
        % save results
        mkdir(fullfile(folderResultCur,Image_Name))
        Results_name = [Image_Name,'.mat'];
        Path_Results_name=fullfile(folderResultCur,Image_Name,Results_name);
        save(Path_Results_name,'X_GT','X_Observed','Omega');
        % save images & show figures
        if Save_Images
            filename =[ pwd '\figs\' setTestCur '\' Image_Name '_SR_' num2str(100*SR)  '_' Method_name{Method_num} '.jpg'];
            imwrite(X_Observed,filename);
        end
        if Show_Figure
            figure(Method_num);
            imshow(X_Observed,'border','tight','initialmagnification','fit');
            set (gcf,'Position',[50,50,50+Figure_Width,50+Figure_Height]);
            title(['Image: ' Image_Name '; SR=' num2str(100*SR) '; ' Method_name{Method_num}]);
        end
        %% FaLRTC
        Method_num = 2;
        if EN_FaLRTC
            %parameter setting
            opts = [ ];
            opts.tol   = 1e-5;
            opts.maxit = 200;
            alpha = [1, 1, 1];
            alpha = alpha/sum(alpha);
            rho = 1e-2;
            
            fprintf('\n');
            disp([Method_name{Method_num}, ':']);
            
            t_Start = tic;
            [X_FaLRTC, Out_FaLRTC] = myHaLRTC(...
                X_Observed, ...                       % a tensor whose elements in Omega are used for estimating missing value
                Omega,...               % the index set indicating the obeserved elements
                alpha,...                  % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_*
                rho,...                      % the initial value of the parameter; it should be small enough
                opts.maxit,...               % the maximum iterations
                opts.tol,...                 % the tolerance of the relative difference of outputs of two neighbor iterations
                X_GT...
                );
            Time_FaLRTC =toc(t_Start);Time_all(Method_num) = Time_FaLRTC;
            % quantitative metrics
            PSNR_FaLRTC = psnr(X_FaLRTC,X_GT);PSNR_all(Method_num) = PSNR_FaLRTC;
            SSIM_FaLRTC = ssim(X_FaLRTC*255,X_GT*255);SSIM_all(Method_num) = SSIM_FaLRTC;
            
            fprintf('PSNR = %.2f, SSIM =%.4f, Time=%.2f\n',PSNR_all(Method_num), SSIM_all(Method_num), Time_all(Method_num))
            fprintf('==================================\n')
            % save results
            Results_name = [Image_Name,'_SR_',num2str(SR),'_SNN_psnr_',num2str(PSNR_FaLRTC,'%.2f'),'_ssim_',num2str(SSIM_FaLRTC,'%.4f'),'.mat'];
            Path_Results_name=fullfile(folderResultCur,Image_Name,Results_name);
            save(Path_Results_name,'X_FaLRTC','opts','Out_FaLRTC');
            if Save_Images
                filename =[ pwd '\figs\' setTestCur '\' Image_Name '_SR_' num2str(100*SR)  '_' Method_name{Method_num} '.jpg'];
                imwrite(X_FaLRTC,filename);
            end
            if Show_Figure
                figure(Method_num);
                imshow(X_FaLRTC,'border','tight','initialmagnification','fit');
                set (gcf,'Position',[50,50,50+Figure_Width,50+Figure_Height]);
                title(['Image: ' Image_Name '; SR=' num2str(100*SR) '; ' Method_name{Method_num}]);
            end
        end
        %% TNN
        Method_num = 3;
        if EN_TNN
            %parameter setting
            opts = [ ];
            opts.tol   = 1e-5;
            opts.maxit = 200;
            opts.Xtrue = X_GT;
            opts.beta = 1e-4;
            opts.rho = 1.1;
            
            disp([Method_name{Method_num}, ':']);
            t_Start = tic;
            [X_TNN, Out_TNN] = LRTC_TNN(X_Observed, Omega, opts);
            Time_TNN =toc(t_Start);Time_all(Method_num) = Time_TNN;
            % calculate PSNR, SSIM and save results
            PSNR_TNN = psnr(X_TNN,X_GT);PSNR_all(Method_num) = PSNR_TNN;
            SSIM_TNN = ssim(X_TNN*255,X_GT*255);SSIM_all(Method_num) = SSIM_TNN;
            
            fprintf('PSNR = %.2f, SSIM =%.4f, Time=%.2f\n',PSNR_all(Method_num), SSIM_all(Method_num), Time_all(Method_num))
            fprintf('==================================\n')
            % save results
            Results_name = [Image_Name,'_SR_',num2str(SR),'_TNN_psnr_',num2str(PSNR_TNN,'%.2f'),'_ssim_',num2str(SSIM_TNN,'%.4f'),'.mat'];
            Path_Results_name=fullfile(folderResultCur,Image_Name,Results_name);
            save(Path_Results_name,'X_TNN','opts','Out_TNN');
            
            if Save_Images
                filename =[ pwd '\figs\' setTestCur '\' Image_Name '_SR_' num2str(100*SR)  '_' Method_name{Method_num} '.jpg'];
                imwrite(X_TNN,filename);
            end
            if Show_Figure
                figure(Method_num);
                imshow(X_TNN,'border','tight','initialmagnification','fit');
                set (gcf,'Position',[50,50,50+Figure_Width,50+Figure_Height]);
                title(['Image: ' Image_Name '; SR=' num2str(100*SR) '; ' Method_name{Method_num}]);
            end
        end
        %% use EN_TNN_3DTV
        Method_num = 4;
        if EN_TNN_3DTV
            %parameter setting
            opts = [ ];
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
                    opts.Xtrue=X_GT;
                    
                    disp([Method_name{Method_num}, ':']);
                    t_Start = tic;
                    [X_TNN_3DTV, Out_TNN_3DTV] = LRTC_tnn3DTV_ADMM_2(X_Observed, Omega, opts);
                    Time_TNN_3DTV =toc(t_Start);Time_all(Method_num) = Time_TNN_3DTV;
                    % calculate PSNR, SSIM and save results
                    PSNR_TNN_3DTV = psnr(X_TNN_3DTV,X_GT);PSNR_all(Method_num) = PSNR_TNN_3DTV;
                    SSIM_TNN_3DTV = ssim(X_TNN_3DTV*255,X_GT*255);SSIM_all(Method_num) = SSIM_TNN_3DTV;
                    
                    fprintf('PSNR = %.2f, SSIM =%.4f, Time=%.2f\n',PSNR_all(Method_num), SSIM_all(Method_num), Time_all(Method_num))
                    fprintf('==================================\n')
                    % save results
                    Results_name = [Image_Name,'_SR_',num2str(SR),'_TNN3DTV_psnr_',num2str(PSNR_TNN_3DTV,'%.2f'),'_ssim_',num2str(SSIM_TNN_3DTV,'%.4f'),'.mat'];
                    Path_Results_name=fullfile(folderResultCur,Image_Name,Results_name);
                    save(Path_Results_name,'X_TNN_3DTV','opts','Out_TNN_3DTV');
                    
                    if Save_Images
                        filename =[ pwd '\figs\' setTestCur '\' Image_Name '_SR_' num2str(100*SR)  '_' Method_name{Method_num} '.jpg'];
                        imwrite(X_TNN_3DTV,filename);
                    end
                    if Show_Figure
                        figure(Method_num);
                        imshow(X_TNN_3DTV,'border','tight','initialmagnification','fit');
                        set (gcf,'Position',[50,50,50+Figure_Width,50+Figure_Height]);
                        title(['Image: ' Image_Name '; SR=' num2str(100*SR) '; ' Method_name{Method_num}]);
                    end
                end
            end
        end
        %% use tSVD-FFDnet
        Method_num = 5;
        if EN_DP3LRTC
            %parameter setting
            opts = [ ];
            opts.tol   = 1e-4;
            opts.maxit = 200;
            opts.Xtrue = X_GT;
            opts.debug = 1;
            opts.sigma = 0.2;
            opts.beta1 = 1.0;
            opts.beta2 = 1.0;
            opts.rho   = 1.01;
            %             opts.maxbeta = 5;
            
            %%%%%
            disp([Method_name{Method_num}, ':']);
            t_Start = tic;
            
            [X_DP3LRTC, Out_DP3LRTC] = admm_DP3LR( X_Observed, Omega, opts );
            Time_DP3LRTC =toc(t_Start);Time_all(Method_num) = Time_DP3LRTC;
            %X_TNN_FFDnet = myNormalization(X_TNN_FFDnet);
            % calculate PSNR, SSIM and save results
            PSNR_DP3LRTC = psnr(X_DP3LRTC,X_GT);PSNR_all(Method_num) = PSNR_DP3LRTC;
            SSIM_DP3LRTC = ssim(X_DP3LRTC*255,X_GT*255);SSIM_all(Method_num) = SSIM_DP3LRTC;
            
            fprintf('PSNR = %.2f, SSIM =%.4f, Time=%.2f\n',PSNR_all(Method_num), SSIM_all(Method_num), Time_all(Method_num))
            fprintf('==================================\n')
            % save results
            Results_name = [Image_Name,'_SR_',num2str(SR),'_DP3LRTC_',num2str(PSNR_DP3LRTC,'%.2f'),'_ssim_',num2str(SSIM_DP3LRTC,'%.4f'),'.mat'];
            Path_Results_name=fullfile(folderResultCur,Image_Name,Results_name);
            save(Path_Results_name,'X_DP3LRTC','opts','Out_DP3LRTC');
            
            if Save_Images
                filename =[ pwd '\figs\' setTestCur '\' Image_Name '_SR_' num2str(100*SR)  '_' Method_name{Method_num} '.jpg'];
                imwrite(X_DP3LRTC,filename);
            end
            if Show_Figure
                figure(Method_num);
                imshow(X_DP3LRTC,'border','tight','initialmagnification','fit');
                set (gcf,'Position',[50,50,50+Figure_Width,50+Figure_Height]);
                title(['Image: ' Image_Name '; SR=' num2str(100*SR) '; ' Method_name{Method_num}]);
            end
        end
        
    end
end