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
Rse_TNN  = RSE(X_TNN(:),Xtrue(:));

SSIMvector=zeros(1,sz(3));
for s=1:1:sz(3)
    SSIMvector(s)=ssim(Xtrue(:,:,s)*256,X_TNN(:,:,s)*256);
end
SSIM_TNN = mean(SSIMvector);

display(sprintf('psnr=%.2f, ssim=%.4f, rse=%.3f, time=%.2f',psnr_TNN, SSIM_TNN, Rse_TNN, Out_TNN.time))
display(sprintf('=================================='))

lastname = [newname,'_SR_',num2str(SR),'_TNN_psnr_',num2str(psnr_TNN,'%.2f'),'_ssim_',num2str(SSIM_TNN,'%.4f'),'.mat'];
imname=fullfile(folderResultCur,newname,lastname);
save(imname,'X_TNN','opts','Out_TNN');

if printFig==1
    filename =[ 'F:\matlab\tSVD FFDnet\results\figs\',newname,'-TNN-',num2str(100*SR), '.jpg'];
    figure,
    imshow(X_TNN,'border','tight','initialmagnification','fit');
    set (gcf,'Position',[0,0,500,500]);
    axis normal;
    print(gcf,'-djpeg',filename);
end