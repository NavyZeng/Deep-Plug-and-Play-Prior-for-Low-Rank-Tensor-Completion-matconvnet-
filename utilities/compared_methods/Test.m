Xtnn = double(imread('starfish-tnn-20.jpg'));
Xtnnffd = double(imread('starfish-tnn-ffd20.jpg'));
Xdp3 = double(imread('starfish-DP3-20.jpg'));
Xtrue = double(imread('starfish.jpg'));
Xmiss = double(imread('starfish-20.jpg'));
Err = Xtnn-Xtrue;
% createHist(Xtrue);
% createHist(Xmiss);
% createHist(Xtnn);
% createHist(Xtnnffd);
% createHist(Xdp3);
createHist(Err);
% Err2 = Xdp3-Xtrue;
% createHist(Err2);
% figure,
% imshow(Err2,[],'border','tight','initialmagnification','fit');
% set (gcf,'Position',[0,0,500,500]);
% axis normal;