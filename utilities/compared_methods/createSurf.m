function createfigure(zdata1)
%CREATEFIGURE(zdata1)
%  ZDATA1:  surface zdata

%  �� MATLAB �� 11-Apr-2019 11:16:34 �Զ�����

% ���� figure
figure1 = figure('Color',[1 1 1]);

% ���� axes
axes1 = axes('Parent',figure1,...
    'Position',[0.0889674078094611 0.0722641509433965 0.765302143304838 0.907245361839955]);
hold(axes1,'on');

% ���� surf
surf(zdata1,'Parent',axes1);
shading interp;
% ���� zlabel
zlabel('PSNR (dB)');

% ���� ylabel
ylabel('\beta');

% ���� xlabel
xlabel('\sigma');

view(axes1,[-37.5 30]);
grid(axes1,'on');
% ������������������
set(axes1,'FontName','Times New Roman','FontSize',12,'FontWeight','bold',...
    'XDir','reverse','XMinorTick','on','XScale','log','XTick',[1 2 3 4 5 6],...
    'XTickLabel',{'10','1','10^{-1}','10^{-2}','10^{-3}','10^{-4}'},'YDir',...
    'reverse','YMinorTick','on','YScale','log','YTickLabel',...
    {'10','1','10^{-1}','10^{-2}','10^{-3}','10^{-4}'});
% ���� colorbar
colorbar('peer',axes1);

