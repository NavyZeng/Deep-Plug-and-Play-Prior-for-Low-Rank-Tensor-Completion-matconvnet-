function createHist(X)
%CREATEFIGURE(data1, data2, data3)
%  DATA1:  histogram data
%  DATA2:  histogram data
%  DATA3:  histogram data

%  �� MATLAB �� 09-Apr-2019 14:22:56 �Զ�����
data1 = X(:,:,1);
data2 = X(:,:,2);
data3 = X(:,:,3);
% ���� figure
figure1 = figure('Color',[1 1 1]);

% ���� axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% ���� histogram
histogram(data1,'DisplayName','R channel','Parent',axes1,'EdgeColor','none',...
    'FaceColor',[1 0 0],...
    'BinMethod','auto');

% ���� histogram
histogram(data2,'DisplayName','G channel','Parent',axes1,'LineWidth',2,...
    'EdgeColor','none',...
    'FaceColor',[0 1 0],...
    'BinMethod','auto');

% ���� histogram
histogram(data3,'DisplayName','B channel','Parent',axes1,'EdgeColor','none',...
    'FaceColor',[0 0 1],...
    'BinMethod','auto');

% ���� ylabel
ylabel('Number of pixels','FontWeight','bold','FontName','Times New Roman');

% ���� xlabel
xlabel('Pixel values','FontWeight','bold','FontName','Times New Roman');



% ȡ�������е�ע���Ա����������� X ��Χ
xlim(axes1,[-255 255]);
box(axes1,'on');
% ������������������
set(axes1,'FontName','Times New Roman','FontWeight','bold','XTick',...
    [-150 -120 -90 -60 -30 0 30 60 90 120 150],'XTickLabel',...
    {'-1','-0.8','-0.6','-0.4','-0.2','0','0.2','0.4','0.6','0.8','1','',''});
% set(axes1,'FontName','Times New Roman','FontWeight','bold','XTickLabel',...
%     {'0','0.2','0.4','0.6','0.8','1'});
% ���� legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northeast');
