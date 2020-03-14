function createHist(X)
%CREATEFIGURE(data1, data2, data3)
%  DATA1:  histogram data
%  DATA2:  histogram data
%  DATA3:  histogram data

%  由 MATLAB 于 09-Apr-2019 14:22:56 自动生成
data1 = X(:,:,1);
data2 = X(:,:,2);
data3 = X(:,:,3);
% 创建 figure
figure1 = figure('Color',[1 1 1]);

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 创建 histogram
histogram(data1,'DisplayName','R channel','Parent',axes1,'EdgeColor','none',...
    'FaceColor',[1 0 0],...
    'BinMethod','auto');

% 创建 histogram
histogram(data2,'DisplayName','G channel','Parent',axes1,'LineWidth',2,...
    'EdgeColor','none',...
    'FaceColor',[0 1 0],...
    'BinMethod','auto');

% 创建 histogram
histogram(data3,'DisplayName','B channel','Parent',axes1,'EdgeColor','none',...
    'FaceColor',[0 0 1],...
    'BinMethod','auto');

% 创建 ylabel
ylabel('Number of pixels','FontWeight','bold','FontName','Times New Roman');

% 创建 xlabel
xlabel('Pixel values','FontWeight','bold','FontName','Times New Roman');



% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes1,[-255 255]);
box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'FontName','Times New Roman','FontWeight','bold','XTick',...
    [-150 -120 -90 -60 -30 0 30 60 90 120 150],'XTickLabel',...
    {'-1','-0.8','-0.6','-0.4','-0.2','0','0.2','0.4','0.6','0.8','1','',''});
% set(axes1,'FontName','Times New Roman','FontWeight','bold','XTickLabel',...
%     {'0','0.2','0.4','0.6','0.8','1'});
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northeast');
