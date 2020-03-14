function [X,his] = LRTC_ffd(B, Omega, opts)
global sigmas

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts,'sigma');       sigma = opts.sigma;                  else  sigma = 0.05;                end
useGPU =1;
sigmas = sigma;
load(fullfile('FFDNet_Clip_color.mat'));
net = vl_simplenn_tidy(net);
if useGPU
    net = vl_simplenn_move(net, 'gpu') ;
end

[w, h, ~] = size(B);
X = B;

input = X;
input = single(input); %
if mod(w,2)==1
    input = cat(1,input, input(end,:,:)) ;
end
if mod(h,2)==1
    input = cat(2,input, input(:,end,:)) ;
end
if useGPU
    input = gpuArray(input);
end
max_in = max(input(:));min_in = min(input(:));
input = (input-min_in)/(max_in-min_in);
sigmas = sigma/(max_in-min_in);
res    = vl_simplenn(net,input,[],[],'conserveMemory',true,'mode','test');
output = res(end).x;
output(output<0)=0;output(output>1)=1;
output = output*(max_in-min_in)+min_in;
if mod(w,2)==1
    output = output(1:end-1,:,:);
end
if mod(h,2)==1
    output = output(:,1:end-1,:);
end
if useGPU
    output = gather(output);
end
X = double(output);
X(Omega) = B(Omega); 
his=[];
%fprintf('TNN ends: total iterations = %d   difference=%f\n', k, res(k));

end
