% Script to carry out method on model data
% Example using SNR = 10 (use to reproduce figures B3, B4)
% ------------------------------------------------------ 
%% Overlap sparsity plots
%%% Model %%%
f = x;
[nrow, ncol] = size(f);
N = floor(log2(ncol));

[C,L] = wavedec(f(1,:),N,'coif3');

% using N = 10
[cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10] = detcoef(C,L,[1 2 3 4 5 6 7 8 9 10]);
detail_vec = zeros(1,N);

for level = 1:N    
    cd = eval(sprintf('cd\%.f',level));
    detail_vec(level) = max(abs(cd))/sum(abs(cd));     
end

% sparsity plot
close all
stem(detail_vec,'Color',[0, 0.4470, 0.7410], 'linewidth', 7)
ylim([0 0.5])
hold on
p1 = plot(detail_vec,'Color',[0, 0.4470, 0.7410], 'linewidth', 7, 'DisplayName', 'Signal');
hold on


%%% Model + Noise %%%
f = new_model_10;
[nrow, ncol] = size(f);
N = floor(log2(ncol));

[C,L] = wavedec(f(1,:),N,'coif3');

% using N = 10
[cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10] = detcoef(C,L,[1 2 3 4 5 6 7 8 9 10]);
detail_vec = zeros(1,N);

for level = 1:N    
    cd = eval(sprintf('cd\%.f',level));
    detail_vec(level) = max(abs(cd))/sum(abs(cd));     
end

% sparsity plot
stem(detail_vec,'Color',[0.8500, 0.3250, 0.0980], 'linewidth', 7)
ylim([0 0.5])
hold on
p2 = plot(detail_vec,'Color',[0.8500, 0.3250, 0.0980], 'linewidth', 7, 'DisplayName', 'Signal + Noise');
hold on

%%% Noise Only %%%
f = new_model_10 - x;
[nrow, ncol] = size(f);
N = floor(log2(ncol));

[C,L] = wavedec(f(1,:),N,'coif3');

% using N = 10
[cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10] = detcoef(C,L,[1 2 3 4 5 6 7 8 9 10]);
detail_vec = zeros(1,N);

for level = 1:N    
    cd = eval(sprintf('cd\%.f',level));
    detail_vec(level) = max(abs(cd))/sum(abs(cd));     
end


% sparsity plot
stem(detail_vec,'Color',[0.9290, 0.6940, 0.1250], 'linewidth', 7)
ylim([0 0.5])
xlim([0 7])
hold on
p3 = plot(detail_vec,'Color',[0.9290, 0.6940, 0.1250], 'linewidth', 7, 'DisplayName', 'Noise Only');

xticks(1:N)
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 22;

legend([p1 p2 p3], 'Location', 'NorthEast')
ax.Legend.FontSize = 14;

saveas(gcf, 'sparsity-snr10.tif');


%% Overlap change in sparsity plots
%%% Model %%%
f = x;
[nrow, ncol] = size(f);
N = floor(log2(ncol));

[C,L] = wavedec(f(1,:),N,'coif3');

% using N = 10
[cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10] = detcoef(C,L,[1 2 3 4 5 6 7 8 9 10]);
detail_vec = zeros(1,N);

for level = 1:N    
    cd = eval(sprintf('cd\%.f',level));
    detail_vec(level) = max(abs(cd))/sum(abs(cd));     
end

change_vec = zeros(1,N);
for level = 2:N    
    change_vec(level) = detail_vec(level)-detail_vec(level-1);
end

min_y1 = min(change_vec);

% change in sparsity plot
close all
stem(change_vec,'Color',[0, 0.4470, 0.7410], 'linewidth', 7)
hold on
p1 = plot(change_vec,'Color',[0, 0.4470, 0.7410], 'linewidth', 7, 'DisplayName', 'Signal');
hold on


%%% Model + Noise %%%
f = new_model_10;
[nrow, ncol] = size(f);
N = floor(log2(ncol));

[C,L] = wavedec(f(1,:),N,'coif3');

% using N = 10
[cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10] = detcoef(C,L,[1 2 3 4 5 6 7 8 9 10]);
detail_vec = zeros(1,N);

for level = 1:N    
    cd = eval(sprintf('cd\%.f',level));
    detail_vec(level) = max(abs(cd))/sum(abs(cd));     
end

change_vec = zeros(1,N);
for level = 2:N    
    change_vec(level) = detail_vec(level)-detail_vec(level-1);
end

min_y2 = min(change_vec);

% change in sparsity plot
stem(change_vec,'Color',[0.8500, 0.3250, 0.0980], 'linewidth', 7)
hold on
p2 = plot(change_vec,'Color',[0.8500, 0.3250, 0.0980], 'linewidth', 7, 'DisplayName', 'Signal + Noise');
hold on

%%% Noise Only %%%
f = new_model_10 - x;
[nrow, ncol] = size(f);
N = floor(log2(ncol));

%[cA,cD] = dwt(f(2,:),'coif3');
[C,L] = wavedec(f(1,:),N,'coif3');

% using N = 10
[cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10] = detcoef(C,L,[1 2 3 4 5 6 7 8 9 10]);
detail_vec = zeros(1,N);

for level = 1:N    
    cd = eval(sprintf('cd\%.f',level));
    detail_vec(level) = max(abs(cd))/sum(abs(cd));     
end

change_vec = zeros(1,N);
for level = 2:N    
    change_vec(level) = detail_vec(level)-detail_vec(level-1);
end

min_y3 = min(change_vec);

% change in sparsity plot
stem(change_vec,'Color',[0.9290, 0.6940, 0.1250], 'linewidth', 7)
ylim([-0.05 0.2])
xlim([0 7])
hold on
p3 = plot(change_vec,'Color',[0.9290, 0.6940, 0.1250], 'linewidth', 7, 'DisplayName', 'Noise Only');

xticks(1:N)
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 22;

legend([p1 p2 p3], 'Location', 'NorthEast')
ax.Legend.FontSize = 14;

saveas(gcf, 'change-snr10.tif');
