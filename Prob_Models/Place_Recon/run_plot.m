function run_plot
% Run reconstruction
[loc, recon, post] = place_recon;
n = length(loc);
dt = 1; % Time in seconds
tail = 10;

% Smoothing kernel
sigma = (3 / 2*sqrt(2*log(2)))/.5;
G_kern = fspecial('gaussian',[7*ceil(sigma) 7*ceil(sigma)], ...
                  sigma);

map = [linspace(0,.7,64)' zeros(64,1) zeros(64,1)];
colormap(map);

% Make the video
figure(1);
xlim([1 70]);
ylim([1 70]);
set(gca,'Color',[0 0 0]);
for i = 1:n
    cla
    hold on
    tmp = post{i};
    tmp(isnan(tmp)) = 0;
    posterior = imfilter(tmp,G_kern);
    imagesc(posterior);
    
    [~,argmax_idx] = max(post{i}(:));
    [x y] = ind2sub(size(post{i}),argmax_idx);
    if i <= tail
        
        plot(loc(1:i,1),loc(1:i,2),'-g');
        scatter(loc(i,1),loc(i,2),'g','filled');
        
        plot(recon(1:i,2),recon(1:i,1),'-r');
        scatter(recon(i,2),recon(i,1),'r','filled');
        
    else
        plot(loc(i-tail:i,1),loc(i-tail:i,2),'-g');
        scatter(loc(i,1),loc(i,2),'g','filled');
        
        plot(recon(i-tail:i,2),recon(i-tail:i,1),'-r');
        scatter(recon(i,2),recon(i,1),'r','filled');
    end
    hold off
    pause(dt/10);
end






