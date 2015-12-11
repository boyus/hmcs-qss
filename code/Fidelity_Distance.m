% input data points
% output fidelity histogram, distance histogram, and fidelity vs. distance

clear all
close all
warning('off','all');

% choose one workspace to load based on your prior density

load('hmc_AA_2qb_primitive_1m_pts.mat');
% load('hmc_AA_2qb_Jeffreys_1m_pts.mat');
% load('hmc_AA_2qb_hedged_1m_pts.mat');

ln = length(rho);
randInd = zeros(1,2);
randMat = zeros(4);
fdl = zeros(1,ln);  % fidelity
dist = zeros(1,ln); % distance
dx = 100;  % number of bins
x_axis = 1/dx:1/dx:1;  % x axis of distribution plot
for i = 1:ln
    % random number between 1 and length of rho
	randInd = round(1+(ln-1).*rand(1,2));
    % randomly choose two matrices and multiply them together
    randMat = rho(:,:,randInd(1))*rho(:,:,randInd(2));
    
    randEig = sqrt(eig(randMat));  % square root of eigenvalues
    fdl(i) = real(sum(randEig));   % sum up to get fidelity
    
    randMat2 = rho(:,:,randInd(1))-rho(:,:,randInd(2));
    randMat2 = randMat2'*randMat2;
    dist(i) = trace(sqrtm(randMat2))/2; % trace distance
end

h = figure;
[f,x] = hist(fdl,dx);
bar(x_axis,f/ln);  %  fidelity distribution
axis([0,1,0,0.04]); xlabel('Fidelity');ylabel('Prior Density');
% title(['Fidelity Distribution (', priorType, ', 1m points)']);
fileName1 = strcat(fileName, '_fidelity');
set(gca,'xlim',[0 1]); set(gca,'FontSize',14);
print(h, '-djpeg', fileName1);

h = figure;
[d,x] = hist(dist,100);
bar(x_axis,d/ln);  %  distance distribution
axis([0,1,0,0.04]); xlabel('Distance');ylabel('Prior Density');
% title(['Distance Distribution (', priorType, ', 1m points)']);
fileName2 = strcat(fileName, '_distance');
set(gca,'xlim',[0 1]); set(gca,'FontSize',14);
print(h, '-djpeg', fileName2);

h = figure;
scatter(dist,fdl,'filled');  % fidelity vs. distance
axis([0,1,0,1]);xlabel('Distance');ylabel('Fidelity');
% title(['Fidelity vs. Distance (', priorType, ', 1m points)']);
fileName3 = strcat(fileName, '_fidelity_vs_distance');
set(gca,'FontSize',14);
print(h, '-djpeg', fileName3);

h = figure;
hist3([dist;fdl].',[50,50]);  % 3d histogram of fidelity and distance
xlabel('Distance'); ylabel('Fidelity');
% title(['Fidelity and Distance (', priorType, ', 1m points)']);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
set(gca,'FontSize',14);
fileName4 = strcat(fileName, '_fidelity_distance_3d_histogram');
print(h, '-djpeg', fileName4);
% fileName5 = strcat(fileName, '_fidelity_distance_2d');
% view(0,90);
% print(h, '-djpeg', fileName5);