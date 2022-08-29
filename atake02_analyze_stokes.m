clear all; close all; 
b_use_averaged = false;

folder_processed = fullfile(pwd,'example_data','[p.Chicken052621][s.SANsChicken_smallV][05-26-2021_13-11-59]', 'Processed');
ind1 = findstr(folder_processed, '[p.');
ind2 = findstr(folder_processed, 'Processed')-2;
fileID = folder_processed (ind1:ind2) %'[p.Chicken052621][s.SANsChicken_smallV][05-26-2021_13-11-59]';
frame = 6;
if b_use_averaged % use original python output (Damon's averaged version)
    s1 = read_stokes(folder_processed, fileID, 'sv1', frame);
    size(s1)
    s2 = read_stokes(folder_processed, fileID, 'sv2', frame);
    alinepos_scale = 1;
else % use purely raw stokes 
    s1 = read_stokes(folder_processed, fileID, 'iquv1', frame);
    size(s1)
    s2 = read_stokes(folder_processed, fileID, 'iquv2', frame);
    alinepos_scale = 5/2;
end

% averaging kernel
dimIn = size(s1);
fwz = 9; 
fwx = 21;
h = construct_kernel(fwz, fwx);

%% 
% s = [i;q;u;v] before averaging
% S = [I;Q;U;V] after averaging

S1 = conv2(s1(:,:),h,'same'); %<<now this is averaging along alines>> 
S2 = conv2(s2(:,:),h,'same');
S1 = reshape(S1, dimIn);
S2 = reshape(S2, dimIn);

% L = squared Euclidian length of Q,U,V
L1 = (sum(S1(:,:,:,2:4).^2,4));
L2 = (sum(S2(:,:,:,2:4).^2,4));

If = mean(S1(:,:,:,1).^2 + S2(:,:,:,1).^2,3); % keep filtered intensity as optional output argument

% computation of DOP
dop = sqrt(mean(L1 + L2,3)./If);
If = sqrt(If);

% final normalization
QUV1 = bsxfun(@rdivide,S1(:,:,:,2:4),sqrt(L1));
QUV2 = bsxfun(@rdivide,S2(:,:,:,2:4),sqrt(L2));
I1 = S1(:,:,:,1);
I2 = S2(:,:,:,1);
quv1 = s1(:,:,:,2:4)./s1(:,:,:,1);
quv2 = s2(:,:,:,2:4)./s2(:,:,:,1);
% clear L1 L2;


%%
qq_bin = 3;
ii1 = S1(:,:,qq_bin,1);
ii2 = S2(:,:,qq_bin,1);

T1 = I1./sqrt(L1);
T2 = I2./sqrt(L2);
a1 = squeeze(T1(:,:,3)); %this is very close to 1

%%

fpos = [550, 0, 560, 1100];

figure(1);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(10*log10(s1(:,:,qq_bin,1))); colorbar; title('i');
subplot(4,1,2); imagesc(quv1(:,:,qq_bin,1), [-1,1]); colorbar; title('q');
subplot(4,1,3); imagesc(quv1(:,:,qq_bin,2), [-1,1]); colorbar; title('u');
subplot(4,1,4); imagesc(quv1(:,:,qq_bin,3), [-1,1]); colorbar; title('v');
sgtitle('(unaveraged) s1=iquv1') 

figure(2);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(10*log10(s2(:,:,qq_bin,1))); colorbar; title('i');
subplot(4,1,2); imagesc(quv2(:,:,qq_bin,1), [-1,1]); colorbar; title('q');
subplot(4,1,3); imagesc(quv2(:,:,qq_bin,2), [-1,1]); colorbar; title('u');
subplot(4,1,4); imagesc(quv2(:,:,qq_bin,3), [-1,1]); colorbar; title('v');
sgtitle('(unaveraged) s2=iquv2') 

figure(11);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(T1(:,:,qq_bin), [0,2]); colorbar; title('t');
subplot(4,1,2); imagesc(QUV1(:,:,qq_bin,1), [-1,1]); colorbar; title('Q');
subplot(4,1,3); imagesc(QUV1(:,:,qq_bin,2), [-1,1]); colorbar; title('U');
subplot(4,1,4); imagesc(QUV1(:,:,qq_bin,3), [-1,1]); colorbar; title('V');
sgtitle('(averaged) QUV1') 

figure(12);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(T2(:,:,qq_bin), [0,2]); colorbar; title('t');
subplot(4,1,2); imagesc(QUV2(:,:,qq_bin,1), [-1,1]); colorbar; title('Q');
subplot(4,1,3); imagesc(QUV2(:,:,qq_bin,2), [-1,1]); colorbar; title('U');
subplot(4,1,4); imagesc(QUV2(:,:,qq_bin,3), [-1,1]); colorbar; title('V');
sgtitle('(averaged) QUV2') 

%% Analyze Stokes (QUV) progression into depth

select_zpos = [550:700]+450;
alinepos = floor(100*alinepos_scale);
nfig = 100;
% select_zpos = [550:700]+480;
% alinepos = floor(800*alinepos_scale);
% nfig = 200;


figure(nfig);
imagesc(log10(ii1)); hold on;
plot(ones(size(select_zpos))*alinepos, select_zpos, 'k'); 
plot(alinepos, select_zpos(1),'kx');
hold off; 
title(fileID, 'Interpreter', 'none');

select_nsv1 = squeeze(QUV1(select_zpos, alinepos, qq_bin, 1:3));
select_nsv2 = squeeze(QUV2(select_zpos, alinepos, qq_bin, 1:3));

hfig = figure(nfig+1); 
[x,y,z] = sphere(16);
hsph = surf(x,y,z);  set(hsph, 'FaceAlpha',0.2, 'LineStyle', 'none'); 
colormap([1 1 0]); hold on;
scatter3(select_nsv1(1,1),select_nsv1(1,2),select_nsv1(1,3),55,'xr'); 
% plot3(select_nsv1(1:end,1),select_nsv1(1:end,2),select_nsv1(1:end,3),'.-r'); 
scatter3(select_nsv1(:,1),select_nsv1(:,2),select_nsv1(:,3),15,autumn(length(select_nsv1)),'filled'); 

scatter3(select_nsv2(1,1),select_nsv2(1,2),select_nsv2(1,3),55, 'xb'); 
% plot3(select_nsv2(1:end,1),select_nsv2(1:end,2),select_nsv2(1:end,3),'.-b'); 
scatter3(select_nsv2(:,1),select_nsv2(:,2),select_nsv2(:,3),15,winter(length(select_nsv2)),'filled'); 
hold off;
title(fileID, 'Interpreter', 'none');

%%
function data = read_stokes(folder_processed, fildID, fkey, frame)
% fkey = 'sv1';
ffile_h5 = fullfile(folder_processed, [fildID, sprintf('.%s.h5',fkey)]);

% h5disp(ffile_h5)
h5_start_vec = [frame,1,1,1,1]; %[frame, bin, stokes, z_ind, x_ind] 
data = h5read(ffile_h5, sprintf('/%s',fkey), h5_start_vec,[1, Inf, Inf, Inf, Inf]);
% size(data)
data = permute(data, [4,5,2,3,1]); 
data = squeeze(data);
% size(data)
end

function kernel = construct_kernel(fwz, fwx)
    % fwx = 15;% filter width x
    nx = (round(fwx*1.5)-1)/2;
    nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
    hx = exp(-nx.^2);
    hx = hx(:)/sum(hx(:));
    % fwz = 1;% filter width z
    nz = (round(fwz*1.5)-1)/2;
    nz = linspace(-nz,nz,round(fwz*1.5))*2*sqrt(log(2))/fwz;
    hz = exp(-nz.^2);
    hz = hz(:)/sum(hz(:));
    kernel = hz*hx'; 
    % figure; imagesc(kernel)
end