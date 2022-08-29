
b_use_averaged = true;

folder_processed = 'D:\OFDIData\user.Ilyas\[p.S3B_11_25_19][s._S_loc3][11-25-2019_12-38-31]\Processed'
fildID = '[p.S3B_11_25_19][s._S_loc3][11-25-2019_12-38-31]';
 
if b_use_averaged % use original python output (Damon's averaged version)
    S1 = read_stokes(folder_processed, fildID, 'sv1');
    size(S1)
    S2 = read_stokes(folder_processed, fildID, 'sv2');
    alinepos_scale = 1;
else % use purely raw stokes 
    S1 = read_stokes(folder_processed, fildID, 'iquv1');
    size(S1)
    S2 = read_stokes(folder_processed, fildID, 'iquv2');
    alinepos_scale = 5/2;
end
% Euclidian length of Q,U,V
L1 = (sum(S1(:,:,:,2:4).^2,4));
L2 = (sum(S2(:,:,:,2:4).^2,4));

If = mean(S1(:,:,:,1).^2 + S2(:,:,:,1).^2,3); % keep filtered intensity as optional output argument

% computation of DOP
dop = sqrt(mean(L1 + L2,3)./If);
If = sqrt(If);

% final normalization
QUV1 = bsxfun(@rdivide,S1(:,:,:,2:4),sqrt(L1));
QUV2 = bsxfun(@rdivide,S2(:,:,:,2:4),sqrt(L2));

% clear L1 L2;

%%
qq_bin = 3;
ii1 = S1(:,:,qq_bin,1);
ii2 = S2(:,:,qq_bin,1);

T1 = S1(:,:,:,1)./sqrt(L1);
T2 = S2(:,:,:,1)./sqrt(L2);
a1 = squeeze(T1(:,:,3)); %this is very close to 1

%%
fpos = [550, 0, 560, 1100];

figure(1);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(T1(:,:,3), [0,2]); colorbar; title('t');
subplot(4,1,2); imagesc(QUV1(:,:,qq_bin,1), [-1,1]); colorbar; title('q');
subplot(4,1,3); imagesc(QUV1(:,:,qq_bin,2), [-1,1]); colorbar; title('u');
subplot(4,1,4); imagesc(QUV1(:,:,qq_bin,3), [-1,1]); colorbar; title('v');
sgtitle('QUV1') 

function data = read_stokes(folder_processed, fildID, fkey)
% fkey = 'sv1';
ffile_h5 = fullfile(folder_processed, [fildID, sprintf('.%s.h5',fkey)]);

% h5disp(ffile_h5)
h5_start_vec = [1,1,1,1,1]; %[frame, bin, stokes, z_ind, x_ind] 
data = h5read(ffile_h5, sprintf('/%s',fkey), h5_start_vec,[1, Inf, Inf, Inf, Inf]);
% size(data)
data = permute(data, [4,5,2,3,1]); 
data = squeeze(data);
% size(data)
end