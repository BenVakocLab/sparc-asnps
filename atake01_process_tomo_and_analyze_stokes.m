% exampls script for codeCleanUp
% saved data structure only supports the legacy software for now.
% last update 20150715 by Stephanie Nam (snam@mit.edu)
addpath('C:\Users\vlab_eye\Documents\local_repo\vlab-ofdi')

% path of the folder containing all acquired data and  
directory = 'C:\Users\vlab_eye\Documents\local_repo\ps-closer-look\data\[p.spinal_cord][11-01-2021_16-09-36]';

b_processed = true;

if ~b_processed 
    slices = 380+ (1:10); %(1:1184);%(300:300);
    zoomFactor = 4;%4;
    monitorBscan = 0;
    monitorAline = 0;
    flipUpDown = 1;
    tomoOpt = v2struct(slices, zoomFactor, monitorBscan, monitorAline, flipUpDown);
    ofdiTOMOGRAPHYclean(directory,tomoOpt);

    tomoOpt_gpu = tomoOpt;
    tomoOpt_gpu.gpuFlag = 0;
    [complexTomoXfilename, complexTomoYfilename] = ofdiTOMOGRAPHYclean_gpu(directory,tomoOpt_gpu); 
%%%>> will need to change REF_lowhigh values when using gpu code.
else
%%
% slices = 1:1; %570:570;%
% meta5tomo = [0 3 1024 1024 10];
% monitorBscan = 1;
% REF_lowhigh = [-10 90];
% 
% structOpt = v2struct(slices,meta5tomo,monitorBscan, REF_lowhigh);
% ofdiSTRUCTUREclean(directory,structOpt);

complexTomoXfilename = '[p.spinal_cord][11-01-2021_16-09-36][0.3.1024.1024.10].scatxX.mgh';
complexTomoYfilename = '[p.spinal_cord][11-01-2021_16-09-36][0.3.1024.1024.10].scatxY.mgh';

tomo_detX = readMgh(fullfile(directory, complexTomoXfilename), 1);
tomo_detY = readMgh(fullfile(directory, complexTomoYfilename), 1);

end

%% Create Stokes vector

tomch1_input1 = tomo_detX(:, 1:2:end);
tomch1_input2 = tomo_detX(:, 2:2:end);
tomch2_input1 = tomo_detY(:, 1:2:end);
tomch2_input2 = tomo_detY(:, 2:2:end);

i1 = abs(tomch1_input1).^2 + abs(tomch2_input1).^2;
q1 = abs(tomch1_input1).^2 - abs(tomch2_input1).^2;
u1 = 2*real(tomch1_input1.*conj(tomch2_input1));
v1 = -2*imag(tomch1_input1.*conj(tomch2_input1));

i2 = abs(tomch1_input2).^2 + abs(tomch2_input2).^2;
q2 = abs(tomch1_input2).^2 - abs(tomch2_input2).^2;
u2 = 2*real(tomch1_input2.*conj(tomch2_input2));
v2 = -2*imag(tomch1_input2.*conj(tomch2_input2));

sv1_to_analyze = cat(3, i1,q1,u1,v1);
sv2_to_analyze = cat(3, i2,q2,u2,v2);

i0 = (sv1_to_analyze(:,:,1)+sv2_to_analyze(:,:,1))/2;

unaveraged_quv1 = sv1_to_analyze(:,:,2:4)./sv1_to_analyze(:,:,1);
unaveraged_quv2 = sv2_to_analyze(:,:,2:4)./sv2_to_analyze(:,:,1);

unaveraged_i1 = sqrt(sum(unaveraged_quv1.^2, 4));
unaveraged_i2 = sqrt(sum(unaveraged_quv2.^2, 4));

fpos = [550, 0, 560, 1100];
figure(3);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(10*log10(sv1_to_analyze(:,:,1)), [-20, 50] ); colorbar; title('i1 [dB]');
subplot(4,1,2); imagesc(unaveraged_quv1(:,:,1), [-1,1]); colorbar; title('q');
subplot(4,1,3); imagesc(unaveraged_quv1(:,:,2), [-1,1]); colorbar; title('u');
subplot(4,1,4); imagesc(unaveraged_quv1(:,:,3), [-1,1]); colorbar; title('v');
sgtitle('unaveraged quv1') 
figure(4);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(10*log10(sv2_to_analyze(:,:,1)), [-20, 50] ); colorbar; title('i2 [dB]');
subplot(4,1,2); imagesc(unaveraged_quv2(:,:,1), [-1,1]); colorbar; title('q');
subplot(4,1,3); imagesc(unaveraged_quv2(:,:,2), [-1,1]); colorbar; title('u');
subplot(4,1,4); imagesc(unaveraged_quv2(:,:,3), [-1,1]); colorbar; title('v');
sgtitle('unaveraged quv2') 


%% Analyze Stokes average
kk = hanning(3)*hanning(11)'; %row col
kk = kk./sum(kk(:))

for ww = 1:size(sv1_to_analyze,3)
    for qq = 1:size(sv1_to_analyze,4)
        sv1_to_analyze(:,:,ww,qq) = imfilter(sv1_to_analyze(:,:,ww,qq), kk);
        sv2_to_analyze(:,:,ww,qq) = imfilter(sv2_to_analyze(:,:,ww,qq), kk);
    end
end

averaged_i1 = sqrt(sum(sv1_to_analyze(:,:,2:4).^2, 3));
averaged_i2 = sqrt(sum(sv2_to_analyze(:,:,2:4).^2, 3));

averaged_quv1 = sv1_to_analyze(:,:,2:4)./averaged_i1;
averaged_quv2 = sv2_to_analyze(:,:,2:4)./averaged_i2;

dop = (averaged_i1+averaged_i2)./(sv1_to_analyze(:,:,1)+sv2_to_analyze(:,:,1));
dop1 = averaged_i1./sv1_to_analyze(:,:,1);
dop2 = averaged_i2./sv2_to_analyze(:,:,1);


figure(5);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(dop1, [0,1]); colorbar; title('dop1');
subplot(4,1,2); imagesc(averaged_quv1(:,:,1), [-1,1]); colorbar; title('q');
subplot(4,1,3); imagesc(averaged_quv1(:,:,2), [-1,1]); colorbar; title('u');
subplot(4,1,4); imagesc(averaged_quv1(:,:,3), [-1,1]); colorbar; title('v');
sgtitle('averaged quv1') 

figure(6);
set(gcf, 'Position', fpos);
subplot(4,1,1); imagesc(dop2, [0,1]); colorbar; title('dop2');
subplot(4,1,2); imagesc(averaged_quv2(:,:,1), [-1,1]); colorbar; title('q');
subplot(4,1,3); imagesc(averaged_quv2(:,:,2), [-1,1]); colorbar; title('u');
subplot(4,1,4); imagesc(averaged_quv2(:,:,3), [-1,1]); colorbar; title('v');
sgtitle('averaged quv2') 


%% Analyze Stokes trajectory along depth
select_zpos = 350+ (1:80); 
select_bin = 1 ;


figure(100);
imagesc(10*log10(i0)); colorbar; title('ii [dB]'); 
colorbar; 
hold on;

for alinepos = [201]% [260] %[600:100:800]%1000:100:1600
    select_aline = alinepos;
figure(100);
plot(ones(size(select_zpos))*select_aline, select_zpos, 'k'); 
plot(select_aline, select_zpos(1),'ko');

select_nsv1 = squeeze(averaged_quv1(select_zpos, alinepos, :));
select_nsv2 = squeeze(averaged_quv2(select_zpos, alinepos, :));

nfig = 101;
hfig = figure(nfig); 
[x,y,z] = sphere(16);
hsph = surf(x,y,z);  set(hsph, 'FaceAlpha',0.2, 'LineStyle', 'none'); 
colormap([1 1 0]); hold on;
plot3(select_nsv1(1,1),select_nsv1(1,2),select_nsv1(1,3), 'or', 'MarkerSize',12); 
plot3(select_nsv1(1:end,1),select_nsv1(1:end,2),select_nsv1(1:end,3),'.-r'); 

plot3(select_nsv2(1,1),select_nsv2(1,2),select_nsv2(1,3),'ob', 'MarkerSize',12); 
plot3(select_nsv2(1:end,1),select_nsv2(1:end,2),select_nsv2(1:end,3),'.-b'); 


end 
figure(100); hold off;
figure(nfig); hold off;

