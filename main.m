darkCurrentPhotons = 5;
conversionFactor = 500;
% darkCurrentPhotons = 0;
% conversionFactor = 500;

dirtemp = ['conversionFactor_',num2str(conversionFactor),'_darkCurrentPhotons_',num2str(darkCurrentPhotons)];
mkdir(dirtemp);


% Define simulation volume [x, y, z]
vol_size = [128, 128, 5];

% Set ground truth parameters
start_point = [5, 5, 1];
end_point = [123, 123, 5];
fwhm_axial = 200; 
fwhm_radial = 2; 

% Generate ground truth
ground_truth = createGaussianBeam3D(vol_size, start_point, end_point, fwhm_axial, fwhm_radial);
ground_truth = double(ground_truth);

saveTiffStackSingleFile(ground_truth,[dirtemp,'/ground_truth.tiff']);

noise_3D = addPMTShotNoise3D(ground_truth, darkCurrentPhotons, conversionFactor);

saveTiffStackSingleFile(normxin(noise_3D),[dirtemp,'/noise_3D_norm.tiff']);
saveTiffStackSingleFile((noise_3D),[dirtemp,'/noise_3D.tiff']);


[Nyy,Nxx,Nzz] = size(noise_3D);
pixel_numel = Nyy*Nxx;

sample_matrix = ones(Nzz,Nzz) - eye([Nzz,Nzz]);
% sample_matrix = eye([Nzz,Nzz]);

ground_truth_line = reshape(ground_truth,[pixel_numel,Nzz]);

sample_result_line = zeros([Nzz,pixel_numel]);

for ii = 1:1:pixel_numel
    pixel_temp = ground_truth_line(ii,:);
    pixel_temp = pixel_temp';
    sample_result_temp = sample_matrix*pixel_temp;
    sample_result_line(:,ii) = sample_result_temp;
end

sample_result = reshape(sample_result_line',[Nyy,Nxx,Nzz]);
sample_result_noisy = addPMTShotNoise3D(sample_result, darkCurrentPhotons, conversionFactor);
saveTiffStackSingleFile(sample_result_noisy,[dirtemp,'/sample_result_noisy.tiff']);
sample_result_noisy_line = reshape(sample_result_noisy,[pixel_numel,Nzz]);
sample_result_noisy_line = sample_result_noisy_line';

resolved_result_line = sample_result_line*0;
resolved_result_line0 = resolved_result_line;
resolved_result_line2 = resolved_result_line;

%% Optimization settings
options = optimoptions('fmincon',...
    'Display', 'iter',...      
    'Algorithm', 'interior-point',... 
    'MaxIterations', 1000);     % Maximum number of iterations

lb = zeros(Nzz, 1);  % x >= 0
ub = [];           

parfor ii = 1:pixel_numel
    pixel_temp = sample_result_noisy_line(:,ii);
    pixel_resolved0 = lsqnonneg(sample_matrix, pixel_temp);
    resolved_result_line0(:,ii) = pixel_resolved0;

    negLogLikelihood = @(pixel_resolved) sum( -pixel_temp'*log(max(sample_matrix*pixel_resolved, eps)) + sum(sample_matrix*pixel_resolved) );
%     pixel_resolved0 = pixel_resolved0*0;
    [pixel_resolved, fval] = fmincon(negLogLikelihood, pixel_resolved0, [], [], [], [], lb, ub, [], options);
    
    resolved_result_line(:,ii) = pixel_resolved;

    resolved_result_line2(:,ii) = sample_matrix \ pixel_temp;

    
end

resoved_result = reshape(resolved_result_line',[Nyy,Nxx,Nzz]);
resoved_result0 = reshape(resolved_result_line0',[Nyy,Nxx,Nzz]);
resoved_result2 = reshape(resolved_result_line2',[Nyy,Nxx,Nzz]);



resoved_result_norm = normxin(resoved_result);
saveTiffStackSingleFile(resoved_result_norm,[dirtemp,'\resoved_result_norm.tiff']);
saveTiffStackSingleFile(resoved_result,[dirtemp,'\resoved_result.tiff']);


resoved_result_norm0 = normxin(resoved_result0);
saveTiffStackSingleFile(resoved_result_norm0,[dirtemp,'\resoved_result_norm0.tiff']);
saveTiffStackSingleFile(resoved_result0,[dirtemp,'\resoved_result0.tiff']);

resoved_result_norm2 = normxin(resoved_result2);
saveTiffStackSingleFile(resoved_result_norm2,[dirtemp,'\resoved_result_norm2.tiff']);
saveTiffStackSingleFile(resoved_result2,[dirtemp,'\resoved_result2.tiff']);
saveTiffStackSingleFile(resoved_result2/maxxin(resoved_result2),[dirtemp,'\resoved_result_norm2_2.tiff']);



psnr_noise = psnr(ground_truth,normxin(noise_3D))
psnr_resoved_result_norm = psnr(ground_truth,resoved_result_norm)
psnr_resoved_result_norm0 = psnr(ground_truth,resoved_result_norm0)
psnr_resoved_result_norm2 = psnr(ground_truth,resoved_result_norm2)


ssim_noise = ssim(ground_truth,normxin(noise_3D))
ssim_resoved_result_norm = ssim(ground_truth,resoved_result_norm)
ssim_resoved_result_norm0 = ssim(ground_truth,resoved_result_norm0)
ssim_resoved_result_norm2 = ssim(ground_truth,resoved_result_norm2)



save([dirtemp,'\alldata.mat']);
