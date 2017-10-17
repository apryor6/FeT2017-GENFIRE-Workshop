noise_level = 5e-6;

noisy_pj = importdata('projections_FePt.mat');
noisy_pj = noisy_pj + randn(size(noisy_pj)) * noise_level;
noisy_pj(noisy_pj<0) = 0;

save('projections_FePt_noisy','noisy_pj')


misaligned_pj = importdata('projections_FePt.mat');
shifts = round(randn(2, size(misaligned_pj,3)) * 1);
for i = 1:size(misaligned_pj, 3)
    misaligned_pj(:, :, i) = circshift(misaligned_pj(:, :, i), [shifts(1,i), shifts(2,i)]);
end

save('projections_FePt_misaligned','misaligned_pj')
