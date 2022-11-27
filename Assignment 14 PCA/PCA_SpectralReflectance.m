%% Reading the Dataset
warning('off','all');
SGref_dataset = readtable("ccsg.xlsx");
SGref = SGref_dataset{:,5:end};

%Compute the principal components of 36 
[coeff,score,latent,tsquared,explained,mu] = pca(SGref);
coeff;

SGwl = (380:10:730)';

figure;
hold on
plot(SGwl,coeff(1:5,:));
xlabel('Wavelength');
ylabel('PCA Coefficient');
xlim([380,730]);
grid on
hold off

coeff_count = 1:1:36;

figure;
hold on
bar(coeff_count,latent,1);
xlabel('Principal Component Number');
ylabel('cumulative percentage of variance');
xticks(1:1:36)
grid on
hold off

%% Question 2
SG_reconst = score*coeff'+mu;

%% Question 3

% Reconstructed CCSG data with 1st principal component
SG_reconst_first = score.*coeff(:,1)'+mu;

[num,txt] = xlsread('all_1nm_data.xlsx');
cmf2 = interp1(num(:,1),num(:,6:8),SGwl);
D65 = interp1(num(:,1),num(:,3),SGwl);
SGXYZ_ref = ((cmf2' * diag(D65) * SGref') ./ (cmf2(:,2)' * D65) )';
SGLab_ref = xyz2lab(SGXYZ_ref);

SGXYZ_reconst_first = ((cmf2' * diag(D65) * SG_reconst_first') ./ (cmf2(:,2)' * D65) )';
SGLab_reconst_first = xyz2lab(SGXYZ_reconst_first);

De00_first = deltaE00(SGLab_reconst_first',SGLab_ref');
De00_first_mean = mean(De00_first);
De00_first_max = max(De00_first);


%% Question 4

% Reconstructed CCSG data with 1st principal component
SG_reconst_second = score.*coeff(:,1:2)+mu;

SGXYZ_reconst_second = ((cmf2' * diag(D65) * SG_reconst_second') ./ (cmf2(:,2)' * D65) )';
SGLab_reconst_second = xyz2lab(SGXYZ_reconst_second);

De00_second = deltaE00(SGLab_reconst_second',SGLab_ref');
De00_second_mean = mean(De00_second);
De00_second_max = max(De00_second);


%% Bonus

% coefficient_number = 1;
% 
% for i=1:36
%     SG_reconst = score.*coeff(:,)'+mu;
% end
% 
% if (De00_max <= 0.5)
%     


