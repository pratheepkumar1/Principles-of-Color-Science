%% Reading the Dataset
warning('off','all');
SGref_dataset = readtable("ccsg.xlsx");
SGref = SGref_dataset{:,5:end};
SGwl = (380:10:730)';

%Compute the principal components of 36 wavelengths using 
[coeff,score,latent,tsquared,explained,mu] = pca(SGref);



figure;
hold on
plot(SGwl,coeff(:,1:5),'LineWidth',1.5);
xlabel('Wavelength');
ylabel('PCA Coefficient');
xlim([380,730]);
legend("e1","e2","e3","e4","e5");
grid on
hold off


%Calculating cumulative sum of the variance

cumsum_var = cumsum(latent);
cumsum_var_percent = (cumsum_var/cumsum_var(end,:))*100;

figure;
hold on
coeff_count = 1:1:36;
plot(coeff_count,cumsum_var_percent,'LineWidth',1.5,'Color','red');
xlabel('Number of Components');
ylabel('Cumulative Variance (%)');
xticks(1:1:36)
grid on
hold off

%% Question 2
SG_reconst = score*coeff'+mu;

%% Question 3

% Reconstructed CCSG data with 1st principal component
SG_reconst_first = score(:,1)*coeff(:,1)'+mu;

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

% Reconstructed CCSG data with 2nd principal component
SG_reconst_second = score(:,1:2)*coeff(:,1:2)'+mu;

SGXYZ_reconst_second = ((cmf2' * diag(D65) * SG_reconst_second') ./ (cmf2(:,2)' * D65) )';
SGLab_reconst_second = xyz2lab(SGXYZ_reconst_second);

De00_second = deltaE00(SGLab_reconst_second',SGLab_ref');
De00_second_mean = mean(De00_second);
De00_second_max = max(De00_second);



%% Bonus

for i = 1:36
    SG_reconst_f = score(:,1:i)*coeff(:,1:i)'+mu; 
    SGXYZ_reconst_f = ((cmf2' * diag(D65) * SG_reconst_f') ./ (cmf2(:,2)' * D65) )';
    SGLab_reconst_f = xyz2lab(SGXYZ_reconst_f);
    De00_f = deltaE00(SGLab_reconst_f',SGLab_ref');
    De00_f_mean = mean(De00_f);
    De00_f_max = max(De00_f);
    if(De00_f_max <= 0.5)
        disp("Number of coefficients to minimize De00 to 0.5 is " + i)
        break;
    end
end

    



%% Functions


