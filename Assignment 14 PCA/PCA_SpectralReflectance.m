%% Reading the Dataset
warning('off','all');
SGref_dataset = readtable("ccsg.xlsx");
SG_ref = SGref_dataset{:,5:end};
SGwl = (380:10:730)';

% Compute the principal components of 36 wavelengths using 
[coeff,score,latent,tsquared,explained,mu] = pca(SG_ref);

% Displaying the principal components
coeff;

% Plot of the first 5 principal components versus wavelength
figure;
hold on
plot(SGwl,coeff(:,1:5),'LineWidth',1.5);
xlabel('Wavelength');
ylabel('PCA Coefficient');
xlim([380,730]);
legend("e1","e2","e3","e4","e5");
grid on
hold off


%Calculating & plotting cumulative percentage of the variance
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

% Reconstruct the 140 spectra of CCSG patches
SG_reconst = score*coeff'+mu;

%Confirm if reconstruct match
isequal(fi(SG_ref),fi(SG_reconst))

%% Question 3

% Compute CCSG reference LAB values for De00 calcualtion
[num,txt] = xlsread('all_1nm_data.xlsx');
cmf2 = interp1(num(:,1),num(:,6:8),SGwl);
D65 = interp1(num(:,1),num(:,3),SGwl);
SGXYZ_ref = ((cmf2' * diag(D65) * SG_ref') ./ (cmf2(:,2)' * D65) )';
SGLab_ref = xyz2lab(SGXYZ_ref);

% Reconstructed CCSG data with 1st PCA coefficient
SG_reconst_first = score(:,1)*coeff(:,1)'+mu;

% Compute CCSG 1st PCA coefficient data's LAB values
SGXYZ_reconst_first = ((cmf2' * diag(D65) * SG_reconst_first') ./ (cmf2(:,2)' * D65) )';
SGLab_reconst_first = xyz2lab(SGXYZ_reconst_first);

% DE00 color difference between the original and reconstructed
De00_first = deltaE00(SGLab_reconst_first',SGLab_ref');

% De00 mean
De00_first_mean = mean(De00_first);

% De00 Max Value
De00_first_max = max(De00_first);


%% Question 4

% Reconstructed CCSG data with first two PCA coefficient
SG_reconst_second = score(:,1:2)*coeff(:,1:2)'+mu;

% Compute CCSG 2nd PCA coefficient data's LAB values
SGXYZ_reconst_second = ((cmf2' * diag(D65) * SG_reconst_second') ./ (cmf2(:,2)' * D65) )';
SGLab_reconst_second = xyz2lab(SGXYZ_reconst_second);

% DE00 color difference between the original and reconstructed
De00_second = deltaE00(SGLab_reconst_second',SGLab_ref');

% De00 mean
De00_second_mean = mean(De00_second);

% De00 Max Value
De00_second_max = max(De00_second);



%% Bonus

% Determine number of pca coefficients to reduce the max DE00 below 0.5.
for i = 1:36
    SG_reconst_f = score(:,1:i)*coeff(:,1:i)'+mu; 
    SGXYZ_reconst_f = ((cmf2' * diag(D65) * SG_reconst_f') ./ (cmf2(:,2)' * D65) )';
    SGLab_reconst_f = xyz2lab(SGXYZ_reconst_f);
    De00_f = deltaE00(SGLab_reconst_f',SGLab_ref');
    De00_f_max = max(De00_f);
    if(De00_f_max <= 0.5)
        disp("Number of coefficients to minimize De00 to 0.5 is " + i)
        disp("The De00 value at " + i + " is " + De00_f_max);
        break;
    end
end

    



%% Functions

function   De00=deltaE00(Lab1, Lab2)

%CIELAB Chroma
C1 = sqrt(Lab1(2,:).^2+Lab1(3,:).^2);
C2 = sqrt(Lab2(2,:).^2+Lab2(3,:).^2);

%Lab Prime
mC = (C1+C2)./2;
G=0.5*(1-sqrt((mC.^7)./((mC.^7)+(25.^7))));
LabP1 = [Lab1(1,:) ; Lab1(2,:).*(1+G) ; Lab1(3,:)];
LabP2 = [Lab2(1,:) ; Lab2(2,:).*(1+G) ; Lab2(3,:)];

%Chroma
CP1 = sqrt(LabP1(2,:).^2+LabP1(3,:).^2);
CP2 = sqrt(LabP2(2,:).^2+LabP2(3,:).^2);

%Hue Angle
hP1t = atan2Deg(LabP1(3,:),LabP1(2,:));
hP2t = atan2Deg(LabP2(3,:),LabP2(2,:));

%Add in 360 to the smaller hue angle if absolute value of difference is > 180
hP1 = hP1t + ((hP1t<hP2t)&(abs(hP1t-hP2t)>180)).*360;
hP2 = hP2t + ((hP1t>hP2t)&(abs(hP1t-hP2t)>180)).*360;

%Delta Values
DLP = LabP1(1,:) - LabP2(1,:);
DCP = CP1 - CP2;
DhP = hP1 - hP2;
DHP = 2*(CP1.*CP2).^(1/2).*sinDeg(DhP./2);

%Arithmetic mean of LCh' values
mLP = (LabP1(1,:)+LabP2(1,:))./2;
mCP = (CP1+CP2)./2;
mhP = (hP1+hP2)./2;

%Weighting Functions
SL = 1+(0.015.*(mLP-50).^2)./sqrt(20+(mLP-50).^2);
SC = 1+0.045.*mCP;
T = 1-0.17.*cosDeg(mhP-30)+0.24.*cosDeg(2.*mhP)+0.32.*cosDeg(3.*mhP+6)-0.2.*cosDeg(4.*mhP-63);
SH = 1+0.015.*mCP.*T;

%Rotation function
RC = 2.*sqrt((mCP.^7)./((mCP.^7)+25.^7));
DTheta = 30.*exp(-((mhP-275)./25).^2);
RT = -sinDeg(2.*DTheta).*RC;

%Parametric factors
kL = 1;
kC = 1;
kH = 1;

De00 = ((DLP./kL./SL).^2+(DCP./kC./SC).^2+(DHP./kH./SH).^2+(RT.*(DCP./kC./SC).*(DHP./kH./SH))).^(1/2);

end

%--------------
function out = sinDeg(in);
out = sin(in.*pi./180);
end

%--------------
function out = cosDeg(in);
out = cos(in.*pi./180);
end

%--------------
function out = atan2Deg(inY,inX);
out = atan2(inY,inX).*180./pi;
out = out+(out<0).*360;

end

