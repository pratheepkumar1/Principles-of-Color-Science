% Reading the Dataset
warning('off','all');
all_patches = readtable("color_science_lab_measurements.csv",'Range','B1:AP17');
blue_patch37_dataset = readtable("Blue37_Patch.csv",'Range','B1:AP11');
source_dataset = readtable("Illuminant Data.xlsx");
xyz_std_obs_two_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="2-degree");


%Formatting the dataset (as needed)
wavelength_info = 360:10:750;
bluish_grey_patches = all_patches{2:6,2:end};
grey_patches = all_patches{7:11,2:end};
orange_patches = all_patches{12:16,2:end};
blue_patch37 = blue_patch37_dataset{2:11,2:end};

xyz_std_obs_two_deg = interp1(xyz_std_obs_two_deg_dataset{:,1},xyz_std_obs_two_deg_dataset{:,2:end},wavelength_info,"linear","extrap");
source_D65 = interp1(source_dataset{:,1},source_dataset.D65,wavelength_info,"linear","extrap");

%% Question 1

% figure(1);
% hold on
% line(wavelength_info,bluish_grey_patches);
% legend('Patch 31','Patch 32','Patch 33','Patch 34','Patch 35')
% xlabel('Wavelength');
% ylabel('Reflectance Factor');
% title('Bluish Grey Patches')
% hold off
% 
% 
% figure(2);
% hold on
% line(wavelength_info,grey_patches);
% legend('Patch 19','Patch 20','Patch 21','Patch 22','Patch 23')
% xlabel('Wavelength');
% ylabel('Reflectance Factor');
% title('Grey Patches')
% hold off
% 
% 
% figure(3);
% hold on
% line(wavelength_info,orange_patches);
% legend('Patch 12','Patch 16','Patch 28','Patch 38','Patch 39')
% xlabel('Wavelength');
% ylabel('Reflectance Factor');
% title('Orange/Yellow Patches')
% hold off


%% Question 2

figure(4);
hold on
line(wavelength_info,blue_patch37);
xlabel('Wavelength');
ylabel('Reflectance Factor');
legend('Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Trial 6','Trial 7','Trial 8','Trial 9','Trial 10')
title('Patch 37 - Blue')
hold off


%% Question 3


% White point value of Source A and Source D65 for 2 degree observer
wp_D65_two_deg_source = calcTristimulusSource(xyz_std_obs_two_deg,source_D65,10);


% Tristimulus value for all patches in InkjetColorChecker in A and D65 Source for 2 degree observer
tristimulus_XYZ_blue_patch37 = calcTristimulus(xyz_std_obs_two_deg,source_D65,blue_patch37',10);


%Calculate L_star, a_star, b_star, C_star for all patches in InkjetColorChecker in A and D65 Source 
[L_blue_patch37, a_blue_patch37, b_blue_patch37, C_blue_patch37] = calcXYZtoCIELAB(tristimulus_XYZ_blue_patch37,wp_D65_two_deg_source);
CIELAB_blue_patch37 = [L_blue_patch37; a_blue_patch37; b_blue_patch37; C_blue_patch37];
table_CIELAB_blue_patch37 = createTableCIELAB(CIELAB_blue_patch37)



%% Functions

%MCDM Calculation
MCDM_blue_patch37 = calcMCDM(CIELAB_blue_patch37(1,:),CIELAB_blue_patch37(2,:),CIELAB_blue_patch37(3,:)) 

function mcdm = calcMCDM(L_star,a_star,b_star)
    L_star_mean = mean(L_star);
    a_star_mean = mean(a_star);
    b_star_mean = mean(b_star);
    n = numel(L_star);
    mcdm = (sum(((L_star-L_star_mean).^2+(a_star-a_star_mean).^2+(b_star-b_star_mean).^2).^(1/2)))./n;
end


%Function to create table and assign CIELAB values
function a2t = createTableCIELAB(values)
    a2t = array2table(values');
    a2t.Properties.VariableNames(1:4) = {'L','a','b','C'};
end


%Function to calculate CIELAB from XYZ
function [L_star,a_star,b_star,C_star] = calcXYZtoCIELAB(XYZ,whitepoint)
    XYZ_Prime = calcXYZPrime(XYZ,whitepoint);
    L_star = (116 * calConstants(XYZ_Prime(2,:))) - 16;
    a_star= 500*(calConstants(XYZ_Prime(1,:)) - calConstants(XYZ_Prime(2,:)));
    b_star = 200*(calConstants(XYZ_Prime(2,:)) - calConstants(XYZ_Prime(3,:)));
    C_star = ((a_star).^2 + (b_star).^2).^(1/2);
end


% Calculate X Prime, Y Prime and Z Prime
function XYZ_Prime = calcXYZPrime(XYZ,whitepoint)
    XYZ_Prime = XYZ./whitepoint';
end


%Function to calculate constants in a and b
function k = calConstants(x)
    if (x > (24/116)^3)
        k = (x).^(1/3);
    else
        k = ((841/108).*x) + (16/116);
    end
end

%Function to calculate Tristimulus values for materials
function t = calcTristimulus(xyz_value,source,material,d_lambda)   

    %Calculate normalizing constant
    k = (100/(source * xyz_value(:,2) * d_lambda))/100;

    s_lambda = diag(source);

    %Calculating tristimulus
    t = k.*((s_lambda*xyz_value)'*material)*d_lambda;

%   t = custom_normalization(t);
end


%Function to calculate Tristimulus values for sources (Whitepoint)
function ts = calcTristimulusSource(xyz_value,source,d_lambda)   

    %Calculate normalizing constant
    k = (100/(source * xyz_value(:,2) * d_lambda))/100;

    %Calculating tristimulus
    ts = k.*((source*xyz_value))*d_lambda;

%   t = custom_normalization(t);
end

