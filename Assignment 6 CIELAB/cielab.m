% Reading the Dataset
warning('off','all');
xyz_std_obs_two_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="2-degree");
xyz_std_obs_ten_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="10-degree");
source_dataset = readtable("Illuminant Data.xlsx");
patch_dataset = readtable("MacbethColorChecker.xlsx");
cc_spectral_locus_dataset = readtable("TwoDegChromaticity.xlsx");
% cc_spectral_locus_cc = interpolateData(cc_spectral_locus_dataset{:,1},cc_spectral_locus_dataset{:,2:end},intp_wavelength_info);


% Create a interpolated data for the datasets
intp_wavelength_info = struct('min',380,'max',780,'range',5);

xyz_std_obs_two_deg = interpolateData(xyz_std_obs_two_deg_dataset{:,1},xyz_std_obs_two_deg_dataset{:,2:end},intp_wavelength_info);
xyz_std_obs_ten_deg = interpolateData(xyz_std_obs_ten_deg_dataset{:,1},xyz_std_obs_ten_deg_dataset{:,2:end},intp_wavelength_info);
patches = patch_dataset{2:end,2:25};

%Calculate change in Wavlength
wavelength = patch_dataset{2:end,1};
d_lambda = mean(diff(wavelength));

%% Question 1

%Get spectral radiance of Source D50
sources_D50 = interpolateData(source_dataset{:,1},source_dataset.D50,intp_wavelength_info);

% Tristimulus value for all patches in D50 Source for 2 and 10 degree observer
tristimulus_XYZ_D50_two_deg = calcTristimulus(xyz_std_obs_two_deg,sources_D50,patches,d_lambda);
tristimulus_XYZ_D50_ten_deg = calcTristimulus(xyz_std_obs_ten_deg,sources_D50,patches,d_lambda);

% chromacity_coordinates for all patches in D50 Source for 2 and 10 degree observer
chromacity_coord_D50_two_deg = calChromacityCoordinates(tristimulus_XYZ_D50_two_deg);
chromacity_coord_D50_ten_deg = calChromacityCoordinates(tristimulus_XYZ_D50_ten_deg);

table_cc_D50 = createTablexyz(chromacity_coord_D50_two_deg);

%% Question 2

%Get spectral radiance of Source A and Source D65
sources_A = interpolateData(source_dataset{:,1},source_dataset.A,intp_wavelength_info);
sources_D65 = interpolateData(source_dataset{:,1},source_dataset.D65,intp_wavelength_info);


% Tristimulus value for all patches in Source A and Source D65 for 2 degree observer
tristimulus_XYZ_A_two_deg = calcTristimulus(xyz_std_obs_two_deg,sources_A,patches,d_lambda);
tristimulus_XYZ_D65_two_deg = calcTristimulus(xyz_std_obs_two_deg,sources_D65,patches,d_lambda);

% chromacity_coordinates for all patches in Source A and Source D65 for 2 degree observer
chromacity_coord_A_two_deg = calChromacityCoordinates(tristimulus_XYZ_A_two_deg);
chromacity_coord_D65_two_deg = calChromacityCoordinates(tristimulus_XYZ_D65_two_deg);

table_cc_A = createTablexyz(chromacity_coord_A_two_deg);
table_cc_D65 = createTablexyz(chromacity_coord_D65_two_deg);


%% Question 3


% Working
figure(1)
hold on
scatter(table_cc_D50{:,1},table_cc_D50{:,2},20,'filled','Color','#8EB5E0');
scatter(table_cc_A{:,1},table_cc_A{:,2},20,'filled','Color','#CEDAE8');
scatter(table_cc_D65{:,1},table_cc_D65{:,2},20,'filled','Color','#4988CF');
axis equal
xlim([0 0.8]);
ylim([0 0.8]);
legend('A','D50','D65');
hold off


% figure(2)
% plotChromaticity
% hold on
% scatter3(table_cc_A{:,1},table_cc_A{:,2},table_cc_A{:,3},24,'filled','Color','#8EB5E0');
% hold on
% scatter3(table_cc_D50{:,1},table_cc_D50{:,2},table_cc_D50{:,3},24,'filled','Color','#8EB5E0');
% hold on
% scatter3(table_cc_D65{:,1},table_cc_D65{:,2},table_cc_D65{:,3},24,'filled','Color','#8EB5E0');
% xlabel('x');
% ylabel('y');
% zlabel('z');
% grid on
% axis equal
% xlim([0 0.8]);
% ylim([0 0.8]);
% legend('A','D50','D65');
% hold off

% wp = calChromacityCoordinates(tristimulus_XYZ_D65_two_deg_source');
% x_whitepoint = wp(1,:);
% y_whitepoint = wp(2,:);
% plotChromaticity
% hold on
% scatter(x_whitepoint,y_whitepoint,36,'filled','black');
% scatter(table_cc_D65{:,1},table_cc_D65{:,2},24,'black');
% axis equal
% hold off

% Working
figure(3)
hold on
plot(cc_spectral_locus_dataset{:,2},cc_spectral_locus_dataset{:,3},'b-x');
scatter(table_cc_D50{:,1},table_cc_D50{:,2},20,'filled','Color','#8EB5E0');
scatter(table_cc_A{:,1},table_cc_A{:,2},20,'filled','Color','#CEDAE8');
scatter(table_cc_D65{:,1},table_cc_D65{:,2},20,'filled','Color','#4988CF');
axis equal
xlim([0 1]);
ylim([0 1]);
legend('A','D50','D65');
hold off


%% Question 6 (Bonus)

%Get spectral radiance of Source A and Source D65
sources_A = interpolateData(source_dataset{:,1},source_dataset.A,intp_wavelength_info);
sources_D65 = interpolateData(source_dataset{:,1},source_dataset.D65,intp_wavelength_info);

% Tristimulus value of Source A and Source D65 for 2 degree observer
wp_A_two_deg_source = calcTristimulusSource(xyz_std_obs_two_deg,sources_A,d_lambda);
wp_D65_two_deg_source = calcTristimulusSource(xyz_std_obs_two_deg,sources_D65,d_lambda);

% Chromaticity Coordinates of Source A and Source D65 for 2 degree observer
cc_wp_A_two_deg = calChromacityCoordinates(wp_A_two_deg_source');
cc_wp_D65_two_deg = calChromacityCoordinates(wp_D65_two_deg_source');

hold on
plot(cc_spectral_locus_dataset{:,2},cc_spectral_locus_dataset{:,3},'-');
scatter(cc_wp_A_two_deg(1,:),cc_wp_A_two_deg(2,:),36,'filled','black');
scatter(table_cc_A{:,1},table_cc_A{:,2},24,'black');
scatter(cc_wp_D65_two_deg(1,:),cc_wp_D65_two_deg(2,:),36,'filled','green');
scatter(table_cc_D65{:,1},table_cc_D65{:,2},24,'green');
axis equal
hold off


%% Question 7


%Calculate L_star, a_star, b_star, C_star for all patches in A and D65 source 
[L_A, a_A, b_A, C_A] = calcXYZtoCIELAB(tristimulus_XYZ_A_two_deg,wp_A_two_deg_source);
[L_D65, a_D65, b_D65, C_D65] = calcXYZtoCIELAB(tristimulus_XYZ_D65_two_deg,wp_D65_two_deg_source);

%Plot a against b under source A
figure(5)
hold on
scatter(a_A',b_A');
scatter(a_D65',b_D65');
axis equal
hold off

figure(6)
hold on
scatter(C_A',L_A');
scatter(C_D65',L_D65');
axis equal
hold off






%% Functions


%Function to interpolate the values to a consistent wavelength
function i = interpolateData(wavelengths,values,intp_wavelength_data)
    %Create the interpolant by passing wavelength and corresponding values to griddedInterpolant.
    GI = griddedInterpolant(wavelengths,values);

    % Create a vector of query points with 5nm spacing for 380 to 780 wavelength.
    wl = intp_wavelength_data.min:intp_wavelength_data.range:intp_wavelength_data.max;

    % Evaluate the interpolant at the each wavelength for each value set 
    i = GI(wl);
end


%Function to calculate Tristimulus values for materials
function t = calcTristimulus(xyz_value,source,material,d_lambda)   

    %Calculate normalizing constant
    k = 100/(source * xyz_value(:,2) * d_lambda);

    s_lambda = diag(source);

    %Calculating tristimulus
    t = k.*((s_lambda*xyz_value)'*material)*d_lambda;

%   t = custom_normalization(t);
end


%Function to calculate Chromacity Coordinate
function cc = calChromacityCoordinates(tristimulus_XYZ)
    sum_XYZ = sum(tristimulus_XYZ,1);
    cc = tristimulus_XYZ./sum_XYZ;
end


%Function to calculate Tristimulus values for sources (Whitepoint)
function ts = calcTristimulusSource(xyz_value,source,d_lambda)   

    %Calculate normalizing constant
    k = 100/(source * xyz_value(:,2) * d_lambda);

    %Calculating tristimulus
    ts = k.*((source*xyz_value))*d_lambda;

%   t = custom_normalization(t);
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


%Function to create table and assign x,y,z
function a2t = createTablexyz(values)
    a2t = array2table(values');
    a2t.Properties.VariableNames(1:3) = {'x','y','z'};
end

% Function to normalize a matrix with its peak value (not using default normalize function)
function n = custom_normalization(x)
    max_value = max(x, [], 'all');
    n = x/max_value;
end




