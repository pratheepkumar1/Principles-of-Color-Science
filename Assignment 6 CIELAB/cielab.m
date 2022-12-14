% Reading the Dataset
warning('off','all');
xyz_std_obs_two_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="2-degree");
xyz_std_obs_ten_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="10-degree");
source_dataset = readtable("Illuminant Data.xlsx");
patch_dataset = readtable("MacbethColorChecker.xlsx");
cc_spectral_locus_dataset = readtable("TwoDegChromaticity.xlsx");
% cc_spectral_locus_cc = interpolateData(cc_spectral_locus_dataset{:,1},cc_spectral_locus_dataset{:,2:end},intp_wavelength_info);


% Create a interpolated data for the datasets
prompt = "Enter min wavelength for interpolation = ";
min_value = input(prompt);
prompt = "Enter max wavelength for interpolation = ";
max_value = input(prompt);
prompt = "Enter wavelength dataset interval range = ";
range_value = input(prompt);
intp_wavelength_info = struct('min',min_value,'max',max_value,'range',range_value);

xyz_std_obs_two_deg = InterpolateData(xyz_std_obs_two_deg_dataset{:,1},xyz_std_obs_two_deg_dataset{:,2:end},intp_wavelength_info);
xyz_std_obs_ten_deg = InterpolateData(xyz_std_obs_ten_deg_dataset{:,1},xyz_std_obs_ten_deg_dataset{:,2:end},intp_wavelength_info);
patches = InterpolateData(patch_dataset{2:end,1},patch_dataset{2:end,2:end},intp_wavelength_info);
% patches = patch_dataset{2:end,2:25};

%Calculate change in Wavlength
% wavelength = patch_dataset{2:end,1};
d_lambda = range_value;

%% Question 1

%Get spectral radiance of Source D50
source_D50 = InterpolateData(source_dataset{:,1},source_dataset.D50,intp_wavelength_info);

% Tristimulus value for all patches in D50 Source for 2 and 10 degree observer
tristimulus_XYZ_D50_two_deg = calcTristimulus(xyz_std_obs_two_deg,source_D50,patches,d_lambda);
tristimulus_XYZ_D50_ten_deg = calcTristimulus(xyz_std_obs_ten_deg,source_D50,patches,d_lambda);

% chromacity_coordinates for all patches in D50 Source for 2 and 10 degree observer
cc_D50_two_deg = calChromacityCoordinates(tristimulus_XYZ_D50_two_deg);
cc_D50_ten_deg = calChromacityCoordinates(tristimulus_XYZ_D50_ten_deg);

table_cc_D50 = createTablexyz(cc_D50_two_deg);

%% Question 2

%Get spectral radiance of Source A and Source D65
source_A = InterpolateData(source_dataset{:,1},source_dataset.A,intp_wavelength_info);
source_D65 = InterpolateData(source_dataset{:,1},source_dataset.D65,intp_wavelength_info);


% Tristimulus value for all patches in Source A and Source D65 for 2 degree observer
tristimulus_XYZ_A_two_deg = calcTristimulus(xyz_std_obs_two_deg,source_A,patches,d_lambda);
tristimulus_XYZ_D65_two_deg = calcTristimulus(xyz_std_obs_two_deg,source_D65,patches,d_lambda);

% chromacity_coordinates for all patches in Source A and Source D65 for 2 degree observer
cc_A_two_deg = calChromacityCoordinates(tristimulus_XYZ_A_two_deg);
cc_D65_two_deg = calChromacityCoordinates(tristimulus_XYZ_D65_two_deg);

table_cc_A = createTablexyz(cc_A_two_deg);
table_cc_D65 = createTablexyz(cc_D65_two_deg);


%% Question 3


% % Working
% hold on
% scatter(table_cc_D50{:,1},table_cc_D50{:,2},32,'MarkerFaceColor','#0F4786');
% scatter(table_cc_A{:,1},table_cc_A{:,2},32,'MarkerEdgeColor','#930796','MarkerFaceColor','#930796');
% scatter(table_cc_D65{:,1},table_cc_D65{:,2},32,'MarkerFaceColor','#089E32');
% plot(cc_spectral_locus_dataset{:,2},cc_spectral_locus_dataset{:,3},'-','LineWidth',2,'Color','black');
% axis equal
% xlim([0 0.9]);
% ylim([0 0.9]);
% xlabel('x');
% ylabel('y');
% title('Chromaticity Coordinates of 24 Patches for D50, A, D65');
% legend('D50','A','D65','Spectrum Locus');
% dcm_obj = datacursormode(h);
% set(dcm_obj,'UpdateFcn',{@custom_label})
% hold off



%% Question 6 (Bonus)

%Get spectral radiance of Source A and Source D65
source_A = InterpolateData(source_dataset{:,1},source_dataset.A,intp_wavelength_info);
source_D65 = InterpolateData(source_dataset{:,1},source_dataset.D65,intp_wavelength_info);

% Tristimulus value of Source A and Source D65 for 2 degree observer
wp_A_two_deg_source = calcTristimulusSource(xyz_std_obs_two_deg,source_A,d_lambda);
wp_D65_two_deg_source = calcTristimulusSource(xyz_std_obs_two_deg,source_D65,d_lambda);

% Chromaticity Coordinates of Source A and Source D65 for 2 degree observer
cc_wp_A_two_deg = calChromacityCoordinates(wp_A_two_deg_source');
cc_wp_D65_two_deg = calChromacityCoordinates(wp_D65_two_deg_source');

% cc_spectral_locus_dataset(1:81,2:end)

% Calculate the purity of each patch
for i = 1:24
    a = norm(cc_wp_A_two_deg(1:2,:)' - cc_D50_two_deg(1:2,i)');
%     b = norm(cc_wp_A_two_deg(1:2,:)' - cc_spectral_locus_dataset(1:81,2:end))
end



%% Question 7

%Calculate L_star, a_star, b_star, C_star for all patches in A and D65 source 
[L_A, a_A, b_A, C_A] = calcXYZtoCIELAB(tristimulus_XYZ_A_two_deg,wp_A_two_deg_source);
[L_D65, a_D65, b_D65, C_D65] = calcXYZtoCIELAB(tristimulus_XYZ_D65_two_deg,wp_D65_two_deg_source);


CIELAB_A = [L_A; a_A; b_A; C_A];
CIELAB_D65 = [L_D65; a_D65; b_D65; C_D65];


%Plot a against b under source A
% h = figure(3);
% hold on
% scatter(a_A',b_A');
% xlabel('a^{*}');
% ylabel('b^{*}');
% title('CIELAB a^{*} vs b^{*} for 24 patches in A illuminant');
% xlim([-100 100]);
% ylim([-100 100]);
% axis equal
% grid on
% dcm_obj = datacursormode(h);
% set(dcm_obj,'UpdateFcn',{@custom_label})
% hold off
% 
% u = figure(4);
% hold on
% scatter(C_A',L_A');
% xlabel('Chroma (C^{*}_{ab})');
% ylabel('Lightness (L^{*}_{ab})');
% title('CIELAB L^{*}_{ab} vs C^{*}_{ab} for 24 patches in A illuminant');
% axis equal
% grid on
% dcm_obj = datacursormode(u);
% set(dcm_obj,'UpdateFcn',{@custom_label})
% hold off
% 
% 
% r = figure(5);
% hold on
% scatter(a_D65',b_D65');
% xlabel('a^{*}');
% ylabel('b^{*}');
% title('CIELAB a^{*} vs b^{*} for 24 patches in D65 illuminant');
% xlim([-100 100]);
% ylim([-100 100]);
% axis equal
% grid on
% dcm_obj = datacursormode(r);
% set(dcm_obj,'UpdateFcn',{@custom_label})
% hold off
% 
% b = figure(6);
% hold on
% scatter(C_D65',L_D65');
% xlabel('Chroma (C^{*}_{ab})');
% ylabel('Lightness (L^{*}_{ab})');
% title('CIELAB L^{*}_{ab} vs C^{*}_{ab} for 24 patches in D65 illuminant');
% axis equal
% grid on
% dcm_obj = datacursormode(b);
% set(dcm_obj,'UpdateFcn',{@custom_label})
% hold off


%% Functions

%Function to calculate Chromacity Coordinate
function cc = calChromacityCoordinates(tristimulus_XYZ)
    sum_XYZ = sum(tristimulus_XYZ,1);
    cc = tristimulus_XYZ./sum_XYZ;
end



%Function to create table and assign x,y,z
function a2t = createTablexyz(values)
    a2t = array2table(values');
    a2t.Properties.VariableNames(1:3) = {'x','y','z'};
end


% Function to customize text of data tips in the plot
function txt = custom_label(~,event_obj)
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Patch: ',num2str(I)]};
end

% Function to normalize a matrix with its peak value (not using default normalize function)
function n = custom_normalization(x)
    max_value = max(x, [], 'all');
    n = x/max_value;
end




