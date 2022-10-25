% Reading the Dataset
warning('off','all');
inkjetColorChecker_dataset = readtable("InkjetColorChecker.xlsx");
source_dataset = readtable("Illuminant Data.xlsx");


%Getting the inkjetcolorchecker values after interpolation
inkjetColorChecker = InterpolateData(inkjetColorChecker_dataset{2:end,1},inkjetColorChecker_dataset{2:end,2:end},intp_wavelength_info);


%Calculate change in Wavlength
wavelength = inkjetColorChecker_dataset{2:end,1};
d_lambda = mean(diff(wavelength));


%% Question 1

% Tristimulus value for all patches in InkjetColorChecker in A and D65 Source for 2 degree observer
tristimulus_XYZ_A_two_deg_inkjet = calcTristimulus(xyz_std_obs_two_deg,source_A,inkjetColorChecker,d_lambda);
tristimulus_XYZ_D65_two_deg_inkjet = calcTristimulus(xyz_std_obs_two_deg,source_D65,inkjetColorChecker,d_lambda);


%Calculate L_star, a_star, b_star, C_star for all patches in InkjetColorChecker in A and D65 Source 
[L_A_inkjet, a_A_inkjet, b_A_inkjet, C_A_inkjet] = calcXYZtoCIELAB(tristimulus_XYZ_A_two_deg_inkjet,wp_A_two_deg_source);
[L_D65_inkjet, a_D65_inkjet, b_D65_inkjet, C_D65_inkjet] = calcXYZtoCIELAB(tristimulus_XYZ_D65_two_deg_inkjet,wp_D65_two_deg_source);


table_A = createTableCIELAB([L_A_inkjet; a_A_inkjet; b_A_inkjet; C_A_inkjet])
table_D65 = createTableCIELAB([L_D65_inkjet; a_D65_inkjet; b_D65_inkjet; C_D65_inkjet])


%Visualizing to view the difference
% hold on
% scatter(a_A',b_A','red');
% scatter(a_A_inkjet',b_A_inkjet','blue');
% axis equal
% grid on
% hold off


%% Question 2



%% Functions

%Function to create table and assign CIELAB values
function a2t = createTableCIELAB(values)
    a2t = array2table(values');
    a2t.Properties.VariableNames(1:4) = {'L','a','b','C'};
end