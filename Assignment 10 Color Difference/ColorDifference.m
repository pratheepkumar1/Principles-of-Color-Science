% Reading the Dataset
warning('off','all');
inkjetColorChecker_dataset = readtable("InkjetColorChecker.xlsx");
source_dataset = readtable("Illuminant Data.xlsx");
xyz_std_obs_two_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="2-degree");
patch_dataset = readtable("MacbethColorChecker.xlsx");


%Getting wavelength information for adjusting the dataset
disp("Note: For this assignment, enter min value = 380, max value = 780, range = 5");
prompt = "Enter min wavelength for interpolation = ";
min_value = input(prompt);
prompt = "Enter max wavelength for interpolation = ";
max_value = input(prompt);
prompt = "Enter wavelength dataset interval range = ";
d_lambda = input(prompt);
wavelength_info = min_value:d_lambda:max_value;


%Interpolating and extrapolating the dataset (as needed)
xyz_std_obs_two_deg = interp1(xyz_std_obs_two_deg_dataset{:,1},xyz_std_obs_two_deg_dataset{:,2:end},wavelength_info,"linear","extrap");
source_A = interp1(source_dataset{:,1},source_dataset.A,wavelength_info,"linear","extrap");
source_D65 = interp1(source_dataset{:,1},source_dataset.D65,wavelength_info,"linear","extrap");
patch_macbeth_ref = interp1(patch_dataset{2:end,1},patch_dataset{2:end,2:end},wavelength_info,"linear","extrap");
inkjetColorChecker = interp1(inkjetColorChecker_dataset{2:end,1},inkjetColorChecker_dataset{2:end,2:end},wavelength_info,"linear","extrap");


% White point value of Source A and Source D65 for 2 degree observer
wp_A_two_deg_source = calcTristimulusSource(xyz_std_obs_two_deg,source_A,d_lambda);
wp_D65_two_deg_source = calcTristimulusSource(xyz_std_obs_two_deg,source_D65,d_lambda);


%% Question 1

% Tristimulus value for all patches in InkjetColorChecker in A and D65 Source for 2 degree observer
tristimulus_XYZ_A_two_deg_inkjet = calcTristimulus(xyz_std_obs_two_deg,source_A,inkjetColorChecker,d_lambda);
tristimulus_XYZ_D65_two_deg_inkjet = calcTristimulus(xyz_std_obs_two_deg,source_D65,inkjetColorChecker,d_lambda);


%Calculate L_star, a_star, b_star, C_star for all patches in InkjetColorChecker in A and D65 Source 
[L_A_inkjet, a_A_inkjet, b_A_inkjet, C_A_inkjet] = calcXYZtoCIELAB(tristimulus_XYZ_A_two_deg_inkjet,wp_A_two_deg_source);
[L_D65_inkjet, a_D65_inkjet, b_D65_inkjet, C_D65_inkjet] = calcXYZtoCIELAB(tristimulus_XYZ_D65_two_deg_inkjet,wp_D65_two_deg_source);

CIELAB_A_inkjet = [L_A_inkjet; a_A_inkjet; b_A_inkjet; C_A_inkjet];
CIELAB_D65_inkjet = [L_D65_inkjet; a_D65_inkjet; b_D65_inkjet; C_D65_inkjet];

table_CIELAB_A_inkjet = createTableCIELAB(CIELAB_A_inkjet)
table_CIELAB_D65_inkjet = createTableCIELAB(CIELAB_D65_inkjet)


%Visualizing to view the difference
% hold on
% scatter(a_A',b_A','red');
% scatter(a_A_inkjet',b_A_inkjet','blue');
% axis equal
% grid on
% hold off


%% Question 2 (Prerequisite Calculations)

% Tristimulus value for all patches in InkjetColorChecker in A and D65 Source for 2 degree observer
tristimulus_XYZ_A_two_deg_macbeth_ref = calcTristimulus(xyz_std_obs_two_deg,source_A,patch_macbeth_ref,d_lambda);
tristimulus_XYZ_D65_two_deg_macbeth_ref = calcTristimulus(xyz_std_obs_two_deg,source_D65,patch_macbeth_ref,d_lambda);

%Calculate L_star, a_star, b_star, C_star for all patches in InkjetColorChecker in A and D65 Source 
[L_A_macbeth_ref, a_A_macbeth_ref, b_A_macbeth_ref, C_A_macbeth_ref] = calcXYZtoCIELAB(tristimulus_XYZ_A_two_deg_macbeth_ref,wp_A_two_deg_source);
[L_D65_macbeth_ref, a_D65_macbeth_ref, b_D65_macbeth_ref, C_D65_macbeth_ref] = calcXYZtoCIELAB(tristimulus_XYZ_D65_two_deg_macbeth_ref,wp_D65_two_deg_source);

CIELAB_A_macbeth_ref = [L_A_macbeth_ref; a_A_macbeth_ref; b_A_macbeth_ref; C_A_macbeth_ref];
CIELAB_D65_macbeth_ref = [L_D65_macbeth_ref; a_D65_macbeth_ref; b_D65_macbeth_ref; C_D65_macbeth_ref];

table_CIELAB_A_macbeth_ref = createTableCIELAB(CIELAB_A_macbeth_ref)
table_CIELAB_D65_macbeth_ref = createTableCIELAB(CIELAB_D65_macbeth_ref)


%% Question 2


% Delta calculation for A
deltaL_A = deltaL(L_A_inkjet,L_A_macbeth_ref(:,13:24));

deltaa_A = deltaa(a_A_inkjet,a_A_macbeth_ref(:,13:24));

deltab_A = deltab(b_A_inkjet,b_A_macbeth_ref(:,13:24));

deltaC_A = deltaC(C_A_inkjet,C_A_macbeth_ref(:,13:24));

deltah_A = deltah(CIELAB_A_inkjet,CIELAB_A_macbeth_ref(:,13:24));

deltaH_A = deltaH(C_A_inkjet,C_A_macbeth_ref(:,13:24),deltah_A);

deltaEab_A = deltaEab(deltaL_A,deltaC_A,deltaH_A);

deltaE00_A = deltaE00(CIELAB_A_inkjet,CIELAB_A_macbeth_ref(:,13:24));

color_diff_A = [deltaL_A; deltaa_A; deltab_A; deltaC_A; deltaH_A; deltaEab_A; deltaE00_A];

table_color_diff_A = createTableColorDiff(color_diff_A)

figure(1)
hold on
x_value = linspace(1,12,12);
scatter(x_value,deltaEab_A','filled','red');
scatter(x_value,deltaE00_A',"filled",'green');
hold off


% Delta calculation for D65

deltaL_D65 = deltaL(L_D65_inkjet,L_D65_macbeth_ref(:,13:24));

deltaa_D65 = deltaa(a_D65_inkjet,a_D65_macbeth_ref(:,13:24));

deltab_D65 = deltab(b_D65_inkjet,b_D65_macbeth_ref(:,13:24));

deltaC_D65 = deltaC(C_D65_inkjet,C_D65_macbeth_ref(:,13:24));

deltah_D65 = deltah(CIELAB_D65_inkjet,CIELAB_D65_macbeth_ref(:,13:24));

deltaH_D65 = deltaH(C_D65_inkjet,C_D65_macbeth_ref(:,13:24),deltah_D65);

deltaEab_D65 = deltaEab(deltaL_D65,deltaC_D65,deltaH_D65);

deltaE00_D65 = deltaE00(CIELAB_D65_inkjet,CIELAB_D65_macbeth_ref(:,13:24));

color_diff_D65 = [deltaL_D65; deltaa_D65; deltab_D65; deltaC_D65; deltaH_D65; deltaEab_D65; deltaE00_D65];

table_color_diff_D65 = createTableColorDiff(color_diff_D65)


figure(2)
hold on
x_value = linspace(1,12,12);
scatter(x_value,deltaEab_D65','filled','red');
scatter(x_value,deltaE00_D65',"filled",'green');
hold off


figure(3)
hold on
bar(deltaEab_D65);
bar(deltaE00_D65);
hold off



%% Functions

%Function to create table and assign CIELAB values
function a2t = createTableCIELAB(values)
    a2t = array2table(values');
    a2t.Properties.VariableNames(1:4) = {'L','a','b','C'};
end

function a2t = createTableColorDiff(values)
    a2t = array2table(values');
    a2t.Properties.VariableNames(1:7) = {'ΔL','Δa','Δb','ΔC','ΔH','ΔEab','ΔE00'};
end

%Function to calculate Delta L
function DL = deltaL(L_bat,L_std)
    DL = L_bat - L_std;
end


%Function to calculate Delta a
function Da = deltaa(a_bat,a_std)
    Da = a_bat - a_std;
end


%Function to calculate Delta b
function Db = deltab(b_bat,b_std)
    Db = b_bat - b_std;
end

%Function to calculate Delta C*
function DC = deltaC(C_bat,C_std)
    DC = C_bat - C_std;
end


%Function to calculate Delta h
function Dh = deltah(bat,std)
    Dh = hueAngle(bat(2,:),bat(3,:)) - hueAngle(std(2,:),std(3,:));
end


%Function to calculate hue angle
function h = hueAngle(a,b)
    h = atan2Deg(b,a);
end


%Function to calculate Delta H
function DH = deltaH(C_bat,C_std,deltah)
    DH = 2*(C_bat.*C_std).^(1/2).*sinDeg(deltah./2);
end


%Function to calculate Delta Eab
function e = deltaEab(deltaL,deltaa,deltab)
    e = sqrt(deltaL.^2 + deltaa.^2 + deltab.^2);
end


function out = atan2Deg(inY,inX)
    out = atan2(inY,inX).*180./pi;
    out = out+(out<0).*360;
end


function out = sinDeg(in)
    out = sin(in.*pi./180);
end






