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

CIELAB_A_inkjet = [L_A_inkjet; a_A_inkjet; b_A_inkjet; C_A_inkjet];
CIELAB_D65_inkjet = [L_D65_inkjet; a_D65_inkjet; b_D65_inkjet; C_D65_inkjet];

table_A = createTableCIELAB(CIELAB_A_inkjet);
table_D65 = createTableCIELAB(CIELAB_D65_inkjet);


%Visualizing to view the difference
% hold on
% scatter(a_A',b_A','red');
% scatter(a_A_inkjet',b_A_inkjet','blue');
% axis equal
% grid on
% hold off


%% Question 2

deltaL_A = deltaL(L_A_inkjet,L_A(:,13:24))

deltaa_A = deltaa(a_A_inkjet,a_A(:,13:24))

deltab_A = deltab(b_A_inkjet,b_A(:,13:24))

deltaC_A = deltaC(C_A_inkjet,C_A(:,13:24))

deltaH_A = deltaH(CIELAB_A_inkjet,CIELAB_A(:,13:24))

deltaEab_A = deltaEab(deltaL_A,deltaa_A,deltab_A) 

deltaE00_A = deltaE00(CIELAB_A_inkjet,CIELAB_A(:,13:24))

% 
% ((a_bat*b_std) - (a_std*b_bat)) / ...
%     ((0.5*((C_bat*C_std) + (a_bat*a_std) + (b_bat*b_std)))^1/2)


%% Functions

%Function to create table and assign CIELAB values
function a2t = createTableCIELAB(values)
    a2t = array2table(values');
    a2t.Properties.VariableNames(1:4) = {'L','a','b','C'};
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


%Function to calculate Delta H
function Dh = deltaH(bat,std)
    Dh = hueAngle(bat(2,:),bat(3,:)) - hueAngle(std(2,:),std(3,:));
end



%Function to calculate hue angle
function h = hueAngle(a,b)
    h = atan2Deg(b,a);
    if(h<0)
        h = (180/pi)*h;
    else
        h = ((180/pi)*h)+360;
    end
end


%Function to calculate Delta Eab
function e = deltaEab(deltaL,deltaa,deltab)
    e = (deltaL.^2 + deltaa.^2 + deltab.^2).^1/2;
end


function out = atan2Deg(inY,inX);
    out = atan2(inY,inX).*180./pi;
    out = out+(out<0).*360;
end








