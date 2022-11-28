%% Reading the Dataset

lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
rgb_spd_dataset = readtable("HW_DisplaySPD_Data.xlsx",Sheet="DisplaySPD");
cie_dataset = readtable("HW_DisplaySPD_Data.xlsx",Sheet="CIE 1931");

%% Question 1 - 3

% Create LMS relative sensitivity matrix
lms_matrix = transpose(lms_dataset{:,[2:4]});

% Creating Spectral Radiance matrix
rgb_spd = transpose(rgb_spd_dataset{:,[2:4]});

%Change in Wavlength
wavelength = lms_dataset{:,1};
d_lambda = mean(diff(wavelength));

%QUESTION 1: LMS values of the displayQs three primaries
LMS_RGB = calc_LMS_source(lms_matrix,transpose(rgb_spd),d_lambda);
LMS_RGB_norm = custom_normalization(LMS_RGB);
% M = max(LMS_RGB, [], 'all');
% LMS_RGB = LMS_RGB/M;


%QUESTION 2a: Color matching functions of the displayQs primaries
CMF_RGB = LMS_RGB\lms_matrix;
CMF_RGB_norm = custom_normalization(CMF_RGB);


%Normalizing CMF by peak value
% C = max(CMF_RGB, [], 'all');
% CMF_RGB = (CMF_RGB/C);


%Plotting the graph
CMF_RGB_norm_transpose = transpose(CMF_RGB_norm);
figure;
colors=['r','g','b'];
labels = ['Red','Green','Blue'];
for i = 1:3
    hold on
    plot(lms_dataset{:,1},CMF_RGB_norm_transpose(:,i),"LineWidth",2,"Color",colors(i))
end
grid on
xlabel('Wavelength (in nm)')
ylabel('Tristimulus Values')
xlim([380 730])
title('Color Matching Functions of Primaries')
hold off


%% Question 4 & 5

%Create Source matrix for
CIE_source = cie_dataset{:,2:4};

%Normalize r,g,b to unit area invidually and create a matrix
r_bar_norm = custom_normalization(CMF_RGB(1,:));
g_bar_norm = custom_normalization(CMF_RGB(2,:));
b_bar_norm = custom_normalization(CMF_RGB(3,:));
rgb_bar_norm = [r_bar_norm;g_bar_norm;b_bar_norm];

%Assigning y-bar value from the dataset to V_lamda
V_lamda = transpose(CIE_source(:,2));


% Least squares used to calculate the scalars of each color-matching function
M_function = transpose(CIE_source)*pinv(rgb_bar_norm);

% Calculating n_lamda to fix the residual error
n_denominator = M_function(2,1)*(rgb_bar_norm(1,:)) + M_function(2,2)*(rgb_bar_norm(2,:)) + M_function(2,3)*(rgb_bar_norm(3,:));
n_lamda = V_lamda/n_denominator;
n_lamda_diagonal = diag(ones(36,1)*n_lamda);

% Calculate xyz-bar for the displays 
xyz_bar = calc_CIE(M_function,rgb_bar_norm,n_lamda_diagonal);
xyz_bar_norm = (custom_normalization(xyz_bar));


%Plotting the graphs
CIE_source_norm = custom_normalization(CIE_source);
xyz_bar_norm_transpose = transpose(xyz_bar_norm);

figure;
colors=['r','g','b'];
for i = 1:3
    hold on
    plot(lms_dataset{:,1},xyz_bar_norm_transpose(:,i),"LineWidth",1,"Color",colors(i));
    plot(lms_dataset{:,1},CIE_source_norm(:,i),"LineWidth",1,"Color",colors(i),"LineStyle","--");
end
grid on
xlabel('Wavelength (in nm)');
ylabel('Tristimulus Values');
xlim([lms_dataset{1,1} lms_dataset{end,1}]);
title('Color Matching Functions of Primaries');
legend('LCD display x','CIE1931 std observer x','LCD display y','CIE1931 std observer y','LCD display z','CIE1931 std observer z');
hold off


% Function to normalize a matrix with its peak value (not using default normalize function)
function n = custom_normalization(x)
    max_value = max(x, [], 'all');
    n = x/max_value;
end

% Function to calculate LMS of light source
function T = calc_LMS_source(t,s,d_lambda)
    T = t*s*d_lambda;
end

%Function to calculate the CIE 1931 standard colorimetric observer
function CIE = calc_CIE(M,rgb_bar,n_lamda)
    CIE = M*(rgb_bar*n_lamda);
end


