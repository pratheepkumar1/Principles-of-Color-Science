%% Question 1 - 3

% Reading the Data Set
lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
rgb_spd_dataset = readtable("HW_DisplaySPD_Data.xlsx",Sheet="DisplaySPD");
cie_dataset = readtable("HW_DisplaySPD_Data.xlsx",Sheet="CIE 1931");

% Create LMS relative sensitivity matrix
lms_matrix = transpose(lms_dataset{:,[2:4]});

% Creating Spectral Radiance matrix
rgb_spd = transpose(rgb_spd_dataset{:,[2:4]});

%Change in wavlength (d lamda)
lambda = mean(diff(lms_dataset{:,1}));

%QUESTION 1: LMS values of the displayQs three primaries
LMS_RGB = (lms_matrix * transpose(rgb_spd))*lambda;

%Normalizing LMS by peak value
M = max(LMS_RGB, [], 'all');
LMS_RGB = LMS_RGB/M;


%QUESTION 2a: Color matching functions of the displayQs primaries
CMF_RGB = inv(LMS_RGB)*lms_matrix;


%Normalizing CMF by peak value
C = max(CMF_RGB, [], 'all');
CMF_RGB = (CMF_RGB/C);
CMF_RGB_transpose = transpose(CMF_RGB);




%Plotting the graph
figure(1)
colors=['r','g','b'];
labels = ['Red','Green','Blue'];
for i = 1:3
    hold on
    plot(lms_dataset{:,1},CMF_RGB_transpose(:,i),"LineWidth",2,"Color",colors(i))
end
grid on
hold off
xlabel('Wavelength (in nm)')
ylabel('Tristimulus Values')
xlim([380 730])
title('Color Matching Functions of Primaries')

%% Question 4 & 5

CIE_source = cie_dataset{:,[2:4]};
CIE_source_norm = normalize(CIE_source,'norm',1)

%Normalize r,g,b to unit area invidually and create a matrix
r_bar = normalize(CMF_RGB(1,:),'norm',1)
g_bar = normalize(CMF_RGB(2,:),'norm',1)
b_bar = normalize(CMF_RGB(3,:),'norm',1)
rgb_bar = [r_bar;g_bar;b_bar]

%Assigning y-bar value from the dataset to V_lamda
V_lamda = transpose(CIE_source(:,2));


% Least squares used to calculate the scalars of each color-matching function
M_function = transpose(CIE_source)*pinv(rgb_bar)


% Calculating n_lamda to fix the residual error
n_denominator = M_function(2,1)*(rgb_bar(1,:)) + M_function(2,2)*(rgb_bar(2,:)) + M_function(2,1)*(rgb_bar(3,:));
n_lamda = V_lamda/n_denominator;


% m_function
xyz_bar = transpose(CIE_source).*((rgb_bar)*diag(n_lamda))

xyz_bar_norm = normalize(xyz_bar,'norm',1)
xyz_bar_norm_transpose = transpose(xyz_bar)

%Plotting the graphs
figure(2)
colors=['r','g','b'];
labels = ['Red','Green','Blue'];
for i = 1:3
    hold on
    plot(lms_dataset{:,1},xyz_bar_norm_transpose(:,i),"LineWidth",2,"Color",colors(i))
     plot(lms_dataset{:,1},CIE_source_norm(:,i),"LineWidth",2,"Color",colors(i),"LineStyle","--")
end

grid on
hold off
xlabel('Wavelength (in nm)')
ylabel('Tristimulus Values')
xlim([380 730])
title('Color Matching Functions of Primaries')






% y_bar = cie_dataset{:,[2]}
