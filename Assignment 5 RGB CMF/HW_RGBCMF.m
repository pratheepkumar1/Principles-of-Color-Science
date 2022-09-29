lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
rgb_spd_dataset = readtable("HW_DisplaySPD_Data.xlsx",Sheet="DisplaySPD");
cie_dataset = readtable("HW_DisplaySPD_Data.xlsx",Sheet="CIE 1931");

% Create LMS relative sensitivity matrix
lms_matrix = transpose(lms_dataset{:,[2:4]});

rgb_spd = rgb_spd_dataset{:,[2:4]};

%Change in wavlength (d lamda)
lambda = mean(diff(lms_dataset{:,1}));

%LMS values of the display’s three primaries
LMS_RGB = transpose((lms_matrix * rgb_spd)*lambda)

%Normalizing LMS by peak value
M = max(LMS_RGB, [], 'all')
LMS_RGB = LMS_RGB/M



%color matching functions of the display’s primaries
CMF_RGB = inv(LMS_RGB)*lms_matrix;

%Normalizing CMF by peak value
C = max(CMF_RGB, [], 'all');
CMF_RGB = (CMF_RGB/C)
CMF_RGB_transpose = transpose(CMF_RGB)



colors=['r','g','b'];
labels = ['Red','Green','Blue']

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



%Converting from CMFs to CIE standard colorimeter observer



% Fixing residual error in calculating
CIE_source = cie_dataset{:,[2:4]};

% x_lamda = transpose(CIE_source(:,1))*pinv(CMF_RGB)
y_lamda = transpose(CIE_source(:,2))*pinv(CMF_RGB)
% z_lamda = transpose(CIE_source(:,3))*pinv(CMF_RGB)
% m_function = [x_lamda;y_lamda;z_lamda]




n_denominator = y_lamda(1)*(CMF_RGB_transpose(:,1)) + y_lamda(2)*(CMF_RGB_transpose(:,2)) + y_lamda(3)*(CMF_RGB_transpose(:,3))

CIE_source(:,2)
n_denominator
% n_lamda = cie_xyz(:,2)*transpose(n_denominator)
% 
% 
% transpose(m_function*(CMF_RGB*n_lamda))

% cie_dataset(:,"y_bar")
% CMF_RGB(:,1)


%Least Squares(
% transpose(cie_xyz(:,2))*pinv(CMF_RGB)






% y_bar = cie_dataset{:,[2]}
