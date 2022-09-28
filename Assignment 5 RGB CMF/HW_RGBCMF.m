lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
rgb_spd_dataset = readtable("HW_DisplaySPD_Data.xlsx",Sheet="DisplaySPD");

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
CMF_RGB = transpose(inv(LMS_RGB)*lms_matrix);

%Normalizing CMF by peak value
C = max(CMF_RGB, [], 'all')
CMF_RGB = CMF_RGB/C


colors=['r','g','b'];
labels = ['Red','Green','Blue']

for i = 1:3
    hold on
    plot(lms_dataset{:,1},CMF_RGB(:,i),"LineWidth",2,"Color",colors(i))
end
grid on
hold off
xlabel('Wavelength (in nm)')
ylabel('Tristimulus Values')
title('Color Matching Functions of Primaries')

