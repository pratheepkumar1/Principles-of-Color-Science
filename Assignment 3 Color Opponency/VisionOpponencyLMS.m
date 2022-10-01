color_checker_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="ColorChecker");
lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
light_sources_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="Sources");

% Create LMS relative sensitivity matrix
lms_matrix = transpose(lms_dataset{:,[2:4]});


[light_sources_numRows,light_sources_numCols] = size(light_sources_dataset); 
normalize_s_factor = 100*ones(light_sources_numRows,1); %normalise the matrix values to 0-1
incandescent_source = light_sources_dataset{:,2}./ normalize_s_factor;
daylight_source= light_sources_dataset{:,3} ./ normalize_s_factor;
n_lamda = mean(diff(light_sources_dataset{:,1}))



% create a diagonal matrix to calculate illuminant's spectral power distribution = S𝜆
s_incandescent_matrix = diag(incandescent_source);
s_daylight_matrix = diag(daylight_source);

% Objectspectral reflectance factor = R𝜆 
[color_checker_numRows,color_checker_numCols] = size(color_checker_dataset);
r_matrix = color_checker_dataset{:,[2:color_checker_numCols]};

% LMS values for all 24 ColorChecker samples
LMS_incandescent_stimuli = LMS_stimuli(lms_matrix, s_incandescent_matrix , r_matrix,n_lamda);
LMS_daylight_stimuli = LMS_stimuli(lms_matrix, s_daylight_matrix, r_matrix, n_lamda);


LMS_incandescent_stimuli
LMS_daylight_stimuli

%-----------------------------------------
% Chromatic Adaptation Transformations (CAT)

% Equi-energy light source is L1
% Incandescent is L2

% LMS_daylight_lightsource
LMS_daylight_lightsource = transpose(lms_matrix * s_daylight_matrix)
LMS_incandescent_lightsource = transpose(lms_matrix * s_incandescent_matrix)




% Getting signals of a stimulus that appears white
white_col_number = find(string(color_checker_dataset.Properties.VariableNames) == "White");
LMS_white_incandescent = transpose(LMS_incandescent_stimuli(white_col_number,:));
LMS_white_daylight = transpose(LMS_daylight_stimuli(white_col_number,:));

% D_factor (complete chromatic adaptation = 1)
degree_of_adaptation = 1;



% m_von_kries_incadescent_factor = diag((degree_of_adaptation*LMS_white_incandescent)+...
%     ((1-degree_of_adaptation)*LMS_white_daylight))./(LMS_white_daylight);
% 
% m_von_kries_daylight_factor = diag((degree_of_adaptation*LMS_white_daylight)+...
%     ((1-degree_of_adaptation)*LMS_white_incandescent))./(LMS_white_incandescent);

LMS_daylight_von_kries = transpose(m_von_kries_daylight_factor * transpose(LMS_incandescent_stimuli));
LMS_incandescent_von_kries = transpose(m_von_kries_incadescent_factor * transpose(LMS_daylight_stimuli));

LMS_incandescent_von_kries
LMS_daylight_von_kries

%--------------------------------------------
%Opponency Signals
m_opponency = [.64 .39 -0.01; 1.12 -1.5 0.34; 0.35 0.15 -0.53];
gamma = 2.4;
opponency_incandescent = transpose(m_opponency * transpose((LMS_incandescent_von_kries).^(1/gamma)));
opponency_daylight = transpose(m_opponency * transpose((LMS_daylight_von_kries).^(1/gamma)));
opponency_incandescent;
opponency_incandescent(19:24,:);
opponency_incandescent(16,:);

opponency_incandescent(:,1)
opponency_daylight

disp("Opponency for neutral patches (19-24)")
opponency_incandescent(19:24,:)
disp("Opponency for yellow (16)")
opponency_incandescent(16,:)

% plot(opponency_incandescent(19:24,:))

% hold on
% L_daylight = plot(LMS_daylight(:,[1]));
% % M_daylight = plot(LMS_daylight(:,[2]));
% % S_daylight = plot(LMS_daylight(:,[3]));
% L_daylight_von_kries = plot(LMS_daylight_von_kries(:,[1]));
% % M_daylight_von_kries = plot(LMS_daylight_von_kries(:,[2]));
% % S_daylight_von_kries = plot(LMS_daylight_von_kries(:,[3]));
% hold off
% grid on


function T_stimuli = LMS_stimuli(t,s,r,n_lamda)
    T_stimuli = transpose(t*s*r*n_lamda)
end
