color_checker_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="ColorChecker");
lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
light_sources_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="Sources");

% Create LMS relative sensitivity matrix
lms_matrix = transpose(lms_dataset{:,[2:4]});


[light_sources_numRows,light_sources_numCols] = size(light_sources_dataset); 
normalize_s_factor = 100*ones(light_sources_numRows,1); %normalise the matrix values to 0-1
incandescent_source = light_sources_dataset{:,2}./ normalize_s_factor;
daylight_source= light_sources_dataset{:,3} ./ normalize_s_factor;


% create a diagonal matrix to calculate illuminant's spectral power distribution = S𝜆
s_incandescent_matrix = diag(incandescent_source);
s_daylight_matrix = diag(daylight_source);

% Objectspectral reflectance factor = R𝜆 
[color_checker_numRows,color_checker_numCols] = size(color_checker_dataset);
r_matrix = color_checker_dataset{:,[2:color_checker_numCols]};

% LMS values for all 24 ColorChecker samples
LMS_incandescent = transpose(lms_matrix * s_incandescent_matrix * r_matrix)
LMS_daylight = transpose(lms_matrix * s_daylight_matrix * r_matrix);

% Chromatic Adaptation Transformations (CAT)
white_col_number = find(string(color_checker_dataset.Properties.VariableNames) == "White");
LMS_white_incandescent = transpose(LMS_incandescent(white_col_number,:));
LMS_white_daylight = transpose(LMS_daylight(white_col_number,:));

degree_of_adaptation = 0.8; %d_factor

CAT_von_kries_incadescent_factor = diag((degree_of_adaptation*LMS_white_incandescent)+...
    ((1-degree_of_adaptation)*LMS_white_daylight))./(LMS_white_daylight);

CAT_von_kries_daylight_factor = diag((degree_of_adaptation*LMS_white_daylight)+...
    ((1-degree_of_adaptation)*LMS_white_incandescent))./(LMS_white_incandescent);

LMS_daylight_von_kries = transpose(CAT_von_kries_daylight__factor * transpose(LMS_incandescent));
LMS_incandescent_von_kries = transpose(CAT_von_kries_incadescent_factor * transpose(LMS_daylight))


hold on
L_daylight = plot(LMS_daylight(:,[1]));
% M_daylight = plot(LMS_daylight(:,[2]));
% S_daylight = plot(LMS_daylight(:,[3]));
L_daylight_von_kries = plot(LMS_daylight_von_kries(:,[1]));
% M_daylight_von_kries = plot(LMS_daylight_von_kries(:,[2]));
% S_daylight_von_kries = plot(LMS_daylight_von_kries(:,[3]));
hold off
grid on