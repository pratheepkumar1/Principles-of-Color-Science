color_checker_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="ColorChecker");
lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
light_sources_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="Sources");
incandescent_source = light_sources_dataset{:,2};
daylight_source= light_sources_dataset{:,3};


lms_matrix = transpose(lms_dataset{:,[2:4]});

% create a diagonal matrix to calculate illuminant's spectral power distribution = S𝜆
s_incandescent_matrix = diag(incandescent_source);
s_daylight_matrix = diag(daylight_source);

% Objectspectral reflectance factor = R𝜆 
[color_checker_numRows,color_checker_numCols] = size(color_checker_dataset);
r_matrix = color_checker_dataset{:,[2:color_checker_numCols]};

LMS_incandescent = lms_matrix*s_incandescent_matrix*r_matrix;
