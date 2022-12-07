%% Reading the Dataset
warning('off','all');
color_checker_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="ColorChecker");
lms_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="LMS");
light_sources_dataset = readtable("HW_Opponency_Data.xlsx",Sheet="Sources");

%Change in Wavelength
wavelength = 380:10:730;
d_lambda = mean(diff(wavelength));

% Create LMS relative sensitivity matrix
lms_matrix = transpose(lms_dataset{:,2:4});

% create a diagonal matrix to calculate illuminant's spectral power distribution = S𝜆
s_incandescent_matrix = (normalize(light_sources_dataset{:,2},'range'));
s_daylight_matrix = (normalize(light_sources_dataset{:,3},'range'));

%%%explain normalization is done to light sources


% Objectspectral reflectance factor = R𝜆 
[color_checker_numRows,color_checker_numCols] = size(color_checker_dataset);
r_matrix = color_checker_dataset{:,2:color_checker_numCols};

% LMS values for all 24 ColorChecker samples
LMS_stimuli_incandescent = LMS_stimuli(lms_matrix,diag(s_incandescent_matrix), r_matrix,d_lambda);
LMS_stimuli_daylight = LMS_stimuli(lms_matrix,diag(s_daylight_matrix), r_matrix, d_lambda);


%% Chromatic Adaptation Transformations (CAT)

% Equi-energy light source is Illuminant E (L2)
s_illuminantE = ones(36,1);


% LMS_lightsource
LMS_source_daylight = LMS_source(lms_matrix,s_daylight_matrix,d_lambda);
LMS_source_incandescent = LMS_source(lms_matrix,s_incandescent_matrix,d_lambda);
LMS_source_illuminantE = LMS_source(lms_matrix,s_illuminantE,d_lambda);


% Computing the CAT (Von Kries M Matrix)
vonkries_incandescent = vonkries(LMS_source_incandescent,LMS_source_illuminantE);
vonkries_daylight = vonkries(LMS_source_daylight,LMS_source_illuminantE);


% LMS values of stimuli after applying von kries matrix
LMS_stimuli_incandescent_vk = vonkries_incandescent * LMS_stimuli_incandescent;
LMS_stimuli_daylight_vk = vonkries_daylight * LMS_stimuli_daylight;


% patches_size = 1:1:24;
% 
% figure;
% line(patches_size,LMS_stimuli_incandescent,'LineStyle','-');
% hold on
% line(patches_size,LMS_stimuli_incandescent_vk,'LineWidth',2);
% hold off


%% Opponency Signal

m_opponency = [0.64 0.39 -0.01; 1.12 -1.50 0.34; 0.35 0.15 -0.53];

gamma = 2.4;

opponency_incandescent = calcOpponency(m_opponency,LMS_stimuli_incandescent_vk,gamma);
opponency_daylight = calcOpponency(m_opponency,LMS_stimuli_daylight_vk,gamma);


%% Functions 


% Function to calculate LMS of the stimuli
function T_stimuli = LMS_stimuli(t,s,r,d)
    T_stimuli = (t*s*r*d);
end


% Function to calculate LMS of the light source
function T_source = LMS_source(t,s,d)
    T_source = (t*s*d);
end


% Function to calculate Chromatic Adaptation Transformation - Von Kries
function M_factor = vonkries(LMS1,LMS2)
%     M_factor = diag([LMS2(1,:)\LMS1(1,1) ; LMS2(:,2)\LMS1(:,2) ;...
%         LMS2(:,3)\LMS1(:,3)]);
    M_factor = diag(LMS2./LMS1);
end


% Function to calculate Opponency
function opp = calcOpponency(m_opponency, LMS_vk, gamma)
    opp = m_opponency * ((LMS_vk).^(1/gamma));
end



function n = custom_normalization(x)
    max_value = max(x, [], 'all');
    n = x/max_value;
end
