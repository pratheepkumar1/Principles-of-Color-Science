% Reading the Dataset
xyz_std_obs_two_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="2-degree");
xyz_std_obs_ten_deg_dataset = readtable("StdObsFuncs.xlsx",Sheet="10-degree");
source_dataset = readtable("Illuminant Data.xlsx");
patch_dataset = readtable("MacbethColorChecker.xlsx");


% Create a interpolated data for the datasets
intp_wavelength_info = struct('min',380,'max',780,'range',5,'values',patch_dataset{2:end,1});

min_wavelength = 380;
max_wavelength = 780;
xyz_std_obs_two_deg = interpolateData(xyz_std_obs_two_deg_dataset{:,1},xyz_std_obs_two_deg_dataset{:,2:end},intp_wavelength_info);
xyz_std_obs_ten_deg = interpolateData(xyz_std_obs_ten_deg_dataset{:,1},xyz_std_obs_ten_deg_dataset{:,2:end},intp_wavelength_info);
sources_d50 = interpolateData(source_dataset{:,1},source_dataset.D50,intp_wavelength_info);
patches = patch_dataset{2:end,2:25};

%Calculate change in Wavlength
wavelength = patch_dataset{2:end,1};
d_lambda = mean(diff(wavelength));

tristimulus_XYZ_two_deg = calcTristimulus(xyz_std_obs_two_deg,sources_d50,patches,d_lambda);
tristimulus_XYZ_ten_deg = calcTristimulus(xyz_std_obs_ten_deg,sources_d50,patches,d_lambda);

% chromacity_coordinates
chromacity_coord_two_deg = calChromacityCoordinates(tristimulus_XYZ_two_deg);
chromacity_coord_ten_deg = calChromacityCoordinates(tristimulus_XYZ_ten_deg);



%Function to interpolate the values to a consistent wavelength
function i = interpolateData(wavelengths,values,intp_wavelength_data)
    %Create the interpolant by passing wavelength and corresponding values to griddedInterpolant.
    GI = griddedInterpolant(wavelengths,values);

    % Create a vector of query points with 5nm spacing for 380 to 780 wavelength.
    wl = intp_wavelength_data.min:intp_wavelength_data.range:intp_wavelength_data.max;

    % Evaluate the interpolant at the each wavelength for each value set 
    i = GI(wl);
end


%Function to calculate Tristimulus of material
function t = calcTristimulus(xyz_value,source,material,d_lambda)   

    %Calculate normalizing constant
    k = 100/(source * xyz_value(:,2) * d_lambda)

    s_lambda = diag(source);

    %Calculating tristimulus
    t = k.*((s_lambda*xyz_value)'*material)*d_lambda;

%   t = custom_normalization(t);
end

%Function to calculate Chromacity Coordinate
function cc = calChromacityCoordinates(tristimulus_XYZ)
    sum_XYZ = tristimulus_XYZ(1,:)+tristimulus_XYZ(2,:)+tristimulus_XYZ(3,:);
    cc = tristimulus_XYZ./sum_XYZ;
end


% Function to normalize a matrix with its peak value (not using default normalize function)
function n = custom_normalization(x)
    max_value = max(x, [], 'all');
    n = x/max_value;
end



