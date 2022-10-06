% Reading the Dataset
xyz_std_obs_two_degree_dataset = readtable("StdObsFuncs.xlsx",Sheet="2-degree");
xyz_std_obs_ten_degree_dataset = readtable("StdObsFuncs.xlsx",Sheet="10-degree");
source_dataset = readtable("Illuminant Data.xlsx");
patch_dataset = readtable("MacbethColorChecker.xlsx");


% Create a interpolated data for the datasets
xyz_std_obs_two_degree = interpolateData(xyz_std_obs_two_degree_dataset{:,1},xyz_std_obs_two_degree_dataset{:,2:end});
xyz_std_obs_ten_degree = interpolateData(xyz_std_obs_ten_degree_dataset{:,1},xyz_std_obs_ten_degree_dataset{:,2:end});
sources_d50 = interpolateData(source_dataset{:,1},source_dataset.D50);
patches = patch_dataset{2:end,2:25};


%Change in Wavlength
wavelength = patch_dataset{2:end,1};
d_lambda = mean(diff(wavelength));













%Function to interpolate the values to a consistent wavelength
function i = interpolateData(wavelengths,values)
    %Create the interpolant by passing wavelength and corresponding values to griddedInterpolant.
    GI = griddedInterpolant(wavelengths,values);

    % Create a vector of query points with 5nm spacing for 380 to 780 wavelength.
    wl = 380:5:780;

    % Evaluate the interpolant at the each wavelength for each value set 
    i = GI(wl);
end




