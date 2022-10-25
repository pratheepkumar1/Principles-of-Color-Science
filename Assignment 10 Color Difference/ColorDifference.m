% Reading the Dataset
warning('off','all');
inkjetColorChecker_dataset = readtable("InkjetColorChecker.xlsx");
source_dataset = readtable("Illuminant Data.xlsx");

% Create a interpolated wavelength info for the datasets
intp_wavelength_info = struct('min',380,'max',750,'range',10);

%Extracting illuminant A and D65 from the source dataset
source_A = interpolateData(source_dataset{:,1},source_dataset.A,intp_wavelength_info);
source_D65 = interpolateData(source_dataset{:,1},source_dataset.D65,intp_wavelength_info);










%Function to interpolate the values to a consistent wavelength
function i = interpolateData(wavelengths,values,intp_wavelength_data)
    %Create the interpolant by passing wavelength and corresponding values to griddedInterpolant.
    GI = griddedInterpolant(wavelengths,values);

    % Create a vector of query points with 5nm spacing for 380 to 780 wavelength.
    wl = intp_wavelength_data.min:intp_wavelength_data.range:intp_wavelength_data.max;

    % Evaluate the interpolant at the each wavelength for each value set 
    i = GI(wl);
end