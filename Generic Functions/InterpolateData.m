%Function to interpolate the values to a consistent wavelength
function i = InterpolateData(wavelengths,values,intp_wavelength_data)
    %Create the interpolant by passing wavelength and corresponding values to griddedInterpolant.
    GI = griddedInterpolant(wavelengths,values);

    % Create a vector of query points with 5nm spacing for 380 to 780 wavelength.
    wl = intp_wavelength_data.min:intp_wavelength_data.range:intp_wavelength_data.max;

    % Evaluate the interpolant at the each wavelength for each value set 
    i = GI(wl);
end
