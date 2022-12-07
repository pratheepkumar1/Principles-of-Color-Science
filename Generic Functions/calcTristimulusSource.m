%Function to calculate Tristimulus values for sources (Whitepoint)
function ts = calcTristimulusSource(xyz_value,source,d_lambda)   

    %Calculate normalizing constant
    k = 100/(source * xyz_value(:,2) * d_lambda);

    %Calculating tristimulus
    ts = k.*((source*xyz_value))*d_lambda;

%   t = custom_normalization(t);
end