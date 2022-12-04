%Function to calculate Tristimulus values for materials
function t = calcTristimulus(xyz_value,source,material,d_lambda)   

    %Calculate normalizing constant
    k = (100/(source * xyz_value(:,2) * d_lambda));

    s_lambda = diag(source);

    %Calculating tristimulus
    t = k.*((s_lambda*xyz_value)'*material)*d_lambda;

%   t = custom_normalization(t);
end


