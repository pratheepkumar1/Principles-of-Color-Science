% prompt = "Enter n_material = ";
% n_material = input(prompt);
% prompt = "Enter n_air = ";
% n_air = input(prompt);
% prompt = "Enter angle of incidence (in deg) = ";
% angle_incidence = input(prompt);

n_material = [1.6410 1.6226 1.6175];
n_air = 1;
angle_incidence = 55;

for i = 1 : length(n_material)

r_parallel = ((cosd(angle_incidence) - ...
    sqrt((n_material(i)/n_air)^2 - sind(angle_incidence)^2)) ...
    / ...
    (cosd(angle_incidence) + ...
    sqrt((n_material(i)/n_air)^2 - sind(angle_incidence)^2)))^2 ;


r_perpendicular = ((((n_material(i)/n_air)^2)*(cosd(angle_incidence)) - ...
    sqrt((n_material(i)/n_air)^2 - sind(angle_incidence)^2)) ...
    / ...
    (((n_material(i)/n_air)^2)*(cosd(angle_incidence)) + ...
    sqrt((n_material(i)/n_air)^2 - sind(angle_incidence)^2)))^2 ;

surface_reflection = ((r_parallel + r_perpendicular ) / 2) * 100;
disp(surface_reflection)

end


% Calculations
% a = cosd(angle_incidence)
% b = (n_material/n_air)^2
% c = sind(angle_incidence^2)
% pl_numerator = a - sqrt(b-c)
% pl_denominator = a + sqrt(b-c)
% r_pl = (pl_numerator/pl_denominator)^2
% 
% pr_numerator = b*a - sqrt(b-c)
% pr_denominator = b*a + sqrt(b-c)
% r_pr = (pr_numerator/pr_denominator)^2
% 
% s_ref = (r_pl + r_pr)/2