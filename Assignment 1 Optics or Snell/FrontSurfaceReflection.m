n_material = 1.46;
n_air = 1;
surface_reflection = zeros(9,1);
angle_incidence = linspace(0,80,9);

for i = 1 : length(angle_incidence)
    r_parallel = ((cosd(angle_incidence(i)) - ...
        sqrt((n_material/n_air)^2 - sind(angle_incidence(i) ^ 2))) ...
        / ...
        (cosd(angle_incidence(i)) + ...
        sqrt((n_material/n_air)^2 - sind(angle_incidence(i) ^ 2))))^2 ;
    
    
    r_perpendicular = ((((n_material/n_air)^2)*(cosd(angle_incidence(i))) - ...
        sqrt((n_material/n_air)^2 - sind(angle_incidence(i) ^ 2))) ...
        / ...
        (((n_material/n_air)^2)*(cosd(angle_incidence(i))) + ...
        sqrt((n_material/n_air)^2 - sind(angle_incidence(i) ^ 2))))^2 ;
    
   surface_reflection(i) = ((r_parallel + r_perpendicular ) / 2) * 100;
end

disp(surface_reflection)
h = figure;
plot(angle_incidence,surface_reflection)

datacursormode(h,'on');
xlabel('Angle of Incidence (in degrees)')
ylabel('Surface Reflection (in %)')
title('Front Surface Reflection')
ylim([0 100])