%% Reading the Dataset
warning('off','all');
colDiff_dataset = readtable("RIT_DuPont_Data_2019.xlsx");

lab_std = (colDiff_dataset{:,5:7})';
lab_trl = (colDiff_dataset{:,8:10})';

deltaV_colDiff = (colDiff_dataset{:,4})';


%% Computing ΔEab, ΔE94, ΔE00 color difference metrics fo
% r all pairs
deltaEab_colDiff = deltaEab(lab_trl,lab_std);

deltaE94_colDiff = deltaE94(lab_trl,lab_std);

deltaE00_colDiff = deltaE00(lab_trl,lab_std);

Color_Difference = table(deltaEab_colDiff',deltaE94_colDiff',deltaE00_colDiff','VariableNames', {'ΔEab', 'ΔE94','ΔE00'})


%Table for mean, min and max value of all color difference metrics
mean_deltaEab = mean(deltaEab_colDiff);
mean_deltaE94 = mean(deltaE94_colDiff);
mean_deltaE00 = mean(deltaE00_colDiff);
mean_deltaE = [mean_deltaEab;mean_deltaE94;mean_deltaE00];

min_deltaEab = min(deltaEab_colDiff);
min_deltaE94 = min(deltaE94_colDiff);
min_deltaE00 = min(deltaE00_colDiff);
min_deltaE = [min_deltaEab;min_deltaE94;min_deltaE00];

max_deltaEab = max(deltaEab_colDiff);
max_deltaE94 = max(deltaE94_colDiff);
max_deltaE00 = max(deltaE00_colDiff);
max_deltaE = [max_deltaEab;max_deltaE94;max_deltaE00];

Color_Difference_Summary = table(mean_deltaE,min_deltaE,max_deltaE, 'VariableNames', {'Mean','Min','Max'}, 'RowNames', {'ΔEab', 'ΔE94','ΔE00'})


%% STRESS for each of the three color difference formulas 

STRESS_colDiff_Eab = CalcSTRESS(deltaEab_colDiff,deltaV_colDiff);

STRESS_DuPont_Eab = CalcSTRESS(deltaEab_colDiff(:,1:156),deltaV_colDiff(:,1:156));

STRESS_colDiff_E94 = CalcSTRESS(deltaE94_colDiff,deltaV_colDiff);

STRESS_DuPont_E94 = CalcSTRESS(deltaE94_colDiff(:,1:156),deltaV_colDiff(:,1:156));

STRESS_colDiff_E00 = CalcSTRESS(deltaE00_colDiff,deltaV_colDiff);

STRESS_DuPont_E00 = CalcSTRESS(deltaE00_colDiff(:,1:156),deltaV_colDiff(:,1:156));

STRESS_Values = table([STRESS_colDiff_Eab;STRESS_colDiff_E94;STRESS_colDiff_E00],[STRESS_DuPont_Eab;STRESS_DuPont_E94;STRESS_DuPont_E00], 'VariableNames', {'All Pairs','DuPont Pairs'}, 'RowNames', {'ΔEab', 'ΔE94','ΔE00'})

%% Visualization

figure(1)
hold on

tiledlayout(3,1)
nexttile
scatter(deltaEab_colDiff',deltaV_colDiff,'filled','red','MarkerEdgeColor','black');
hold on
scatter(deltaEab_colDiff(:,1:156)',deltaV_colDiff(:,1:156),'filled','green','MarkerEdgeColor','black');
title('Stress: ΔEab vs ΔV')
legend('All Pairs','DuPont');
xlabel('ΔEab');
ylabel('ΔV');
xlim([0,4]);
ylim([0,4]);
grid on


nexttile
scatter(deltaE94_colDiff',deltaV_colDiff,'filled','red','MarkerEdgeColor','black');
hold on
scatter(deltaE94_colDiff(:,1:156)',deltaV_colDiff(:,1:156),'filled','green','MarkerEdgeColor','black');
title('Stress: ΔE94 vs ΔV');
legend('All Pairs','DuPont');
xlabel('ΔE94');
ylabel('ΔV');
xlim([0,4]);
ylim([0,4]);
grid on

nexttile
scatter(deltaE00_colDiff(:,157:end)',deltaV_colDiff(:,157:end),'filled','red','MarkerEdgeColor','black');
hold on
scatter(deltaE00_colDiff(:,1:156)',deltaV_colDiff(:,1:156),'filled','green','MarkerEdgeColor','black');
title('Stress: ΔE00 vs ΔV');
legend('All Pairs','DuPont');
xlabel('ΔE00');
ylabel('ΔV');
xlim([0,4]);
ylim([0,4]);
grid on

hold off



%% Functions

%Function to calculate Delta Eab
function Deab = deltaEab(lab_bat,lab_ref)
    
    %Chroma
    C_star_ref = C_star(lab_ref(2,:),lab_ref(3,:));
    C_star_bat = C_star(lab_bat(2,:),lab_bat(3,:));

    %Delta Values
    dL = deltaL(lab_bat(1,:),lab_ref(1,:));
    dC = deltaC(C_star_bat,C_star_ref);
    dh = deltah(lab_bat,lab_ref);
    dH = deltaH(C_star_bat,C_star_ref,dh);

    Deab = sqrt(dL.^2 + dC.^2 + dH.^2);
end


%Function to calculate Delta E94
function De94 = deltaE94(lab_bat,lab_ref)


    %Chroma
    C_star_ref = C_star(lab_ref(2,:),lab_ref(3,:));
    C_star_bat = C_star(lab_bat(2,:),lab_bat(3,:));
    C_Star_ab = sqrt(C_star_bat.*C_star_ref);

    %Delta Values
    dL = deltaL(lab_bat(1,:),lab_ref(1,:));
    dC = deltaC(C_star_bat,C_star_ref);
    dh = deltah(lab_bat,lab_ref);
    dH = deltaH(C_star_bat,C_star_ref,dh);

    %Weighting Functions
    SL = 1;
    SC = 1+(0.045.*C_Star_ab);
    SH = 1+(0.015.*C_Star_ab);

    %Parametric factors
    kL = 1;
    kC = 1;
    kH = 1;

    De94 = ((dL./ (kL.*SL)).^2+(dC./ (kC.*SC)).^2+(dH./ (kH.*SH)).^2).^(1/2);

end

%Function to calculate STRESS
function st = CalcSTRESS(delE,delV)
    F = sum(delE.^2)./sum(delE.*delV);
    st = 100*(sum((delE-(F.*delV)).^2)./sum((F.^2).*(delV.^2))).^(1/2);
end

%Function to calculate Delta L
function DL = deltaL(L_bat,L_std)
    DL = L_bat - L_std;
end


%Function to calculate Delta a
function Da = deltaa(a_bat,a_std)
    Da = a_bat - a_std;
end


%Function to calculate Delta b
function Db = deltab(b_bat,b_std)
    Db = b_bat - b_std;
end


function C = C_star(a_star,b_star)
    C = sqrt((a_star).^2 + (b_star).^2);
end


%Function to calculate Delta C*
function DC = deltaC(C_bat,C_std)
    DC = C_bat - C_std;
end


%Function to calculate Delta h
function Dh = deltah(bat,std)
    Dh = hueAngle(bat(2,:),bat(3,:)) - hueAngle(std(2,:),std(3,:));
end


%Function to calculate hue angle
function h = hueAngle(a,b)
      h = atan2(b,a).*180./pi;
      h = h+(h<0).*360;
end


%Function to calculate Delta H
function DH = deltaH(C_bat,C_std,deltah)
    DH = 2*(C_bat.*C_std).^(1/2).*sinDeg(deltah./2);
end


function out = sinDeg(in)
    out = sin(in.*pi./180);
end







