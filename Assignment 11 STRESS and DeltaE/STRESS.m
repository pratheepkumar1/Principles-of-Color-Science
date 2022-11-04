% Reading the Dataset
warning('off','all');
colDiff_dataset = readtable("RIT_DuPont_Data_2019.xlsx");

lab_std = (colDiff_dataset{:,5:7})';
lab_trl = (colDiff_dataset{:,8:10})';

C_star_std = C_star(lab_std(2,:),lab_std(3,:));
C_star_trl = C_star(lab_trl(2,:),lab_trl(3,:));

deltaV_colDiff = (colDiff_dataset{:,4})';

deltaL_colDiff = deltaL(lab_trl(1,:),lab_std(1,:));

deltaa_colDiff = deltaa(lab_trl(2,:),lab_trl(2,:));

deltab_colDiff = deltab(lab_trl(3,:),lab_trl(3,:));

deltaC_colDiff = deltaC(C_star_trl,C_star_std);

deltah_colDiff = deltah(lab_trl,lab_std);

deltaH_colDiff = deltaH(C_star_trl,C_star_std,deltah_colDiff);



deltaEab_colDiff = deltaEab(deltaL_colDiff,deltaC_colDiff,deltaH_colDiff);

deltaE00_colDiff = deltaE00(lab_trl,lab_std);

deltaE94_colDiff = deltaE94(lab_trl,lab_std,C_star_trl,C_star_std);



STRESS_colDiff_E00 = CalcSTRESS(deltaE00_colDiff,deltaV_colDiff)

STRESS_DuPont_E00 = CalcSTRESS(deltaE00_colDiff(:,1:156),deltaV_colDiff(:,1:156))

STRESS_colDiff_Eab = CalcSTRESS(deltaEab_colDiff,deltaV_colDiff)

STRESS_DuPont_Eab = CalcSTRESS(deltaEab_colDiff(:,1:156),deltaV_colDiff(:,1:156))

STRESS_colDiff_E94 = CalcSTRESS(deltaE94_colDiff,deltaV_colDiff)

STRESS_DuPont_E94 = CalcSTRESS(deltaE94_colDiff(:,1:156),deltaV_colDiff(:,1:156))


% r= figure(1);
% hold on
% scatter(deltaL_colDiff(:,1:156)',lab_std(1,1:156)');
% hold off


%% Functions

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
    h = atan2Deg(b,a);
end


%Function to calculate Delta H
function DH = deltaH(C_bat,C_std,deltah)
    DH = 2*(C_bat.*C_std).^(1/2).*sinDeg(deltah./2);
end

function out = atan2Deg(inY,inX)
    out = atan2(inY,inX).*180./pi;
    out = out+(out<0).*360;
end


function out = sinDeg(in)
    out = sin(in.*pi./180);
end


%Function to calculate Delta Eab
function e = deltaEab(deltaL,deltaa,deltab)
    e = sqrt(deltaL.^2 + deltaa.^2 + deltab.^2);
end


%Function to calculate Delta E00
function De94 = deltaE94(lab_bat,lab_ref,C_star_bat,C_star_ref)

    %Delta Values
    deltaL_F = deltaL(lab_bat(1,:),lab_ref(1,:));
    deltaa_F = deltaa(lab_bat(2,:),lab_bat(2,:));
    deltab_F = deltab(lab_bat(3,:),lab_bat(3,:));
    deltaC_F = deltaC(C_star_bat,C_star_ref);
    deltah_F = deltah(lab_bat,lab_ref);
    deltaH_F = deltaH(C_star_bat,C_star_ref,deltah_F);

    C_Star_ab = sqrt(C_star_bat.*C_star_ref);

    %Weighting Functions
    SL = 1;
    SC = 1+(0.045.*C_Star_ab);
    SH = 1+(0.015.*C_Star_ab);

    %Parametric factors
    kL = 1;
    kC = 1;
    kH = 1;

    De94 = ((deltaL_F./ (kL.*SL)).^2+(deltaC_F./ (kC.*SC)).^2+(deltaH_F./ (kH.*SH)).^2).^(1/2);

end


function st = CalcSTRESS(delE,delV)
    F = sum(delE.^2)./sum(delE.*delV);
    st = 100*(sum((delE-(F.*delV)).^2)./sum((F.^2).*(delV.^2))).^(1/2);
end







