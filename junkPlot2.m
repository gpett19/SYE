%This plots the effective refractive indices and fill factors of cones and
%pyramids from cone_tf_sweep_1 and pyr_tf_sweep_1. (Note you'll have to
%rename your ffs variables to keep them all loaded!)

close
yyaxis left

plot(ffs_cone, '-b')
hold on
%plot(ffs_pyr, '-b')

xlabel("Layer Number")
ylabel("Silica Fill Factor")

%Now, get the neffs

nGrating = 1.444;

for i = 1:length(ffs_cone)
    neff1_cone(i) = ffs_cone(i)*nGrating + (1-ffs_cone(i))*1;
    
    neff2_cone(i) = sqrt(ffs_cone(i)*nGrating^2 + (1^2)*(1-ffs_cone(i)));
end


%for i = 1:length(ffs_pyr)
%    neff1_pyr(i) = ffs_pyr(i)*nGrating + (1-ffs_pyr(i))*1;
%    
%    neff2_pyr(i) = sqrt(ffs_pyr(i)*nGrating^2 + (1^2)*(1-ffs_pyr(i)));
%end

yyaxis right

plot(neff1_cone, '--b')
hold on
%plot(neff1_pyr, '--r')

plot(neff2_cone, '-.b')
%hold on
%plot(neff2_pyr, '-.r')

ylabel("Effective Refractive Index")

%legend("Fill Factor (cones)", "Fill Factor (pyramids)", "n_{eff} (cones, Approx. 1)", "n_{eff} (pyrs, Approx. 1)", "n_{eff} (cones, Approx. 2)", "n_{eff} (pyrs, Approx. 2)")


%legend("Fill Factor (pyramids)", "n_{eff} (pyrs, Approx. 1)", "n_{eff} (pyrs, Approx. 2)")
legend("Fill Factor", "Index Approximation", "Permittivity Approximation")


title("Fill Factor and n_{eff} vs. Layer Index")

%%

%This section plots the actual structure sim (R_list) vs the thin film
%sims (R_list1 & R_list2). For theta.
close all
figure()
hold on
lsts = {R_list, R_list1, R_list2};

%Get the fresnel reflections:
n1 = 1;
n2 = 1.444;
Rs = abs(((n1*cosd(theta)-n2*sqrt(1-((n1/n2)*sind(theta)).^2))./(n1*cosd(theta)+n2*sqrt(1-(n1/n2*sind(theta)).^2)))).^2;
Rp = abs((n1*sqrt(1-(n1/n2 * sind(theta)).^2)-n2*cosd(theta))./(n1*sqrt(1-(n1/n2.*sind(theta)).^2)+n2.*cosd(theta))).^2;
Rtot = (Rs + Rp) ./ 2;

for i = 1:length(lsts)

    R_list_pre = preprocess_R_list(lsts{i});
    R_list_pre = {R_list_pre{3, 1, 1, 1, :}};
    Rmats{i} = Rl_to_mat(R_list_pre);

    %Plot after dividing out the base reflectance!
    plot(theta, Rmats{i} ./ Rtot)
    %plot(theta, Rmats{i}, '--')
end
%plot(theta, Rtot, '-.')

xlabel("Incidence Angle (Deg)")

ylabel("Reflectance (Normalized to Fresnel Reflectance)")

title("Normalized Reflectance vs Angle")

legend("Structure Simulation", "Linear Approximation", "Quadratic Approximation")

%%

%This section plots the actual structure sim (R_list) vs the thin film
%sims (R_list1 & R_list2). For Everything Else!.
close all
figure(1)
hold on
figure(2)
hold on
lsts = {R_list, R_list1, R_list2};

%Get the fresnel reflections:
n1 = 1;
n2 = 1.444;
%Rs = abs(((n1*cosd(theta)-n2*sqrt(1-((n1/n2)*sind(theta)).^2))./(n1*cosd(theta)+n2*sqrt(1-(n1/n2*sind(theta)).^2)))).^2;
%Rp = abs((n1*sqrt(1-(n1/n2 * sind(theta)).^2)-n2*cosd(theta))./(n1*sqrt(1-(n1/n2.*sind(theta)).^2)+n2.*cosd(theta))).^2;
%Rtot = (Rs + Rp) ./ 2;

for i = 1:length(lsts)

    R_list_pre = preprocess_R_list(lsts{i});
    R_list_pre = {R_list_pre{3, 1, 1, 1, :}};
    Rmats{i} = Rl_to_mat(R_list_pre);

    %Plot after dividing out the base reflectance!
    figure(1)
    plot(theta, Rmats{i})
    figure(2)
    if i ~= 1 %Don't plot the percent difference with itself
        plot(theta, 200*(Rmats{i} - Rmats{1})./(Rmats{i} + Rmats{1}))
    end
    %plot(theta, Rmats{i}, '--')
end
%plot(theta, Rtot, '-.')
figure(1)
xlabel("Wavelength (\mu m)")

ylabel("Reflectance")

title("Normalized Reflectance vs Wavelength")

legend("Structure Simulation", "Linear Approximation", "Quadratic Approximation")

figure(2)
xlabel("Angle (Deg)")

ylabel("Percent Difference")

title("Percent Difference vs Structure Simulation")

legend("Linear Approximation", "Quadratic Approximation")

%% Plots the fill factors for pyramids and bullets.

for l1 = 1:L1
    ff1(l1) = ((2/L1^3)*l1^3 - (3/L1^2)*l1^2 + 1);
    ff2(l1) = ((L1-l1+0.5)/L1)^2;
end

close
plot(1:L1, ff1)
hold on
plot(1:L1, ff2)

xlabel("Layer Index")
ylabel("Fill Factor")
title("Fill Factor vs. Layer Index")
legend("Bullet-shaped Structures", "Linear Pyramid Structures")


%%
%This section plots the actual structure sim (R_list) for a bunch of
%different structure shapes but same parameters!
close all
figure()
hold on
lsts = {R_list_bul, R_list_cone, R_list_pyr};

%Get the fresnel reflections:
n1 = 1;
n2 = 1.444;
Rs = abs(((n1*cosd(theta)-n2*sqrt(1-((n1/n2)*sind(theta)).^2))./(n1*cosd(theta)+n2*sqrt(1-(n1/n2*sind(theta)).^2)))).^2;
Rp = abs((n1*sqrt(1-(n1/n2 * sind(theta)).^2)-n2*cosd(theta))./(n1*sqrt(1-(n1/n2.*sind(theta)).^2)+n2.*cosd(theta))).^2;
Rtot = (Rs + Rp) ./ 2;

for i = 1:length(lsts)

    R_list_pre = preprocess_R_list(lsts{i});
    R_list_pre = {R_list_pre{1, 1, :}};
    Rmats{i} = Rl_to_mat(R_list_pre);

    %Plot after dividing out the base reflectance!
    plot(theta, Rmats{i} ./ Rtot)
    %plot(theta, Rmats{i}, '--')
end
%plot(theta, Rtot, '-.')

xlabel("Incidence Angle (Deg)")

ylabel("Reflectance (Normalized to Fresnel Reflectance)")

title("Normalized Reflectance vs Angle")

legend("Bullets", "Cones", "Pyramids")

%%

close all
plot(bullet_ffs, "Color", "#0072BD")
hold on
plot(pyr_ffs, "Color", "#D95319")
plot(cone_ffs, "Color", "#7E2F8E")
%plot(cyl_ffs)

xlabel("Layer Number")
ylabel("Silica Fill Factor")
title("Fill Factors vs Layer Index")

yyaxis right

plot(neff1_bullet, '--', "Color", "#0072BD")
plot(neff1_pyr, '--', "Color", "#D95319")
plot(neff1_cone, '--', "Color", "#7E2F8E")

ylabel("Effective Refractive Index")
legend("Bullets", "Pyramids", "Cones", "Bullet Effective Index", "Pyramid Effective Index", "Cone Effective Index")


%%
%neffs
nGrating = sqrt(1.444^2);
%Then, calculate neff
%Assuming n = 1 for air refractive index
neff1_bullet = bullet_ffs.*nGrating + (1-bullet_ffs)*1;
neff2_bullet = sqrt(bullet_ffs.*nGrating^2 + (1^2)*(1-bullet_ffs));

neff1_cone = cone_ffs.*nGrating + (1-cone_ffs)*1;
neff2_cone = sqrt(cone_ffs.*nGrating^2 + (1^2)*(1-cone_ffs));

neff1_pyr = pyr_ffs.*nGrating + (1-pyr_ffs)*1;
neff2_pyr = sqrt(pyr_ffs.*nGrating^2 + (1^2)*(1-pyr_ffs));
