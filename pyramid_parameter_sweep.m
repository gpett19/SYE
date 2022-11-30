clear
close all

disp(' ')
disp('test_pyramids.m')

% Define parameters for grating structure and incident field:
d1 = 2; % grating period in x2 direction
d2 = 2; % grating period in x3 direction
%Radius of pyrs
r = 1;


h = 0.1:0.1:10; % grating height
wlr = 14;
%h = 1;
%wavelength = 1.55;
phi = 0; % incidence polar angle from x1 axis, deg
psi = 0; % incidence azimuthal angle around x1 axis (zero on x2 axis), deg
grating_pmt = 3.42^2; % grating permittivity
L1 = 50; % number of grating strata
m_max = 2; % maximum diffraction order index




% Specify which diffracted orders are to be retained in the calculations,
% Eq's 4.12-13 in GD-Calc.pdf.
order = [];
for m2 = -m_max:m_max
    order(end+1).m2 = m2; %#ok<SAGROW> 
    order(end).m1 = -m_max:m_max;
end

R_list = cell(length(h), length(d1));
T_list = cell(length(h), length(d1));

%Count percentages...
count = 1;
%Total number of possible combos
totNum = length(h) * length(d1) * length(r) * length(wlr);

for i = 1:length(h)
    for j = 1:length(d1)
        for k = 1:length(r)
% Construct grating.
clear grating
grating.pmt = {1,grating_pmt}; % grating material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
grating.d21 = d1(j); % first grating period: x2 projection
grating.d31 = 0; % first grating period: x3 projection
grating.d22 = 0; % second grating period: x2 projection
grating.d32 = d2(j); % second grating period: x3 projection
grating.stratum = {};
clear stratum
stratum.type = 2; % biperiodic stratum
stratum.thick = h(i)/L1; % stratum thickness
% The following h11 ... h22 spec indicates that the stratum's period
% vectors match the grating periods (GD-Calc.pdf, Eq. 3.20).
stratum.h11 = 1;
stratum.h12 = 0;
stratum.h21 = 0;
stratum.h22 = 1;
clear stripe
stripe.type = 0; % first stripe is homogeneous
stripe.pmt_index = 1; % first stripe's permittivity index
stratum.stripe{1} = stripe;
clear stripe
stripe.type = 1; % second stripe is inhomogeneous
stripe.block{1}.pmt_index = 1; % first block's permittivity index
stripe.block{2}.pmt_index = 2; % second block's permittivity index
stratum.stripe{2} = stripe;
clear stripe
for l1 = 1:L1
    % Set first and second stripes' boundaries (c1) on positive side:
    stratum.stripe{1}.c1 = -((L1-l1+0.5)/L1)/2*r(k);
    stratum.stripe{2}.c1 = -stratum.stripe{1}.c1;
    % Set first and second block' boundaries (c2) on positive side:
    stratum.stripe{2}.block{1}.c2 = -((L1-l1+0.5)/L1)/2*r(k);
    stratum.stripe{2}.block{2}.c2 = -stratum.stripe{2}.block{1}.c2;
    grating.stratum{end+1} = stratum;
end
clear stratum

                for l = 1:length(wlr)

                    wavelength = wlr(l);

% Define the indicent field.
clear inc_field
inc_field.wavelength = wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector, Eq's 4.1-2 in
% GD-Calc.pdf. The grating-normal (x1) projection is implicitly
%   f1 = -cosd(phi)/wavelength
inc_field.wavelength = wavelength;
inc_field.f2 = sind(phi)*cosd(psi)./wavelength;
inc_field.f3 = sind(phi)*sind(psi)./wavelength;


% Run the diffraction calculations.
tic
[~,scat_field,inc_field] = gdc(grating,inc_field,order,false);
toc


% Compute the diffraction efficiencies.
[R,T] = gdc_eff(scat_field,inc_field);

R_list{i,j, k, l} = R;
T_list{i,j, k, l} = T;

pct = count / totNum * 100;

disp("Progress = " + string(pct) + "%")
count = count + 1;

            end

        end
    end
end

