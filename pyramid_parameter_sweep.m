clear
close all

disp(' ')
disp('test_pyramids.m')

% Define parameters for grating structure and incident field:
d1 = 1.13; % grating period in x2 direction
d2 = 1.13; % grating period in x3 direction
%Radius of pyrs
r = 1;


h = 0.9; % grating height
wlr = 1.4;
%h = 1;
%wavelength = 1.55;
phi = 0:2:80; % incidence polar angle from x1 axis, deg
psi = 0; % incidence azimuthal angle around x1 axis (zero on x2 axis), deg
grating_pmt = 1.444^2; % grating permittivity
L1 = 100; % number of grating strata
m_max = 4; % maximum diffraction order index






R_list = cell(length(h), length(d1));
T_list = cell(length(h), length(d1));

%Count percentages...
count = 1;
%Total number of possible combos
totNum = length(h) * length(d1) * length(r) * length(wlr);

% Specify which diffracted orders are to be retained in the calculations,
% Eq's 4.12-13 in GD-Calc.pdf.
order = [];
for m2 = -m_max:m_max
    order(end+1).m2 = m2; %#ok<SAGROW> 
    order(end).m1 = - m_max:m_max;
end


for i = 1:length(h)
    for j = 1:length(d1)
        for k = 1:length(r)

%List of possible slant scales
sSL = 1;
%Loop over said slant scales.
for sS = 1:length(sSL)


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

%Scales for the slating of the structure
slantToggle = [1,0]; %Toggle for slanting; set to 1 if want slanting in that direction, 0 otherwise
slantScale = [sSL(sS), sSL(sS)]; %The slant scales themselves; 0.5 gives vertical boundaries.
slantScale = slantScale .* slantToggle; %"Deactivates" the unwanted slants
for l1 = 1:L1
    % Set first and second stripes' boundaries (c1) on positive side:
    stratum.stripe{1}.c1 = -(((L1-l1+0.5)/L1)/2 - slantScale(1)*l1/L1)*r(k);
    stratum.stripe{2}.c1 = -(stratum.stripe{1}.c1 - 2*slantScale(1)*l1/L1); %Need 2x factor to effectively shift it to the right.
    % Set first and second block' boundaries (c2) on positive side:
    stratum.stripe{2}.block{1}.c2 = -(((L1-l1+0.5)/L1)/2 - slantScale(2)*l1/L1)*r(k);
    stratum.stripe{2}.block{2}.c2 = -(stratum.stripe{2}.block{1}.c2 - 2*slantScale(2)*l1/L1);
    grating.stratum{end+1} = stratum;
end
%OLD LOOP
%{
for l1 = 1:L1
    % Set first and second stripes' boundaries (c1) on positive side:
    stratum.stripe{1}.c1 = -((L1-l1+0.5)/L1)/2*r(k);
    stratum.stripe{2}.c1 = -stratum.stripe{1}.c1;
    % Set first and second block' boundaries (c2) on positive side:
    stratum.stripe{2}.block{1}.c2 = -((L1-l1+0.5)/L1)/2*r(k);
    stratum.stripe{2}.block{2}.c2 = -stratum.stripe{2}.block{1}.c2;
    grating.stratum{end+1} = stratum;
end
%}
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
%[~,scat_field,inc_field] = gdc(grating,inc_field,order,false);
toc


% Compute the diffraction efficiencies.
%[R,T] = gdc_eff(scat_field,inc_field);

%R_list{i,j, k, l} = R;
%T_list{i,j, k, l} = T;

pct = count / totNum * 100;

disp("Progress = " + string(pct) + "%")
count = count + 1;

            end
end
        end
    end
end

% Plot the grating.
clear pmt_display
pmt_display(1).name = '';

pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.75,.75,.75];
pmt_display(2).alpha = 1;
x_limit = [-0.5*h,-1.5*d1,-1.5*d2;1.5*h,1.5*d1,1.5*d2];
h_plot = gdc_plot(grating,1,pmt_display,x_limit);
view(-25,15)
