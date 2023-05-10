clear
close all

disp(' ')
disp('test_pyramids.m')

% Define parameters for grating structure and incident field:
d1 = 1.09; % grating period in x2 direction
d2 = 1.09; % grating period in x3 direction
%Radius of pyrs
r = 1;


h = 0.9; % grating height
wlr = 1.5;
%h = 1;
%wavelength = 1.55;
phi = 0:2:80; % incidence polar angle from x1 axis, deg
psiRng = 0:60:360; % incidence azimuthal angle around x1 axis (zero on x2 axis), deg
grating_pmt = 1.444^2; % grating permittivity
L1 = 40; % number of grating strata
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


for rn = 1
%Right, we now want to have several different "rows" and/or columns that all have different
%heights, with heights defined as a percentage of the maximum height,
%specified by parameter h.

%We'll do this by defining each stratum as a rowNum x rowNum "grid", and
%inserting the appropriate heights.

rowNum = rn;
heights = linspace(0, 1, rowNum + 1);
heights = heights(2: end);


if (length(heights) ~= rowNum)
    error("Error: Different number of heights than number of rows.");
end


for i = 1:length(h)
    for j = 1:length(d1)
        for k = 1:length(r)

% Construct grating.
clear grating
grating.pmt = {1,grating_pmt}; % grating material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
grating.d21 = d1(j)*rowNum; % first grating period: x2 projection
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

%We define the borders between each pyramid structure.
borders = linspace(-0.5, 0.5, rowNum + 1);

%The centers are in halfway between each border.
for b = 1:length(borders)-1
    centers(b) = (borders(b) + borders(b+1) )/2;
end


for l1 = 1:L1

    stripeIdx = 1;
    for ri = 1:rowNum
        hMax = floor(L1*heights(ri));
        %Define the stripes
        %Air stripe
        clear stripe
        stripe.type = 0; % first stripe is homogeneous
        stripe.pmt_index = 1; % first stripe's permittivity index
        stratum.stripe{stripeIdx} = stripe;
        clear stripe
        %Block stripe. We index the blocks by stripeIdx as well...
        stripe.type = 1; % second stripe is inhomogeneous
        stripe.block{1}.pmt_index = 1; % first block's permittivity index
        stripe.block{2}.pmt_index = 2; % second block's permittivity index
        stratum.stripe{stripeIdx + 1} = stripe;
        clear stripe
        % Set first and second stripes' boundaries (c1) on positive side:
        rad = ((hMax-l1+0.5)/hMax)/(2*rowNum);
        if rad < 0
            rad = 0;
            %We'll also set the block stripe material to air
            %TODO This is dumb and inefficient.
            stratum.stripe{stripeIdx + 1}.block{2}.pmt_index = 1;            
        end

        stratum.stripe{stripeIdx}.c1 = centers(ri)-rad;
        stratum.stripe{stripeIdx + 1}.c1 = centers(ri)+rad; 
        % Set first and second block' boundaries (c2) on positive side:
        blockRad = ((hMax-l1+0.5)/hMax)/2;
        if blockRad < 0
            blockRad = 0;
        end
        stratum.stripe{stripeIdx + 1}.block{1}.c2 = -blockRad;
        stratum.stripe{stripeIdx + 1}.block{2}.c2 = -stratum.stripe{stripeIdx + 1}.block{1}.c2;
        stripeIdx = stripeIdx + 2;
     end
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
for psi = 1:length(psiRng)
inc_field.wavelength = wavelength;
inc_field.f2 = sind(phi)*cosd(psiRng(psi))./wavelength;
inc_field.f3 = sind(phi)*sind(psiRng(psi))./wavelength;

%{'
% Run the diffraction calculations.
tic
%[~,scat_field,inc_field] = gdc(grating,inc_field,order,false);
toc


% Compute the diffraction efficiencies.
%[R,T] = gdc_eff(scat_field,inc_field);

%R_list{i,j, k, l, rn, psi} = R;
%T_list{i,j, k, l, rn, psi} = T;

pct = count / totNum * 100;

disp("Progress = " + string(pct) + "%")
count = count + 1;
end
%}'


% Plot the grating.
clear pmt_display
pmt_display(1).name = '';

pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.75,.75,.75];
pmt_display(2).alpha = 1;
x_limit = [-0.5*h,-4*d1,-1.5*d2;1.5*h,4*d1,1.5*d2];
h_plot = gdc_plot(grating,1,pmt_display,x_limit);
view(-25,15)


                end
        end
    end
end
end
