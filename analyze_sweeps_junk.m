%Really junky way of analyzing sweeps just to get some pretty graphs...
%Gonna do pyramids first:
%clearvars -except R_list T_list
close all


%The first thing we'll do is look at the zeroth diffraction order of the
%scattered field, and "average" the polarization to get a sense of how this
%structure would behave under natural light. This may not be the best
%approach in the long run, since the long-term goal is to fabricate these
%on the back of one of the WGPs, so we'll necessarily have light of a given
%polarization incident on these structures.
%   -> Note that this could imply that we might even want a radially
%   asymmetric structure! This would require much further thought, though.

[szx, szy] = size(R_list); %Get the size of the cell array to find how many we need to iterate over
Rl_no_pol = cell(szx, szy);

for i = 1:szx
    for j = 1:szy
            %Just get the zeroth diffraction efficiency for plotting...
            
            tmp = R_list{i,j}(extractfield(R_list{i,j}, 'm1') == 0);
            tmp2 = tmp(extractfield(tmp, 'm2') == 0);
            %Average to get reflectance for natural light.
            Rl_no_pol{i,j} = (tmp2.eff1 + tmp2.eff2)./2;
    end
end


%Alright, we finally have the data in a good way! Thank god...
%So, we'll actually make some preliminary plots now, and get the better
%plotting function working later...

%Plot reflectance vs angle, with lines denoting structure height.
plot(0:5:90, cell2mat(Rl_no_pol(:, 5)));
colororder(turbo(10))
xlabel("Angle (degrees)")
ylabel("Reflectance")
title("Reflectance vs Angle for varying height (d = 0.41 \mu m)")
lgd = legend(string(0.01:0.1:1), "Location","northwest");
title(lgd, "Structure Height (\mu m)")

%Some weirdness with the pyramids is that changing the grating period did
%not change the behaviour at all...
%Not sure why this is... Probably something with the grating definition.

%Anyways, now that we've got that working with the pyramids, let's do it
%with the cylinder data...

%%
%clear
close all;

%load("cyl_param_sweeo_2.mat")

[szx, szy, szz] = size(R_list); %Get the size of the cell array to find how many we need to iterate over
Rl_no_pol = cell(szx, szy, szz);

%Again, loop through an extract the zeroth-order efficiency at averaged
% polarization...
for i = 1:szx
    for j = 1:szy
        for k = 1:szz
            %Just get the zeroth diffraction efficiency for plotting...
            %Handle null/zero values.
            if ~isa(R_list{i,j,k}, 'double')
            tmp = R_list{i,j,k}(extractfield(R_list{i,j,k}, 'm1') == 0);
            tmp2 = tmp(extractfield(tmp, 'm2') == 0);
            Rl_no_pol{i,j,k} = (tmp2.eff1 + tmp2.eff2)./2;
            end
        end
    end
end

%Note that I kind of borked the data saving procedure, so our radii are
%actually *reversed* from the r array...


%Alright, we finally have the data in a good way! Thank god...
%So, we'll actually make some preliminary plots now, and get the better
%plotting function working later...

%Plot reflectance vs angle, with lines denoting structure radius.
%Height of 0.51, distance of 0.91 
plot(0:5:90, squeeze(cell2mat(Rl_no_pol(5, 10, :))));
colororder(turbo(10))
xlabel("Angle (degrees)")
ylabel("Reflectance")
title("Reflectance vs Angle for varying radius (h = 0.41 \mu m, d = 0.91 \mu m, \lambda = 1550nm)")
lgd = legend(string(r{10}), "Location","northwest");
title(lgd, "Structure Radius (\mu m)")

%Angle graphs don't really give much information, since the angle
%dependency doesn't much change with the params... (Of course,
%*comparative* angle graphs might be good!)

%So, plot reflectance at normal incidence vs radius...
%Something smarter in analyzing the data would have been to choose radius
%as a set of proportions, so we get the same number of radii every time, 

%Convert the Rl_no_pol cell array to a 3-D Matrix, with maximal distance (0.91 um)
%cell2mat is acting up, so just doing it the old-fashioned way...\
rlnps = squeeze(Rl_no_pol(:,10,:));
rlnpm = [];
[szx, szy] = size(rlnps);
for i = 1:szx
    for j = 1:szy
        for k = 1:length(rlnps{i,j})
            rlnpm(i,j,k) = rlnps{i,j}(k);
        end
    end
end

%So, we're now indexing by (height, radius, angle)
%Let's keep the angle constant at 0

%Plot the reflectance vs radius, lines denoting height
figure();

plot(r{10}, rlnpm(:,:,1))
colororder(turbo(10))
xlabel("Radius (\mu m)")
ylabel("Reflectance")
title("Reflectance vs Radius for varying height (d = 0.91 \mu m), \theta = 0^\circ, \lambda = 1550nm")
lgd = legend(string(h), "Location","southeast");
title(lgd, "Structure Height (\mu m)")

%Now, plot the opposite (ref vs h, lines denote r)...
figure()

plot(h, rlnpm(:,:,1))
colororder(turbo(8))
xlabel("Height (\mu m)")
ylabel("Reflectance")
title("Reflectance vs Height for varying radius(d = 0.91 \mu m), \theta = 0^\circ, \lambda = 1550nm")
lgd = legend(string(r{10}), "Location","southeast");
title(lgd, "Structure Radius (\mu m)")


%%
R_list_cpy = R_list;
[szx, szy, szz] = size(R_list_cpy);
for i = 1:szx
    for j = 1:szy
        for k = 1:szz
            %Check for empty 
            if ~isa(R_list_cpy{i,j,k}, 'double')
                lst = expandRstruct(R_list_cpy{i,j,k});
                for l = 1:length(lst)
                    R_list_cpy{i,j,k,l} = lst{l};R
                end 
            end
        end
    end
end



%This function takes in a single R struct and exapands it out to a cell
%array of structs.
function cellvec = expandRstruct(R)
    %Initial parameter definition:
    numStruct = length(R(1).eff1);
    cellvec = cell(numStruct,1);
    %tempR = struct(length(R));

    %This loop will extract only the ith efficiencies but otherwise
    %basically copy the struct.
    %Loop through each vector entry
    for i = 1:numStruct
        %...for each diffraction order.
        for j = 1:length(R)
            tempR(j).m1 = R(j).m1;
            tempR(j).m2 = R(j).m2;
            %Repeating self but eh
            tempR(j).eff1 = R(j).eff1(i);
            tempR(j).eff2 = R(j).eff2(i);
            tempR(j).eff3 = R(j).eff3(i);
            tempR(j).eff4 = R(j).eff4(i);
        end
        cellvec{i} = tempR;
        
    end


end



