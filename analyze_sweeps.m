%Writing a script to analyze the height/spacing/period data we get from the
%simulation
%We have two outputs, R_list and T_list, cell arrays denoting the
%transmittion and reflection efficiencies respectively.
% Each cell array is indexed by the varied parameters
% For pyramids, this is height and grating period
% For cylinders, this is height, cylinder spacing, and cylinder radius.
%
% So, we'll have a toggle for the type of grating being used.
% Each cell in the array is a single "R (T)" struct, containing the
% efficiencies for the different diffraction orders. Since we did angle sweeps
% we need to extract that angle information...
%   Here, we see the effects of that tradeoff between performance and
%   usability. Parameterizing the angle information through GD-Calc makes
%   the functionality much better, but it also yields yucky data.
% The problem we have now is that we have 3-4 variables that we're
% changing, so we can't really plot things in a good way.
% In any case, we want to plot reflectance vs something
% Our options are angle, height, period, and (radius)
% 
%Again, the problem is that the reflectance is dependent on all of these,
%so how do we find a metric???
%I guess just choose the middle value for each and plot w/ that?

%A toggle for the type of structure we're analyzing
% 1 for pyramid
% 2 for cylinder
% 3 for cone
struc = 1;


%This extracts a single zero-order diffracted mode from the reflectivity
%list...
%tmp = R_list{1,1}(extractfield(R_list{1,1}, 'm1') == 0);
%tmp2 = tmp(extractfield(tmp, 'm2') == 0);


if struc == 1

    disp("Choices:")
    disp("1) Height")
    disp("2) Period")
    disp("3) Angle")
    xax = input("Select x-axis");
    lin = input("Select lines");

    %Now, we've defined which two variables we're dependent on...
    %For the other two, we'll take the first indexed position...
    %So, find the unselected choices:
    choices = [1,2,3];
    unselected = choices\[xax, lin];


    %Handling the angle data will probably be difficult...
    %BUT, the gist of it is that we have two values that we want to extract
    %all of (yielding a 2-D matrix), and one/two that we want to set
    %constant.
    %One approach that could be good is to eschew polarization information
    %altogether; if we look at reflectivity regardless of polarization
    %angle, it will simplify our data extraction process.
    %Obviously, this will not be a long-term solution, but it will work for
    %now.
    %What this will allow us to do is, for each cell in R_list, replace the
    %contained R struct with a list of the "averaged" reflectivity. 
    %The conversion for "natural light" is 0.5*(eff1 + eff2)

    [szx, szy] = size(R_list);
    Rl_no_pol = cell(szx, szy);
    
    for i = 1:szx
        for j = 1:szy
            Rl_no_pol{i,j} = (R_list{i,j}.eff1 + R_list{i,j}.eff2)./2;
        end
    end
    %So, now we have a 3-D "cube" structure indexed by {height, period}(angle)
    %We know that one of these will not be used, so we'll set that to 1...
    %Instead of using defined indices, we'll just use i1, i2, i3 for the
    %above...
    %If we want to select all, then we'll set it to 1; otherwise, it'll be
    %set to 0, and we just select the first index...
    
    i1 = 


%Rs = abs(((n1*cosd(theta)-n2*sqrt(1-((n1/n2)*sind(theta)).^2))./(n1*cosd(theta)+n2*sqrt(1-(n1/n2*sind(theta)).^2)))).^2;
%Rp = abs((n1*sqrt(1-(n1/n2 * sind(theta)).^2)-n2*cosd(theta))./(n1*sqrt(1-(n1/n2.*sind(theta)).^2)+n2.*cosd(theta))).^2;

    %We have several different cases...
    %We first choose the 
    switch lin



    end

    

end