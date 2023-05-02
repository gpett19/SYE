%As described in my notebook entry for 2022-10-13, this is a function to
%plot arbitrary values for arbitrary R_lists.
%TODO: Pase that comment in here!

%Arguments explanation:
%R_list is the output from the simulation
%gd_param is a boolean value denoting whether or not the GD-Calc simulation
%was parameterized. If it was, set this to true.
%preprocessed is also a boolean value that describes if R_list was
%preprocessed or not... True means that we don't need to process it here.
%param1 is the parameter to be plotted on the x-axis as a string, same as
%that in paramList
%param2 is the parameter to be plotted as multiple lines with a legend.
%varargin takes all of the parameters swept over. Each parameter will be
%presented in a vector, of the form [paramName, unit, range,
%default_value], e.g. ["Angle", "Degrees", 0:5:90, 50].

%An example function call:
%analyze_sweeps_2(R_list, true, "Height", "Radius", {"Height", "um", 0.01:0.1:1, 0.5}, {"Distance", "um", 0.01:0.1:1, 0.91}, {"Radius", "um", 0.001:0.05:0.351, 0.51}, {"Angle", "Deg", 0:5:90, 0})

function analyze_sweeps_2(R_list, gd_param, preprocessed, param1, param2, varargin)

    arguments
        R_list cell
        gd_param logical
        preprocessed logical
        param1 string
        param2 string
    end
    arguments (Repeating)
        varargin
    end
    paramList = varargin;
    %Process the arguments! 

    %First, a sanity check; if the number of params specified does not
    %equal the number of dimensions in R_list, or that plus 1, we have an
    %issue!
    if length(paramList) ~= ndims(R_list) && length(paramList) ~= ndims(R_list) + 1
        error("Invalid number of parameters specified!")
    end


    %As detailed, the first thing we need to do is to preprocess R_list in
    %order to easily index into it.
    %Of course, this is only the case if we parameterized inside GD-Calc;
    %if not, we're already all set!
    if gd_param && ~preprocessed
        R_list = preprocess_R_list(R_list);
    end
    %Alright, we've cleaned up R_list now, so we can deal with all the
    %parameters with the same indexing!
    %Note that each R_list cell still is a massive struct with the efficiency orders
    %and efficiencies. However, if we can isolate the 2-D array of structs
    %we want *first*, we can then process that much more easily.
    %Now, we just need to process the input parameter data, and figure
    %out what exactly we're isolating

    %The first thing we need to do is deal with the parameter definitions,
    %and find which indices correspond to which parameters.

    p1idx = -1;
    p2idx = -1;

    for i = 1:length(paramList)
        if any(strcmp(paramList{i}{1}, param1))
            p1idx = i;
        elseif any(strcmp(paramList{i}{1}, param2))
            p2idx = i;
        end
        if p1idx ~= -1 && p2idx ~= -1
            break;
        end
    end

    %So, we know the indices which we're selecting *all* of; now, we need
    %to reduce the *other* indices to their defaults.

    idxs = ones(1,length(size(R_list)));
    idxs(p1idx) = 0;
    idxs(p2idx) = 0;

    %We do a similar thing to above, where we define a cell array of the
    %given length, populate it, and index into R_list with it!

    defaults = paramList(idxs == 1);
    index = cell(length(size(R_list)),1);
    index{p1idx} = ':';
    index{p2idx} = ':';

    %Wack as fuck but should work
    for i = 1:length(defaults)
        tmp = defaults{i};
        value = tmp{4};
        %Due to floating point errors, have to set tolerance value...
        %idx = find(tmp{3}==value);
        idx = find(abs(tmp{3} - value) < 0.0001);
        if isempty(idx)
            error("Invalid 'default' value")
        end
        %Map the default value to an index.
        defaults{i} = idx;
    end

    j = 0;
    for i = 1:length(paramList)
        if i ~= p1idx && i ~= p2idx
            index{i} = defaults{i + j};
        else
            j = j - 1;
        end

    end


    %Now that we've done this, we can extract the 2-D cell matrix of structs!
    
    Rmat = R_list(index{:});
    
    %For now, just deal with the natural zeroth diffracted order. However, this 
    % can (and should!) be changed later!
    
    [szx, szy] = size(Rmat);
    
    for i = 1:szx
        for j = 1:szy
            %Again, existence check...
            if ~isa(Rmat{i,j}, 'double')
                tmp = Rmat{i,j}(extractfield(Rmat{i,j}, 'm1') == 0);
                tmp2 = tmp(extractfield(tmp, 'm2') == 0);
                Rmat{i,j} = (tmp2.eff1 + tmp2.eff2)./2;
            end
        end
    end
    
    %disp("Here")

    Rmat = squeeze(cell2mat(Rmat));
    
    %Right, we should just have a plottable matrix now!
    %The last major issue we should have to solve is which dimension to plot on
    %the axis and which to plot on the lines
    %According to the plot documentation, the dimensions are first matched; if
    %the matrix is square, each column is a separate line!
    %Just to handle the square matrix case:
    if szx == szy
        %p1idx is the columns; this is what we don't, so we flip!
        if p1idx > p2idx
            Rmat = Rmat';
        end
    end

    %Do a bit of filtering; if the first parameter is angle, we want the
    %option to normalize it to the fresnel reflection!

    %Boolean on whether to normalize angle or not
    normAngle = 1;

    if(strcmp(paramList{p1idx}{1}, "Angle") && normAngle == 1)
        %Compute the fresnel reflection!
        %FIXME: This currently is hardcoded, because we don't usually save
        %this information. REMEMBER TO CHANGE if needed!
        n1 = 1;
        n2 = 1.444; %Fused silica
        theta = paramList{p1idx}{3};
        Rs = abs(((n1*cosd(theta)-n2*sqrt(1-((n1/n2)*sind(theta)).^2))./(n1*cosd(theta)+n2*sqrt(1-(n1/n2*sind(theta)).^2)))).^2;
        Rp = abs((n1*sqrt(1-(n1/n2 * sind(theta)).^2)-n2*cosd(theta))./(n1*sqrt(1-(n1/n2.*sind(theta)).^2)+n2.*cosd(theta))).^2;
        Rtot = (Rs + Rp) ./ 2;

        Rmat = Rmat./ Rtot;
    end
    
    close all;
    plot(paramList{p1idx}{3}, Rmat)
    colororder(turbo(length(paramList{p2idx}{3})))
    xlabel(strcat(paramList{p1idx}{1}, "(", paramList{p1idx}{2}, ")"));
    ylabel("Reflectance");
    %Include other parameters in the title:
    titleStr = "";
    for i = 1:length(paramList)
        if i == p1idx || i == p2idx
            %This is an actual param, so don't include it in the title.
            continue;
        else
            %Deal with percentaged values....
            %I.e. a radius value of 1 means that 100% of available radius
            %space was filled.
            if(strcmp(paramList{i}{2}, "%"))
                paramList{i}{4} = paramList{i}{4}*100;
            end


            titleStr = strcat(titleStr, paramList{i}{1}, "=", string(paramList{i}{4}), paramList{i}{2}, ", ");
        end
    end
    %Remove trailing comma:
    if titleStr ~= ""
        titleStr = extractBefore(titleStr, strlength(titleStr)-1);
    end

    title(strcat("Reflectance vs. ", paramList{p1idx}{1}, " (", titleStr ,")"))
    lgd = legend(string(paramList{p2idx}{3}));
    title(lgd, strcat(paramList{p2idx}{1}, "(", paramList{p2idx}{2},")"))
    
    %Change the labels if normalized
    if(strcmp(paramList{p1idx}{1}, "Angle") && normAngle == 1)
        title(strcat("Normalized Reflectance vs. ", paramList{p1idx}{1}, " (", titleStr ,")"))
        ylabel("Normalized Reflectance");
    end
    
end

