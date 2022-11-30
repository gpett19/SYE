%This function works to find the minima of a given preprocessed simulation
%output matrix.
%We first convert the cell array of structs into a matrix summarizing the
%essential information.
%For each struct, we'll first find the "unpolarized" scattered efficiency
%for each diffraction order.
%We'll then add those orders together in order to get a complete picture of
%the behaviour. This process could use some further analysis to see if it's
%viable or not, but we'll roll with it for now.
%Next, we find the minima, and consider using a gaussian blur-type
%function.

%{"Height", "\mu m", h, 0.7, 0.05}, {"Distance", "\mu m", d, 0.51, 0.05}, {"Radius", "%", 0.001:0.08:0.8, 0.48, 0.05}, {"Wavelength", "\mu m", wlr, 1.5, 0.05}, {"Angle", "Deg", theta, 0, 0}

%R_list is either a preprocessed or "raw" list containing the simulation
%data
%preprocesssed is a boolean flag denoting whether or not the R_list has
%been preprocessed for consistent indexing or not.
%k is the number of minima we want to analyze!
%varargin is the parameter names, ranges, etc... (same formatting as
%analyze_sweeps_2)
%Each parameter will be
%presented in a vector, of the form [paramName, unit, range,
%default_value, st_dev], e.g. ["Angle", "Degrees", 0:5:90, 50, 30].
function results = minimize_sweeps(R_list, preprocessed, kmins, varargin)
    arguments
        R_list
        preprocessed
        kmins
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


    %Preprocess the list if we need to
    if ~preprocessed
        R_list = preprocess_R_list(R_list);
    end

    %Now, we extract the requisite information.
    nv = ndims(R_list);
    v = [repmat({1}, 1, nv)];
    vLim = size(R_list);

    ready = false;
    while ~ready
        %Actual operation!
        %Check for existence
        if ~isa(R_list{v{:}}, 'double')
            s = R_list{v{:}};
            %Sum up all the diffraction orders for unpolarized light
            tmp = (extractfield(s, 'eff1') + extractfield(s, 'eff2')) / 2;
            R_list{v{:}} = sum(tmp);
        end

        %Update index vector
        ready = true;
        for k = 1:nv
            v{k} = v{k} + 1;
            if v{k} <= vLim(k)
                ready = false;
                break;
            end
            v{k} = 1;
        end
    end
    
    %Convert the resultant cell array into a double matrix
    R_list = cell2mat(R_list);
    
%-----------------------------------------------------------------------

    %Now, we want to find the k minima of the matrix, plus their indices...
    %This is easy enough when we *know* the dimensions of the array, but we
    %have to keep this extensible to arbitrary dimensions...

    [M, I] = mink(R_list(:), kmins);
    minIdx = cell(size(size(R_list)));
    [minIdx{:}] = ind2sub(size(R_list), I);

        
    %We also convert the minIdx cell array into a simple matrix (rows are
    %the different minima, columns are the indices.
    minIdx = cell2mat(minIdx);


    %For these minima, we want to then determine a score that takes into
    %account the nearby points.

    %Alright, we've got the indices for the k best dimensions in this cell
    %array now!
    
    %Now, we need to apply the gaussian blurring to each parameter, and
    %multiply the results...

    %First, compute the gaussian-weighted sum for each of the dimensions..
    %We'll want an "epsilon," past which we'll terminate the gaussian computation
    % in that direction...
    gEps = 0.0001;
    % Note that we'll have to deal with out-of-bounds errors...

    %So, I'm thinking for each dimension, we want to evaluate the gaussian
    %around that point with the specified SD. Since we're using the SD and
    %center relative to the actual sweep size, we can just plug in the
    %relevant point; if it's less than gEps, then we stop going in that
    %direction.


    %Get the indices we want:

    scores = zeros(1, length(minIdx(:,1)));
    
    %Loop through the indices!!
    for idxVar = 1:length(minIdx(:,1))
        idxs = minIdx(idxVar, :);
    %Now, we don't necessarily want arbitrary looping; rather, we want to
    %start at 'idxs', and continue within the hypercube defined by the SD
    %parameters...

    %So, the first order of business would be to find those boundaries!
    %Rather than doing comparisons with gEps for each sample, we can simply
    %evaluate the inverse of the gaussian we're looking at (special case I
    %believe).

    %However, we need to read in the SDs for each parameter first:
    sds = zeros(size(paramList));
    %Also have a cell array that will contain the vector bounds for each parameter
    bounds = cell(size(paramList));
    %Finally, we'll have a bounds array that contains the bounds for 2
    %standard deviations... this will be for warning out-of-bounds errors
    warnBounds = cell(size(paramList));
    for i = 1:length(paramList)
        sds(i) = paramList{i}{5};
        bounds{i} = inv_gauss_val(paramList{i}{3}(idxs(i)), sds(i), gEps);
        warnBounds{i} = inv_gauss_val(paramList{i}{3}(idxs(i)), sds(i), sds(i)*2);
    end

    %So, we now have the hypercube under which we're evaluating the
    %gaussian...

    %Necessarily, the exact points at which we evaluate the gaussian are
    %determined by the param ranges
    %My rudimentary knowledge of signal processing says we should take a
    %simple sum of the weighted values of all these scores...
    
    %This is the R_list that we will truncate from:
    R_list_trunc = R_list;

    %Cell array that will contain the gaussian values for each parameter
    gVals = cell(size(paramList));

    %We find the bound indices by looking at the ranges for each parameter
    %separately...
    for i = 1:length(paramList)
        rng = paramList{i}{3};
        %Find the logical array of indices that are greater than the lower bound
        lRng = rng > bounds{i}(1);
        %If this equals the trivial array, we know the lower bound is out
        %of bounds:
        if (lRng == (rng == rng))
            %So, check if we need to warn:
            lwRng = rng > warnBounds{i}(1);
            if (lwRng == lRng)
                disp(strcat("Warning: Edge of simulated range within 2 SD of point for parameter ", paramList{i}{1}, ", minimum ", string(idxVar), ". (Lower bound)"));
            end
        end

        %Now, do the upper bound:
        rRng = rng < bounds{i}(2);
        %Then, union the two; we store in rRng for stupid design reasons.
        rRng = and(lRng, rRng);
        %Check if union equals the lower bound array:
        if(rRng == lRng)
            %Check for warning bounds:
            rwRng = lRng < warnBounds{i}(2);
            if(rwRng == rRng)
                disp(strcat("Warning: Edge of simulated range within 2 SD of point for parameter ", paramList{i}{1}, ", minimum ", string(idxVar), ". (Upper bound)"));
            end
        end

        %No matter how we get here, we know that rRng is the logical array
        %we want to restrict that dimension to!
        tv = [repmat({':'}, size(paramList))];

        %Quick error check: if rRng is empty (i.e. we just want to look at
        %the center point, then just stick that index in there!
        if(~any(rRng))
            rRng(idxs(i)) = 1;
        end
        tv{i} = rRng;
        R_list_trunc = R_list_trunc(tv{:});

        %Now, recall we specified the center and SD in terms of the range
        %values; hence, we must find the value of the gaussian for each
        %true value in rRng.
        
        rVals = rng(rRng);
        for j = 1:length(rVals)
            rVals(j) = gaussian_value(paramList{i}{3}(idxs(i)), sds(i), rVals(j));
        end

        gVals{i} = rVals;
    end

    %Now, we should have a truncated R_list that contains only the values
    %we want to look at. Additionally, we should have a cell array that
    %contains the gaussian values for the relevant parameter sweeps (from
    %low to high param values) (e.g. gVals{1} = [0.5, 1, 0.5])
    
    %Now, for *every* value in the truncated R_list, we want to mutliply
    %said value by the relevant gaussian scaling values for the given
    %index. For example, if we're looking at R_list_trunc(1,1,1,1,1), we
    %want to calculate R_list_trunc(1,1,1,1,1)*gVals{1}(1)*gVals{2}(1)
    %I would like to cell2mat this, but I can't really do that due to the
    %fact that it wouldn't be rectangular...
    
    %So, arbitrary looping again!!!
    
    %

    nv = ndims(R_list_trunc);
    v = [repmat({1}, 1, nv)];
    vLim = size(R_list_trunc);

    R_list_score = zeros(size(R_list_trunc));

    ready = false;
    while ~ready
        %Actual operation!
        prod = 1;
        %We length by v instead of gVals because of the annoying special
        %case where the last gVal is a scalar 1.
        for i = 1:nv
            prod = prod*gVals{i}(v{i});
        end
        R_list_score(v{:}) = R_list_trunc(v{:})*prod;
        

        %Update index vector
        ready = true;
        for k = 1:nv
            v{k} = v{k} + 1;
            if v{k} <= vLim(k)
                ready = false;
                break;
            end
            v{k} = 1;
        end
    end
    
    scores(idxVar) = sum(R_list_score(:))/length(R_list_score(:));

    end




    results = {minIdx, scores, M, R_list};
    
end


%Returns the value of a 1-d gaussian function set to a height of 1, with mean
%given by center, standard deviation given by sd, and at the desired point.
function val = gaussian_value(center, sd, point)
    %If sd is 0, return 1 -> just get the center
    if(sd == 0)
        val = 1;
    else
        val = exp(-((point-center)^2/(2*sd^2)));
    end
end

%Returns the two x-values that have a "height" of value for a 1d gaussian
%function
function vals = inv_gauss_val(center, sd, value)
    tmp = sqrt(2*sd^2*log(1/value));
    vals(1) = center - tmp;
    vals(2) = center + tmp;
end
