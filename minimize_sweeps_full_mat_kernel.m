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

% This function actually computes a gaussian kernel, and then convolves
% that with the data matrix.
% ---- NOTE: IMPORTANT: The fact that we're using a kernel means that the
% spacing for a given parameter must be CONSTANT! For example, angles of
% [0, 10, 30, 70] will give junk data, since the kernel will be computed
% for one of these spacings, but will be convolved with the entire data
% matrix.


%Example of a set of parameters for varargin
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
function results = minimize_sweeps_full_mat_kernel(R_list, preprocessed, kmins, varargin)
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
    
    %TODO: The following is catastrophically stupid, and should be fixed ASAP.
    % DON'T DO THIS!!!

    %Get the axis along which the angle measurements are stored.
    angIdx = -1;

    for i1 = 1:length(paramList)
        if any(strcmp(paramList{i1}{1}, "Angle"))
            angIdx = i1;
        end
        if angIdx ~= -1 
            break;
        end
    end

    %Cell array that contains the indices...
    aVec = [repmat({':'}, 1, nv)];
    aVec{angIdx} = 2:vLim(angIdx);
    R_list_doubled = flip(R_list(aVec{:}), angIdx);
    R_list_doubled = cat(angIdx, R_list_doubled, R_list);

    clear  v vLim

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

    %So, I'm thinking for each dimension, we want to evaluate the gaussian
    %around that point with the specified SD. Since we're using the SD and
    %center relative to the actual sweep size, we can just plug in the
    %relevant point; if it's less than gEps, then we stop going in that
    %direction.

    gKernel = compute_kernel(paramList);
    
    %Get the indices we want:

    scores = convn(R_list_doubled, gKernel, 'valid');


%    scores = zeros(size(R_list));
%     nvo = ndims(R_list);
%     vo = [repmat({1}, 1, nvo)];
%     vLimo = size(R_list);
% 
%     readyo = false;
%     while ~readyo
%         %Actual operation!
%        
%             score = gaussify(vo);
% 
%             scores(vo{:}) = score;
% 
%         %Update index vector
%         readyo = true;
%         for k = 1:nvo
%             vo{k} = vo{k} + 1;
%             if vo{k} <= vLimo(k)
%                 readyo = false;
%                 break;
%             end
%             vo{k} = 1;
%         end
%     end

    %scores = flip(scores, angIdx);

    %Since we convolved with the "doubled" matrix, we have to truncate the
    %scores matrix, as well as divide by the number of entries in the
    %kernel to have comparable results to the previous gaussify function.
    %Truncates the scores matrix along the angle dimension
    scrSz = size(scores);
    scrMidpt = ceil(scrSz(angIdx) / 2);
    aVec = [repmat({':'}, 1, nv)];
    aVec{angIdx} = scrMidpt:scrSz(angIdx);
    scores = scores(aVec{:});
    %Then, divides by kernel size:
    scores = scores ./ length(gKernel(:));

    %Pad the scores array with NaNs to have consistent indexing with the
    %R_list array. 
    %First, we need to calculate the padding needed for each dimension.
    scrSz = size(scores);
    kerSz = size(gKernel);
    %This is the amount truncated on each end of each dimension by the
    %convolution
    padSz = floor(kerSz ./ 2);

    %The amount we need to add; generally, this would be padSz, BUT recall
    %that we DON'T want to pad the angles (because we kept the first bits)
    %Note that size(R_list) is just scrSz + 2*padSz for each dimension
    %(EXCEPT ANGLE). Ultimately, we want the size of the padded array to be
    %size(R_list) - padSz for each dimension.
    padAmt = size(R_list) - scrSz - padSz;

    scores = padarray(scores, padAmt, NaN, 'pre');


    
    
    [Ms, Is] = mink(scores(:), kmins);
    minIdxS = cell(size(size(scores)));
    [minIdxS{:}] = ind2sub(size(scores), Is);
    minIdxS = cell2mat(minIdxS);



    results = {minIdx, M, minIdxS, Ms, scores, R_list};
    
%
%
% - - - - - - -- - - - -- - - -- - - -- - - -- - - -- - - -- - - - -
%
%
    function score = gaussify(vo)


    %Compute the gaussiam matrix:
    %To do this, we need the standard deviation, *plus* the points at which
    %we actually want to compute!
    %I.e. in a typical matrix the spacing between points is 1, so we can
    %evaluate the gaussian one to the left or right of center. But, here,
    %the spacing is determined by the dimension we're looking at!




    %-----------------------------------------------------
    %Loop through the indices!!
    %for idxVar = 1:length(minIdx(:,1))
        %idxs = minIdx(idxVar, :);
        idxs = [vo{:}];
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
    %Finally, we'll have a bounds array that contains the bounds for 3
    %standard deviations... this will be for warning out-of-bounds errors
    warnBounds = cell(size(paramList));
    for i = 1:length(paramList)
        sds(i) = paramList{i}{5};
        bounds{i} = inv_gauss_val(paramList{i}{3}(idxs(i)), sds(i), gEps);
        %warnBounds{i} = bounds{i};
        warnBounds{i} = inv_gauss_val(paramList{i}{3}(idxs(i)), sds(i), sds(i)*3);
    end

    %So, we now have the hypercube under which we're evaluating the
    %gaussian...

    %Necessarily, the exact points at which we evaluate the gaussian are
    %determined by the param ranges
    %My rudimentary knowledge of signal processing says we should take a
    %simple sum of the weighted values of all these scores...
    
    %This is the R_list that we will truncate from:
    R_list_trunc = R_list_doubled;

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
                %Instead of displaying, we don't want to include this; just
                %return!
                %disp(strcat("Warning: Edge of simulated range within 2 SD of point for parameter ", paramList{i}{1}, ", minimum ", string(idxVar), ". (Lower bound)"));
                %Now, we'll just check if this is one of the minima; if it
                %is, just pass a warning!
                if(ismember(idxs, minIdx, 'rows'))
                    disp(strcat("Warning: Minimum is within lower truncating range for indices ", strjoin(string(idxs)), ". Consider widening simulation parameters."))
                end
                score = NaN;
                return;
            end
        end

        %Now, do the upper bound:
        rRng = rng < bounds{i}(2);
        %Then, union the two; we store in rRng for stupid design reasons.
        rRng = and(lRng, rRng);
        %Check if union equals the lower bound array:
        if(rRng == lRng)
            %Check for warning bounds:
            rwRng = rng < warnBounds{i}(2);
            rwRng = and(lRng, rwRng);
            if(rwRng == rRng)
                %Similarly, just return!
                %disp(strcat("Warning: Edge of simulated range within 2 SD of point for parameter ", paramList{i}{1}, ", minimum ", string(idxVar), ". (Upper bound)"));
                if(ismember(idxs, minIdx, 'rows'))
                    disp(strcat("Warning: Minimum is within upper truncating range for indices ", strjoin(string(idxs)), ". Consider widening simulation parameters."))
                end
                
                score = NaN;
                return;
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
        %TODO: Remove when refactoring angles

        %if(~(i == angIdx))
        if (true)
            tv{i} = rRng;
        else
            rRng = cat(2,false(1,length(rRng)-1),rRng);
            tv{i} = rRng;
        end
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
    
    %end
    %----------------------------------

    score = sum(R_list_score(:))/length(R_list_score(:));

    end

    
end



%This function computes the gaussian kernel, given the standard deviations
%and spacings provided.
function ker = compute_kernel(paramList)
    %First order of business: read in the SDs and spacings:
    sds = zeros(1,length(paramList));
    spacings = zeros(1,length(paramList));
    for i = 1:length(paramList)
        %Get the user-provided SDs.
        sds(i) = paramList{i}{5};
        %Get the spacings! This will impose the restriction of constant
        %spacing
        %Check to make sure we have at least two
        if(length(paramList{i}{3}) <= 1)
            spacings(i) = 0;
        else
            spacings(i) = paramList{i}{3}(2) - paramList{i}{3}(1);
        end
    end
    %So, we now have spacings and standard deviations. For each dimension
    %now, we need to compute the normalized gaussian, and then compile all
    %of those into a matrix.
    %So, do the gaussians first:
    
    gaussians = cell(1,length(paramList));
    %For each parameter...
    for i = 1:length(paramList)
        %Go in the +'ve direction until hit three sigma:
        j = 0;
        jIdx = 0;
        %Preallocate the array for gaussian values.
        %Basically, we see how many spacings can fit in one half of the
        %gaussian, double it, and add one for the center point
        vals = zeros(1,floor(3*sds(i)/spacings(i))*2 + 1);
        %Get the midpoint:
        mdpt = ceil(length(vals) / 2);
        while j <= 3*sds(i)
            vals(mdpt + jIdx) = gaussian_value(0, sds(i), j);
            %Copy it over to the other side
            vals(mdpt - jIdx) = vals(mdpt + jIdx);
            j = j + spacings(i);
            jIdx = jIdx + 1;

        end
        %Should now have a vector of gaussain values for that parameter.
        gaussians{i} = vals;
    end

    %Right, now we have n different gaussian values; we just need to
    %compute the matrix now!
    %This is a bit tricky, but we can do this algorithmically, building up
    %from one vector.
    %Say we have an N-1 dimensional array, and want to add the gaussian
    %values for vector N (i.e. we have [0.5, 1, 0.5 ; 1, 2, 1; 0.5, 1,
    %0.5] and want to add a third dimension). Basically, we consider the
    %first matrix as a "frame"; for each value "i" in the vector, we multiply
    %that value by the entire first matrix, and the resultant is the "ith"
    %value for the Nth dimension in the new N-dimensional matrix.
    %Var for the kernel matrix:
    kMat = gaussians{1};
    for i = 2:length(gaussians)
        gauss = gaussians{i};
        %Wacky cell indexing:
        kmtidx = [repmat({':'}, 1, i - 1)];
        %For each value in the vector
        for j = 1:length(gauss)
            kMT(kmtidx{:}, j) = kMat*gauss(j);
        end
        kMat = kMT;
    end

ker = kMat;

end

%Returns the value of a 1-d gaussian function set to a height of 1, with mean
%given by center, standard deviation given by sd, and at the desired point.
function val = gaussian_value(center, sd, point)
    %If sd is 0, return 1 -> just get the center
    if(sd == 0)
        val = 1;
    else
        val = (1/sqrt(2*pi)*sd)*exp(-((point-center)^2/(2*sd^2)));
    end
end

%Returns the two x-values that have a "height" of value for a 1d gaussian
%function
function vals = inv_gauss_val(center, sd, value)
    tmp = sqrt(2*sd^2*log(1/value));
    vals(1) = center - tmp;
    vals(2) = center + tmp;
end
