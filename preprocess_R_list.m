%This function preprocesses the R_list output from a simulation into a cell
%array of structs. Basically, it ensures consistent indexing even when
%GD-Calc was parameterized.

function result = preprocess_R_list(R_list) %#codegen

    %Arbitrary looping taking from:
    %https://www.mathworks.com/matlabcentral/answers/447469-looping-through-a-matrix-of-unknown-dimensions
    nv = ndims(R_list);
    v = [repmat({1}, 1, nv)];
    vLim = size(R_list);

    ready = false;
    while ~ready
        %Actual operation!
        %Check for existence
        if ~isa(R_list{v{:}}, 'double')
            lst = expandRstruct(R_list{v{:}});
            for i = 1:length(lst)
                R_list{v{:},i} = lst{i};
            end
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
    
    result = R_list;

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
        tempR = struct();
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