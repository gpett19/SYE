%Simple function to convert a preprocessed R_list to a double matrix!

function Rmat = Rl_to_mat(R_list)

    Rmat = cell(size(R_list));
    for i = 1:length(R_list(:))
         if ~isa(R_list{i}, 'double')
                tmp = R_list{i}(extractfield(R_list{i}, 'm1') == 0);
                tmp2 = tmp(extractfield(tmp, 'm2') == 0);
                Rmat{i} = (tmp2.eff1 + tmp2.eff2)./2;
        end

    end
    Rmat = cell2mat(Rmat);

end
