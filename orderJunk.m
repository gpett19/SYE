%For each cell
F = {};
for i = 1:length(E)
   %Get the zero order thing.
   
   tmp = E{i};
   if(isempty(tmp))
        continue;
   end
   tmp2 = tmp(and([tmp.m2] == 0, [tmp.m1] == 0));

   F{i} = (tmp2.eff1 + tmp2.eff2 )/ 2;

end

F = cell2mat(F);

%Once we have the reflections, we just calculate the percent differences.
pDiff = [];

for i = 2:length(F)
    pDiff = [pDiff, 2*(F(i) - F(i-1))/(F(i) + F(i-1)) * 100];

end

