function x = remmin(x)
for i = 1:numel(x)
    if ~isempty(symvar(x(i)))
         x(i) = -1e100;
    end
end

x = sort(x);
x= x(end);

  
