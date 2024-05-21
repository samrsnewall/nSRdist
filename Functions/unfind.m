function lgc = unfind(idx, N)
%This function returns a logical from a set of indices. It is like the
%reverse of the find function.

    %Go from indicies into logical (for vectors only)
    lgc = false(N,1);
    lgc(idx) = true;
end