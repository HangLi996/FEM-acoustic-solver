function dim = getDimension(eltype)
%   Input:  eltype  ... type of the respective element(s)
%   Output: dim     ... dimension of the domain
%
%   Get dimension of the domain.
%
%   Written by Lennart Moheit (Lennart.Moheit@tum.de)
%   11/08/2014
%
    if any(ismember(eltype, [4 5 6 7 11 12 13 14 17 18 19])) % any 3D elements?
        dim = 3;
        return
    end
    if any(ismember(eltype, [2 3 9 10 16 21])) % any 2D elements?
        dim = 2;
        return
    end
    if any(ismember(eltype, [1 8])) % any 1D elements?
        dim = 1;
    elseif any(ismember(eltype, [15])) % at least a simple dot? Come on ;-)
        dim = 0;
    else % Nope, nothing!
        disp('! WARNING: Invalid element type.')
    end
end