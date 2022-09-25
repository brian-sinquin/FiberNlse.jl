
#! TODO better derivation
function derivate(y::AbstractVector, x::AbstractVector)
    function centraldiff(v::AbstractVector)
        dv = diff(v) / 2 # half the derivative
        a = [dv[1]; dv] # copies first element
        a .+= [dv; dv[end]] # copies last element, add both results to compute average
        return (a)
    end
    return centraldiff(y) ./ centraldiff(x)
end
