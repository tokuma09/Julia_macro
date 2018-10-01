###########################################################
# This is my module that contains some wrapper functions.
# Author: Tokuma Suzuki
###########################################################


using LinearAlgebra, SparseArrays

function speye(size, dtype=Float64)
    """
    Return sparse identity matrix.

    This function needs LinearAlgebra and SparseArrays.
    Arguments:
    size: Int64
    dtype: keywords, Float64 or Int64
    """

    return sparse(Matrix{dtype}(I, size, size))
end


function linspace(start,stop,size=50)
    """

    Return linealy spaced 1d array

    Arguments:
    start: Real
    end: Real
    size: Int
    """

    return collect(LinRange(start,stop,size))
end

function interp1(xpt, ypt, x; method="linear", extrapvalue=nothing)

    """
    This function needs Interpolations and works like interp1 in MATLAB.
    If you want to input array{Float64,1} as xpt, choose "linear"

    Input
        xpt: grid
        ypt: vector that we want to interpolate
        x: evaluation point
        method: "linear" or "cubic"
        extravalue: default is nothing

    Return
        y: interpolated function
    """

    if extrapvalue == nothing
        y = zeros(x)
        idx = trues(x)
    else
        y = extrapvalue*ones(x)
        idx = (x .>= xpt[1]) .& (x .<= xpt[end])
    end

    if method == "linear"
        intf = interpolate((xpt,), ypt, Gridded(Linear()))
        y[idx] = intf[x[idx]]

    elseif method == "cubic"
        itp = interpolate(ypt, BSpline(Cubic(Natural())), OnGrid())
        intf = Interpolations.scale(itp, xpt)
        y[idx] = [intf[xi] for xi in x[idx]]
    end

    return y
end
