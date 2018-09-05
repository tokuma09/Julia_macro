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