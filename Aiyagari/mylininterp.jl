@inline function mylininterp(xpt, ypt, x)
    """
    This is just a wrapper of interpolations for convenience.
    This function do linear interpolation with extrapolation.

    Input
        xpt: grid
        ypt: vector that we want to interpolate
        x: evaluation point


    Return
        y: interpolated function
    """
    y = zeros(x)
    idx = trues(x)
  
    intf = interpolate((xpt,), ypt, Gridded(Linear()))
    y[idx] = intf[x[idx]]

    return y
end