function tauchen(nz,ρ,σ,s::Int=3)
    """
    need package Distributions
    consider log(z_t) = ρlog(z_{t-1}) + ε
    define x_t = log(z_t) and we consider
    this version.

    """

    """
    Arguments:
    nz: grid size
    ρ: persistence of AR(1) process
    σ: standard diviation of random components
    s: decide maximum and minimum value of node

    return:
    z: discretized grid array{float64,1}
    Π: transition matrix of x array{float64,2}
    """
    d = Normal() # standard normal distribution
    σ_x = sqrt(σ^2/(1.0-ρ^2))
    xmin = -s*σ_x
    xmax = s*σ_x
    x = linspace(xmin,xmax,nz)
    w = x[2] -x[1] #width of each grid

    # Get transition matrix

    Π = zeros(nz,nz)

    for row in 1:nz

        # Do end points first

        @inbounds Π[row, 1] = cdf(d,(x[1] - ρ*x[row] + w/2) / σ)

        @inbounds Π[row, nz] = ccdf(d,(x[nz] - ρ*x[row] - w/2) / σ) #1.0- cdf((x[N] - ρ*x[row] - w/2) / σ)
    end

    # fill in the middle columns
    for row in 1:nz
        for col in 2:nz-1

            @inbounds Π[row, col] = (cdf(d,(x[col] - ρ*x[row] + w/2) / σ) -

                           cdf(d,(x[col] - ρ*x[row] - w/2) / σ))
        end
    end

    #we consider x=log(z_t) case. so I want to return z_t
    z = exp.(x)
    return z,Π
end
