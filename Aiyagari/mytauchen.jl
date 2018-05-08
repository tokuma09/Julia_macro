function mytauchen(nz,ρ,σ,s::Int=3)
    """
    need package Distributions
    Discretize the following AR(1) process.
    z_t = ρz_{t-1} + ϵ_t

    Arguments:
    nz::Int64 :grid size # step 1
    ρ::Float64 :persistence of AR(1) process
    σ::Float64 :standard diviation of random components
    s::Int64 :decide maximum and minimum value of node often used s=3.

    return:
    z::array{float64,1} :discretized grid 
    Π::array{float64,2} :transition matrix of z
    """
    d = Normal() # standard normal distribution
    
    ###########################################
    # step 2. decide max and min value of grid
    ##########################################
    
    σ_z = sqrt(σ^2/(1.0-ρ^2))
    zmin = -s*σ_z
    zmax = s*σ_z
    z = collect(linspace(zmin,zmax,nz))
    
    ###########################################
    # step 3. store width of between grid points
    ###########################################
    
    w = z[2] -z[1] #width of between grid points

    ###########################################
    # step 4. Get transition matrix
    ###########################################
    
    Π = zeros(nz,nz) # transition matrix

    for row in 1:nz
        # Do end points first
        @inbounds Π[row, 1] = cdf(d,(z[1] - ρ*z[row] + w/2) / σ)

        @inbounds Π[row, nz] = ccdf(d,(z[nz] - ρ*z[row] - w/2) / σ) 
                               #1.0- cdf((z[nz] - ρ*z[row] - w/2) / σ)
    end

    # fill in the middle columns
    for row in 1:nz
        for col in 2:nz-1
            @inbounds Π[row, col] = (cdf(d,(z[col] - ρ*z[row] + w/2) / σ) 
                                      -cdf(d,(z[col] - ρ*z[row] - w/2) / σ))
        end
    end

    return z,Π
end