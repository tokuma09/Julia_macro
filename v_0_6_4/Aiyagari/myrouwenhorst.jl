function myrouwenhorst(N::Integer,ρ::Real,σ::Real)
    
    """
    Consider the follwoing process
    z_{t+1} = ρz_t  + ϵ_{t+1}

    Arguments:
    N::Int64 :grid size 
    ρ::Float64 :persistence of AR(1) process
    σ::Float64 :standard diviation of random components
    
    Return:
    z::array{float64,1} :discretized grid 
    Π::array{float64,2} :transition matrix 
    
    """
    
    
    #####################################################
    # step 2. create state grid
    #####################################################
    σ_z = σ/sqrt(1.0-ρ^2)
    ψ = sqrt(N-1)* σ_z # end points
    # create shock grid
    zgrid = collect(linspace(-ψ,ψ,N))
    
    #####################################################
    # step 3. compute trandition matrix
    #####################################################
    
    p = (1.0+ρ)/2.0
    Π   =  [p 1-p;1-p p]
    @inbounds for i in 3:N
        zero_v = zeros(i-1,1)
        zero_v_long = zeros(1,i)
        Π = p.*[Π zero_v; zero_v_long] +
             (1-p).*[zero_v Π; zero_v_long] +
             (1-p).*[zero_v_long; Π zero_v] +
             p.*[zero_v_long; zero_v Π]
    ######################################################
    # step 4. devide matrix by 2 except the top and bottom
    ######################################################
        @views Π[2:end-1,:] ./=  2.0 
    end
            
    return zgrid, Π
end