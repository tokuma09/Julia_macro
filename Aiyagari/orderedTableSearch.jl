#= This module is based on Numerical Recipes. =#
module orderedTableSearch
    @inline function locate(x, grid, btm=0, up=0)
        # result is same as searchsortedlast
        #grid must be monotonically increasing.

    if up == 0
        up = length(grid) + 1
    end

    N = length(grid)

    ascnd = (grid[end] > grid[1])

        # evaluation is out of range
    if x < grid[1]
        return 0
    elseif x > grid[end]
        return N
    end

    mid = 0

    @inbounds while up - btm > 1
        mid = div((up+btm),2)#floor(Int, (up + btm) / 2)
        if  ascnd == (x >= grid[mid])
            btm = mid
        else
            up = mid
        end
    end

    return btm
end

    @inline function hunt(x, grid, init_btm::Int)
    N = length(grid)
    ascnd = (grid[end] > grid[1])
    btm = init_btm
    up = 0

    if init_btm < 1 || btm > N
        return locate(x, grid)
    else
        inc = 1 # increment
        if (x > grid[btm]) == ascnd
            @inbounds while true
                up = btm + inc
                if up > N - 1
                    up = N
                    break
                elseif (x > grid[up]) == ascnd
                    btm = up
                    inc += inc
                else
                    break
                end
            end
        else
            up = btm
            @inbounds while true
                btm = up - inc
                if btm < 1
                    btm = 1
                    break
                elseif (x < grid[btm]) == ascnd
                    up = btm
                    inc += inc
                else
                    break
                end

            end
        end
    end
    return locate(x, grid, btm, up)
end
end
