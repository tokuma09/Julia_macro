function Aggregate_Labor(zgrid, Π)

    """
  Arguments:
  zgrid: Array{Float64,1}, shock grid
  Π: Array{Float64,2}, transition matrix

  Return:
  L: Float64, Aggregate Labor

  """

    # use eigenvector method to compute aggregate labor
    vals, vecs  = eig(Π')
    _, ind_val  = findmin(abs.(vals - 1.0))
    @views probst = vecs[:, ind_val] / sum(vecs[:, ind_val])
    L = dot(probst, zgrid)
    return L
end
