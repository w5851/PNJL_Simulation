module init

using FastGaussQuadrature
export gauleg
# 1. gauleg 保持不变
function gauleg(a, b, n)
    t_nodes, t_weights = gausslegendre(n)
    nodes   = @. (b - a)/2 * t_nodes + (a + b)/2
    weights = @. (b - a)/2 * t_weights
    return nodes, weights
end


end  # module init
