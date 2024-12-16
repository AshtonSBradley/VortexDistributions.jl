# Auxiliary Functions
function get_unit_vector(vec)
    vec_norm = norm(vec)
    return vec_norm == 0.0 ? vec : vec / vec_norm
end

# CCMA Class
struct CCMA
    w_ma::Int
    w_cc::Int
    distrib_ma::String
    distrib_cc::String
    rho_ma::Float64
    rho_cc::Float64
    weights_ma::Vector{Vector{Float64}}
    weights_cc::Vector{Vector{Float64}}
    w_ccma::Int

    function CCMA(; w_ma=5, w_cc=3, distrib="hanning", distrib_ma=nothing, distrib_cc=nothing, rho_ma=0.95, rho_cc=0.95)
        distrib_ma = distrib_ma === nothing ? distrib : distrib_ma
        distrib_cc = distrib_cc === nothing ? distrib : distrib_cc
        weights_ma = _get_weights(w_ma, distrib_ma, rho_ma)
        weights_cc = _get_weights(w_cc, distrib_cc, rho_cc)
        w_ccma = w_ma + w_cc + 1
        new(w_ma, w_cc, distrib_ma, distrib_cc, rho_ma, rho_cc, weights_ma, weights_cc, w_ccma)
    end
end

function _get_weights(w, distrib, rho)
    weight_list = []

   
    function get_hanning_kernel(window_size)
        window_size += 2
        hanning_kernel = (0.5 * (1 .- cos.(2 * π * (0:(window_size - 1)) / (window_size - 1))))[2:end-1]
        return hanning_kernel / sum(hanning_kernel)
    end

    for w_i in 0:w
        push!(weight_list, get_hanning_kernel(w_i * 2 + 1))
    end



    return weight_list
end

function _get_3d_from_2d(points)
    return hcat(points, zeros(size(points, 1)))
end

function _get_ma_points(points, weights)
    return hcat(conv(points[:, 1], weights)[length(weights):end - length(weights) + 1],
                conv(points[:, 2], weights)[length(weights):end - length(weights) + 1],
                conv(points[:, 3], weights)[length(weights):end - length(weights) + 1])
end

function _get_curvature_vectors(points)
    curvature_vectors = zeros(size(points))

    for i in 2:(size(points, 1) - 1)
        p0, p1, p2 = points[i - 1, :], points[i, :], points[i + 1, :]
        v1 = p1 - p0
        v2 = p2 - p1
        cross_p = cross(v1, v2)
        cross_norm = norm(cross_p)

        if cross_norm != 0.0
            radius = norm(v1) * norm(v2) * norm(p2 - p0) / (2 * cross_norm)
            curvature = 1.0 / radius
        else
            curvature = 0.0
        end

        curvature_vectors[i, :] = curvature * get_unit_vector(cross_p)
    end

    return curvature_vectors
end

function _get_alphas(points, curvatures)
    alphas = zeros(size(points, 1))

    for idx in 2:(size(points, 1) - 1)
        curvature = curvatures[idx]
        if curvature != 0.0
            radius = 1 / curvature
            dist_neighbors = norm(points[idx + 1, :] - points[idx - 1, :])
            a_sin_arg = (dist_neighbors / 2) / radius
            if a_sin_arg > 1.0
                a_sin_arg = 1.0
            end
            alphas[idx] = asin(a_sin_arg) ## Need to fix DomainError(1.0000000000000002, "asin(x) is not defined for |x| > 1.")
        else
            alphas[idx] = 0.0
        end
    end

    return alphas
end

function _get_normalized_ma_radii(alphas, w_ma, weights)
    radii_ma = zeros(length(alphas))

    for idx in 2:(length(alphas) - 1)
        radius = 1.0 * weights[w_ma + 1]

        for k in 1:w_ma
            radius += 2 * cos(alphas[idx] * k) * weights[w_ma + k + 1]
        end

        radii_ma[idx] = max(0.35, radius)
    end

    return radii_ma
end

function _get_descending_width(ccma::CCMA)
    descending_width_list = []
    w_ma_cur = ccma.w_ma
    w_cc_cur = ccma.w_cc

    while !(w_ma_cur == 0 && w_cc_cur == 0)
        if w_cc_cur >= w_ma_cur
            w_cc_cur -= 1
        else
            w_ma_cur -= 1
        end

        push!(descending_width_list, Dict("w_ma" => w_ma_cur, "w_cc" => w_cc_cur))
    end

    return descending_width_list
end

function _filter(ccma::CCMA, points, w_ma, w_cc, cc_mode)
    w_ccma = w_ma + w_cc + 1
    points_ma = _get_ma_points(points, ccma.weights_ma[w_ma + 1])

    if !cc_mode
        return points_ma
    end

    curvature_vectors = _get_curvature_vectors(points_ma)
    curvatures = mapslices(norm, curvature_vectors, dims=2)
    alphas = _get_alphas(points_ma, curvatures)
    radii_ma = _get_normalized_ma_radii(alphas, w_ma, ccma.weights_ma[w_ma + 1])
    points_ccma = zeros(size(points, 1) - 2 * w_ccma, 3)

    for idx in 1:(size(points, 1) - 2 * w_ccma)
        unit_tangent = get_unit_vector(points_ma[w_cc + idx + 2, :] - points_ma[w_cc + idx, :])
        shift = zeros(3)

        for idx_cc in 1:(2 * w_cc + 1)
            if curvatures[idx + idx_cc] == 0.0
                continue
            end

            u_vec = get_unit_vector(curvature_vectors[idx + idx_cc, :])
            weight = ccma.weights_cc[w_cc + 1][idx_cc]
            shift_magnitude = (1 / curvatures[idx + idx_cc]) * (1 / radii_ma[idx + idx_cc] - 1)
            shift += u_vec * weight * shift_magnitude
        end

        points_ccma[idx, :] = points_ma[idx + w_cc + 1, :] + cross(unit_tangent, shift)
    end

    return points_ccma
end

function filter(ccma::CCMA, points::Matrix{Float64}; mode="padding", cc_mode=true)
    if !(mode in ["none", "padding", "wrapping", "fill_boundary"])
        throw(ArgumentError("Invalid mode! Got :: $mode. Expected :: none | padding | wrapping | fill_boundary."))
    end

    if size(points, 1) < 3
        throw(ArgumentError("At least 3 points are necessary for the CCMA-filtering"))
    end

    if mode == "padding"
        n_padding = cc_mode ? ccma.w_ccma : ccma.w_ma
        points = vcat(repeat(points[1, :]', n_padding, 1), points, repeat(points[end, :]', n_padding, 1))
    end

    if size(points, 1) < ccma.w_ccma * 2 + 1
        throw(ArgumentError("Not enough points are given for complete filtering!"))
    end

    if mode == "wrapping"
        n_padding = cc_mode ? ccma.w_ccma : ccma.w_ma
        points = vcat(points[end-n_padding+1:end, :], points, points[1:n_padding, :])
    end

    is_2d = size(points, 2) == 2
    if is_2d
        points = _get_3d_from_2d(points)
    end

    if mode != "fill_boundary"
        if is_2d
            return _filter(ccma, points, ccma.w_ma, ccma.w_cc, cc_mode)[:, 1:2]
        end
        return _filter(ccma, points, ccma.w_ma, ccma.w_cc, cc_mode)
    else
        dim = is_2d ? 2 : 3

        if cc_mode
            points_ccma = zeros(size(points, 1), dim)
            descending_width_list = reverse(_get_descending_width(ccma))

            points_ccma[1, :] = points[1, 1:dim]
            points_ccma[end, :] = points[end, 1:dim]

            points_ccma[ccma.w_ccma+1:end-ccma.w_ccma, :] = _filter(ccma, points, ccma.w_ma, ccma.w_cc, cc_mode)[:, 1:dim]

            for width_set in descending_width_list
                w_ccma = width_set["w_ma"] + width_set["w_cc"] + 1
                use_ma_1 = width_set["w_ma"] == 0 && ccma.w_ma != 0

                points_ccma[width_set["w_ma"] + 1, :] = _filter(
                    ccma, points[1:width_set["w_ma"] + w_ccma + 1, :],
                    use_ma_1 ? 1 : width_set["w_ma"],
                    width_set["w_cc"],
                    !use_ma_1)[:, 1:dim]

                points_ccma[end-width_set["w_ma"], :] = _filter(
                    ccma, points[end-width_set["w_ma"]-w_ccma+1:end, :],
                    use_ma_1 ? 1 : width_set["w_ma"],
                    width_set["w_cc"],
                    !use_ma_1)[:, 1:dim]
            end

            return points_ccma
        else
            points_ma = zeros(size(points, 1), dim)
            descending_width_list = collect(0:ccma.w_ma-1)

            points_ma[ccma.w_ma+1:end-ccma.w_ma, :] = _filter(ccma, points, ccma.w_ma, 0, false)[:, 1:dim]

            for idx in descending_width_list
                points_ma[idx+1, :] = _filter(ccma, points[1:2*idx+1, :], idx, 0, false)[:, 1:dim]
                points_ma[end-idx, :] = _filter(ccma, points[end-2*idx:end, :], idx, 0, false)[:, 1:dim]
            end

            return points_ma
        end
    end
end