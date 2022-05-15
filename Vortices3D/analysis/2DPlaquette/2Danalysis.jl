function errorTest(num_vorts, n, version)
    Lx = 200; Ly = Lx;
    Nx = n; Ny = n;
    # x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1]; ##periodic 
    x = LinRange(-Lx/2,Lx/2, Nx); y = LinRange(-Ly/2,Ly/2, Ny);

    psi0 = one.(x*y') |> complex
    psi = Torus(psi0,x,y)
    rand_vorts = rand_pointvortex(num_vorts, psi)
    vortex!(psi, rand_vorts)

    vf = findvortices(psi)
    if version == 1
        return findMinAvgError(psi, rand_vorts, vf)/(x[2]-x[1])
    elseif version == 2
        return findMinAvgError2(psi, rand_vorts, vf)/(x[2]-x[1])
    else
        return [findMinAvgError(psi, rand_vorts, vf)/(x[2]-x[1]), findMinAvgError2(psi, rand_vorts, vf)/(x[2]-x[1])]
    end
end

function findMinAvgError2(psi, vorts, vfound)
    @unpack ψ, x, y = psi;
    # @assert length(vorts) == length(vfound)
    # @assert length(vorts) == length(vfound)
    Lx = x[end]-x[1]; Ly = y[end]-y[1];

    if length(vorts) == length(vfound)

        vcopy = copy(vorts)
        vfcopy = copy(vfound)
        diff = []

        while length(vcopy) > 0
            v_curr = pop!(vcopy)
            min = Lx+Ly; min_idx = 0;
            for i in 1:length(vfcopy)
                curr_min = euclidVorts(v_curr, vfcopy[i])
                if curr_min < min
                    min = curr_min; min_idx = i;
                end
            end
            @assert min_idx != 0;
            push!(diff, min)
            deleteat!(vfcopy, min_idx)
        end
        return mean(diff)
    else
        print(Lx+Ly)
        return Lx+Ly

    end
end

function findMinAvgError(psi, vorts, vfound)
    @assert length(vorts) == length(vfound)
    @unpack ψ,x,y = psi
    Lx = x[end]-x[1]; Ly = y[end]-y[1];

    ## vorts and vfound are N element vectors
    p = collect(permutations(vfound))
    error = Lx+Ly
    for i in 1:length(p)
        temp = [euclidVorts(vorts[j], p[i][j]) for j in 1:length(vorts)]
        temp_error = mean(temp)
        if temp_error < error
            error = temp_error
        end
    end
    return error
end

function benchmark2Dsuite(maxSec, num_vorts, Nstart, Nend)
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = maxSec

    suite = BenchmarkGroup();
    
    suite["naive"] = BenchmarkGroup();
    suite["optim"] = BenchmarkGroup();
    suite["optimInterp"] = BenchmarkGroup();

    N = [2^x for x=Nstart:Nend];
    for n in N
        psi = setupTestArea(num_vorts, n)
        suite["naive"][string(n)] = @benchmarkable naivePlaquette($psi)
        suite["optim"][string(n)] = @benchmarkable remove_vortices_edge(findvortices_jumps($psi), $psi)
        suite["optimInterp"][string(n)] = @benchmarkable findvortices($psi)
    end
    return suite
end

function benchmark2D(num_vorts, Nstart, Nend)
    N = [2^x for x=Nstart:Nend]
    benchmarks = []
    for n in N
        println("Running tests for: " * string(n))
        psi = setupTestArea(num_vorts, n)
        benchNaive = @benchmark naivePlaquette($psi)
        benchOptimised = @benchmark remove_vortices_edge(findvortices_jumps($psi), $psi)
        benchOptimisedInterp = @benchmark findvortices($psi)
        push!(benchmarks, [benchNaive, benchOptimised, benchOptimisedInterp, n])
    end
    return benchmarks
end

function setupTestArea(num_vorts, n)
    Lx = 200; Ly = Lx;
    Nx = n; Ny = n;
    # x = LinRange(-Lx/2,Lx/2, Nx+1)[1:end-1]; y = LinRange(-Ly/2,Ly/2, Ny+1)[1:end-1]; ##periodic 
    x = LinRange(-Lx/2,Lx/2, Nx); y = LinRange(-Ly/2,Ly/2, Ny);

    psi0 = one.(x*y') |> complex
    psi = Torus(psi0,x,y)
    rand_vorts = rand_pointvortex(num_vorts, psi)
    vortex!(psi, rand_vorts)
    return psi
end


function euclid(x, y)
    return sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2)
end

function euclidVorts(v1, v2)
    return euclid([v1.xv, v1.yv], [v2.xv, v2.yv])
end

function naivePlaquetteInterp(psi, X, ϵ)
    x = X[1]; y = X[2];
    dx = x[2]-x[1]; dy=y[2]-y[1];
    N = 1;

    x_itp = interpolate(x, BSpline(Linear()));
    y_itp = interpolate(y, BSpline(Linear()));
    psi_itp = interpolate(psi, BSpline(Quadratic(Periodic(OnCell()))))
    x_range = LinRange(1,length(x),N*(length(x)))
    y_range = LinRange(1,length(y),N*(length(y)))
    # psi_itp = interpolate((x, y), psi, Gridded(Linear()))
    
    psi_itp_arr = psi_itp[x_range, y_range]

    vorts = []
    phase = angle.(psi_itp_arr)
    for i in 1:length(phase[:, 1])-1
        for j in 1:length(phase[1, :])-1
            # start from top left corner and work anti-clockwise
            square = [phase[i, j], phase[i+1, j], phase[i+1, j+1], phase[i, j+1], phase[i, j]]
            for k in 2:5
                # diff = square[k]-square[k-1]
                while (diff = square[k]-square[k-1]) > π
                    square[k] += -2*pi
                end
                while (diff = square[k]-square[k-1]) < -π
                    square[k] += 2*pi
                end
            end
            phase_diff = square[end] - square[1]
            m = phase_diff/(2*π)
            if abs(m) > ϵ
                ## vort in box
                push!(vorts, PointVortex(x_itp[x_range[i]] + (dx/(2*N)), y_itp(y_range[j])+(dy/(2*N)), round(m)))
            end
        end
    end
    return vorts
end

function naivePlaquette(psi)
    @unpack ψ,x,y = psi
    dx = x[2]-x[1]; dy=y[2]-y[1];

    vorts = []
    phase = angle.(ψ)
    for i in 1:length(x)-1
        for j in 1:length(y)-1
            # start from top left corner and work anti-clockwise
            square = [phase[i, j], phase[i+1, j], phase[i+1, j+1], phase[i, j+1], phase[i, j]]
            for k in 2:5
                # diff = square[k]-square[k-1]
                while (diff = square[k]-square[k-1]) > π
                    square[k] += -2*pi
                end
                while (diff = square[k]-square[k-1]) < -π
                    square[k] += 2*pi
                end
            end
            phase_diff = square[end] - square[1]
            m = phase_diff/(2*π)
            if abs(m) > 0
                ## vort in box
                push!(vorts, PointVortex(x[i] + dx/2, y[j]+dy/2, round(m)))
            end
        end
    end
    return vorts
end

function benchmarkDict(result)
    newResult = Dict{String,Dict{String, Any}}()
    typeDict = Dict([("min", minimum(result)), ("mean", mean(result)), ("median", median(result))])
    
    
    algos = collect(keys(result))
    types = ["min", "mean", "median", "std"]
    for alg in algos
        typeDict = Dict{String, Any}([("min", minimum(result[alg])), ("mean", mean(result[alg])), ("median", median(result[alg]))])
        for type in types[1:3]
            temp = leaves(typeDict[type])
            temp = [(parse.(Int, temp[i][1][1]), temp[i][2].time) for i in 1:length(temp)]
            temp = hcat(collect.(temp)...)'
            temp = sortslices(temp,dims=1,by=x->x[1],rev=false)
            typeDict[type] = temp
        end
        ## std
        std_keys = collect(keys(result[alg]))
        temp = [(parse.(Int, key), std(result[alg][key]).time) for key in std_keys]
        temp = hcat(collect.(temp)...)'
        temp = sortslices(temp,dims=1,by=x->x[1],rev=false)
        typeDict["std"] = temp
        newResult[alg] = typeDict
    end
    return newResult
end
