using JLD2
using VortexDistributions
using GLMakie

## Load data
@load "src/VortexDetection3D/sim_examples/box_rings.jld2" psi_ring1 psi_ring2 psi_ring3 X

## Plot iso surface for select ring 
function plot_iso(psi)
    density = abs2.(psi)
    pmax = maximum(density)
    density = density/pmax
    volume(density, algorithm = :iso, show_axis=true)
end
plot_iso(psi_ring3)

psi = psi_ring3 # Use "psi" as current sim you're working on
grad = gradient_3D_cent(psi, X)
grad_i = grad_itp(grad, X)
x = X[1]; y = X[2]; z = X[3]

dz = z[2] - z[1]

vorts_3d = find_vortices3D_v2(psi, X)
label = 1;
for zidx in 2:length(z)-2
    vorts_slice = vorts_3d[zidx];
    num_vorts = length(vorts_slice) # Number of vorts on z-slice
    if num_vorts != 0
        num_vorts = length(vorts_slice[:, 1])
        for v in 1:num_vorts
            if (vorts_slice[v, 5] == 0) # Check hasn't been labeled
                vorts_3d[zidx][v, 5] = label;
                label += 1;
            end
            
            vorts_sliceup = vorts_3d[zidx+1];
            
            if length(vorts_sliceup) != 0
                wps = wps_int(grad_i, [vorts_slice[v, 1], vorts_slice[v, 2], z[zidx]]);
                wps = wps .* (dz/abs(wps[3]));
                
                if vorts_slice[v, 4] > 0
                    approx_vort = vorts_slice[v, 1:2] + wps[1:2]
                else 
                    approx_vort = vorts_slice[v, 1:2] - wps[1:2]
                end
                vort_index = find_closest_tuple(vorts_sliceup[:, 1:2], [approx_vort[1], approx_vort[2]]);
                if vorts_sliceup[vort_index, 4] == vorts_slice[v, 4]
                    vorts_sliceup[vort_index, 5] = vorts_slice[v, 5] # Assign same label
                end

            end

        end
    end
end
vorts_3d

# To see which layers have vortices found 
# for i in 1:64
#     if length(vorts_3d[i]) > 0
#         println(i)
#     end
# end


vorts_label = Array{Float64, 2}[];
for i in 1:label-1
    current_vort = [1000 0 0 0 0];

    for zidx in 1:length(z)

        if length(vorts_3d[zidx]) != 0
            for v in 1:length(vorts_3d[zidx][:, 1])

                if vorts_3d[zidx][v, 5] == i
                    if current_vort[1] == 1000
                        current_vort = vorts_3d[zidx][v, :]';
                    else
                        current_vort = vcat(current_vort, vorts_3d[zidx][v, :]');
                    end
                end
            end 
        end

    end
    push!(vorts_label, current_vort)
end

NUM_VORTS = length(vorts_label)
vorts_label[1]

dx = x[2] - x[1];
dy = y[2] - y[1]; 


vort_linked = Array{Float64, 2}[];
while (length(vorts_label) > 0)
    
    v = vorts_label[1]
    if !vortex_boundary_bottom(v, X, dx, dy, dz)
        vortex_link_bottom(1, vorts_label, X)
    end
    if !vortex_boundary_top(v, X, dx, dy, dz)
        vortex_link_top(1, vorts_label, X)
    end


    vorts_label = sort_vorts_label(vorts_label);
    reverse!(vorts_label);
    label = vorts_label[end][1, 5];
    while ((length(vorts_label) > 0) && (vorts_label[end][1, 5] == label))
        pop = pop!(vorts_label);
        push!(vort_linked, pop);
    end
    reverse!(vorts_label)
    println(length(vorts_label))
    # Pop all v in vorts_label with same label as v
end

vort_linked[1]
vort_linked[2]
using VortexDistributions