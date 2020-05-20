function findvortmask(ψ,x,y,R)
ψ = circmask(ψ,x,y,1.1*R)
nt,np,nn,vortices = find_vortices(x,y,ψ)
# remove vortices found outside mask boundary
for (i,xv) in enumerate(vortices[:,1]), yv in vortices[i,2]
    (norm([xv,yv]) > R) && (vortices[i,:] = [0. 0. 0.])
end

vortfound = Array{Float64, 2}[]
for (i,sv) in enumerate(vortices[:,3])
    (sv!=0) && push!(vortfound,vortices[i,:]')
end
    #currently will fail for more that one vortex
return vortfound[1]
end
