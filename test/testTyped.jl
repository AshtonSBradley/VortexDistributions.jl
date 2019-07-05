# test typed version
include("typeversion.jl")
v1 = Vortex(.1,.2,1)
v2 = Vortex(randn(),randn(),rand([1,-1]))

function randvortex()
    return Vortex(randn(),randn(),rand([1,-1]))
end

randvortex()

varray = Vortex.(randn(2),randn(2),[1; -1])

v3 = [randn(2) randn(2) [1;-1]]
varray2 = Vortex.(v3[:,1],v3[:,2],v3[:,3])

varray3 = Vortex.(v3[1,:]...)

vort=Vortex(.1,.2,-1)
