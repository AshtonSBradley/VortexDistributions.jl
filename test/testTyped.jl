# test typed version
include("typeversion2.jl")
v1 = Vortex2(.1,.2,1)
v2 = Vortex2(randn(),randn(),rand([1,-1]))

function randvortex()
    return Vortex2(randn(),randn(),rand([1,-1]))
end

randvortex()

varray = Vortex2.(randn(2),randn(2),[1; -1])

v3 = [randn(2) randn(2) [1;-1]]
varray2 = Vortex2.(v3[:,1],v3[:,2],v3[:,3])

varray3 = Vortex2.(v3[1,:]...)
