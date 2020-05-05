
function found_near(n)
    near = true
    for j in 1:n
        psi,vort = randVortexField(1)
        vortfound = findvortices(psi)
        vfdata = rawData(vortfound)
        vdata = rawData(vort)
        near *= isapprox(vdata,vfdata,rtol = 0.2)
    end
    return near
end

# unitary?
v0 = [.2 .4 1]
w0 = PointVortex(v0)
@test rawData(w0) == v0

v1 = [.2 .4 1;0.7 1.5 -1;-.3 1.2 1]
w1 = PointVortex(v1)
@test rawData(w1) == v1

w3 = randPointVortex(1000)
v3 = rawData(w3)
@test v3 == rawData(w3)

# single vortex creation and detection
@test found_near(1)
@test found_near(30)
