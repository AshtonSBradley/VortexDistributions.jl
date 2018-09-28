function findnz(A)
I = findall(!iszero,A)
v = A[I]
ix = [I[i][1] for i in eachindex(I)]
iy = [I[i][2] for i in eachindex(I)]
return ix,iy,v
end
