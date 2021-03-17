module DifferentialGeometry
using ForwardDiff
using LinearAlgebra
export ∇, G, Jac, ∂, gauss, I, II, Curv, ν
""" 
Basic differential operators:
"""
∇(f, x) = ForwardDiff.gradient(f, x)
Jac(f, x) = ForwardDiff.jacobian(f, x)
∂(f, x, i) = ∇(f, x)[i]

"""
Let f: U -> C ⊂ S ⊂ IR^n be the parametrisation of an open hypersurface chunk C and x in U.
"""
G(f,x) = begin 
    A = Jac(f,x)
    return transpose(A)*A
end
H(f,x) = -transpose(Jac(gauss(f), x))*Jac(f, x)
L(f,x) = (G(f,x)^(-1)) * H(f,x)

"""
The Gauss map:
"""
ν(f, x) = begin
    A = Jac(f, x)
    B = [det([A _e(length(x) + 1, i)]) for i in 1:(length(x) + 1)]
    return B/norm(B)
end
gauss(f) = x -> ν(f,x)

"""
Fundamental forms I and II and Gaussian curvature:
"""
I(f,x,X,Y) = X⋅(G(f,x)*Y)
II(f,x,X,Y) = X⋅(H(f,x)*Y)
Curv(f,x) = det(L(f,x))

# Help functions
_e(n, i) = begin
    A = zeros(n)
    A[i] = 1
    return A
end

end # module
