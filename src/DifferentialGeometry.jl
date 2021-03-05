module DifferentialGeometry
using ForwardDiff
""" 
Basic differential operators:
"""
∇(f, v) = ForwardDiff.gradient(f, v)
Jac(f, v) = ForwardDiff.jacobian(f, v)

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

end # module