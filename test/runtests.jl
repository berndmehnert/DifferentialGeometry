using Test
using DifferentialGeometry
using LinearAlgebra

f(t)  = [cos(t[1]) sin(t[1])]
@test norm(ν(f,[0]) + [1; 0]) < 10^(-15)

k(a) = [cos(a) -sin(a); sin(a) cos(a)]
_g(a,t) = k(a)*[t 0;0 1/t]*k(-a)
g(x)= begin
    A = _g(x[1],x[2])
    return [A[1,1] A[1,2] A[2,2]]
end
h(x)=[1.0 0.0 2.0]⋅g(x)