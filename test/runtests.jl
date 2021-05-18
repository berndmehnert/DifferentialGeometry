using Test
using DifferentialGeometry
using LinearAlgebra

@testset "general tests" begin
    f(t)  = [cos(t[1]); sin(t[1])]  
    @test norm(ν(f,[0]) + [1; 0]) ≈ 0 atol = 1e-10
    @test det([Jac(f,[0]) ν(f,[0])]) > 0
    @test Curv(f,[0]) == 1.0
end

"""
k(a) = [cos(a) -sin(a); sin(a) cos(a)]
_g(a,t) = k(a)*[t 0;0 1/t]*k(-a)
g(x)= begin
    A = _g(x[1],x[2])
    return [A[1,1] A[1,2] A[2,2]]
end
h(x)=[1.0 0.0 2.0]⋅g(x)
p(x) = exp(1im * h(x))*<X,N>*ω(x) !!! TODO

for i in -π:0.1:π
    for j in 0.01:0.01:0.2
    push!(A,[i,j])
    end
end
B = map(x -> Curv(g,x), A)
for i in -π:0.1:π
    for j in 0.01:0.01:0.9
    push!(C,[i,j,Curv(g, [i,j])])
    end
end
s(x,y) = Curv(g, [x,y])
plot(-π:0.1:π, 0.01:0.01:0.2, s,st=:surface)
"""
