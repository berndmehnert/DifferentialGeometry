using Test
using DifferentialGeometry
using LinearAlgebra

@testset "general tests" begin
    f(t)  = [cos(t[1]); sin(t[1])]  
    @test norm(ν(f,[0]) + [1; 0]) ≈ 0 atol = 1e-10
    @test det([Jac(f,[0]) ν(f,[0])]) > 0
    @test Curv(f,[0]) == 1.0
end

@testset "fancier example" begin
    f(x, γ) = [x[1], x[2], (γ+x[2]^2)/x[1]]
    N(x, γ) = [(γ+x[2]^2)/(x[1]^2), -2*x[2]/x[1], 1]
    H(x, γ) = [2*(γ+x[2]^2)/(x[1]^3) -2*x[2]/x[1]; -2*x[2]/x[1] 2/x[1]]/norm(N(x,γ))
    K(x, γ) = 4*γ/((x[1]^4)*norm(N(x,γ))^4)
    @test K([1.0, 2.0], 2.0) ≈ Curv(x -> f(x, 2.0), [1.0, 2.0]) atol = 1e-10
    @test K([1.1, 7.2], 1.2) ≈ Curv(x -> f(x, 1.2), [1.1, 7.2]) atol = 1e-10
    @test K([5.7, 2.2], 0.4) ≈ Curv(x -> f(x, 0.4), [5.7, 2.2]) atol = 1e-10
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
