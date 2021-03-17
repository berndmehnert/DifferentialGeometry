using Test
using DifferentialGeometry
using LinearAlgebra

f(t)  = [cos(t[1]) sin(t[1])]
@test norm(ν(f,[0]) + [1; 0]) < 10^(-15)