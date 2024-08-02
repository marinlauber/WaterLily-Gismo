using WaterLily
using StaticArrays
using Plots
include("../../../WaterLily/examples/TwoD_plots.jl")

WaterLily.L₂(ml::MultiLevelPoisson) = WaterLily.L₂(ml.levels[1])
WaterLily.L∞(ml::MultiLevelPoisson) = WaterLily.L∞(ml.levels[1])
using WaterLily: L∞

function WaterLily.solver!(ml::MultiLevelPoisson;tol=2e-4,itmx=32,L2=true)
    p = ml.levels[1]
    WaterLily.residual!(p); r₂ = L2 ? L₂(p) : L∞(p)
    nᵖ=0
    while nᵖ<itmx
        WaterLily.Vcycle!(ml)
        WaterLily.smooth!(p); r₂ = L2 ? L₂(p) : L∞(p)
        nᵖ+=1; r₂<tol && break
    end
    WaterLily.perBC!(p.x,p.perdir)
    push!(ml.n,nᵖ);
end

function test_1(pois::AbstractPoisson;tol=10eps(T),itmx=32)
    # source term for the solution
    # u := cos(2πx/L) cos(2πy/L)
    uₑ = copy(pois.x); apply!(x->cos(2π*x[1]/L)*cos(2π*x[2]/L),uₑ)
    # ∂u²/∂x² + ∂u²/∂y² = f = -8π²/L² cos(2πx/L) cos(2πy/L)
    apply!(x->-8π^2/L^2*cos(2π*x[1]/L)*cos(2π*x[2]/L),pois.z)

    # solver
    solver!(pois;tol=tol,itmx=itmx)
    
    # L₂ and L∞ nor of the error
    L₂(uₑ[R]-pois.x[R])/length(uₑ), maximum(abs.(uₑ[R]-pois.x[R]))
end

function test_2(pois::AbstractPoisson;tol=10eps(T),itmx=32)
    # source term for the solution
    # u := 
    uₑ = copy(pois.x); apply!(x->,uₑ)
    # ∂u²/∂x² + ∂u²/∂y² = f = 
    apply!(x->,pois.z)

    # solver
    solver!(pois;tol=tol,itmx=itmx)
    
    # L₂ and L∞ nor of the error
    L₂(uₑ[R]-pois.x[R])/length(uₑ), maximum(abs.(uₑ[R]-pois.x[R]))
end

# domain and fields
L = 32
N,D,T = (L,L),2,Float64
Ng = N .+ 2
x,z = zeros(T, Ng) |> Array, zeros(T, Ng) |> Array
μ⁰ = ones(T, (Ng..., D)) |> Array

# apply zero Neumann BC
BC!(μ⁰,zeros(2))

MG = true
L2 = false
@show MG,L2

# construct Poisson problem
pois = MG ? MultiLevelPoisson(x,μ⁰,z) : Poisson(x,μ⁰,z)
R = inside(pois.x)

test_1(pois)
