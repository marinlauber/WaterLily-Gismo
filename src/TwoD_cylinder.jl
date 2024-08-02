using WaterLily
using Gismo
using BenchmarkTools
using StaticArrays
using Plots

# function non_alloc(object::Gismo.EigenMatrix)
#     return SVector{3}([1 for _ in 1:3])
# end

# function asSVector(object::Gismo.EigenMatrix,T)
#     return SVector{3,T}(unsafe_wrap(Array, Gismo.data(object), (Gismo.rows(object),Gismo.cols(object)); own = false))
# end

function asVector(object::Gismo.EigenMatrix)
    ptr = Gismo.data(object)
    SA{eltype(ptr)}[unsafe_load(ptr,1) unsafe_load(ptr,2) unsafe_load(ptr,3)]
end

# function asVector!(a,object::Gismo.EigenMatrix)
#     return unsafe_wrap(a, Gismo.data(object), (Gismo.rows(object),Gismo.cols(object)); own = false)
# end



struct GismoBody{T,G<:Gismo.Geometry,F<:Function} <: AbstractBody
    scale::T
    geom::G
    map::F
    function GismoBody(fname::String,map=(x,t)->x,scale=1,T=Float64)
        geom = Gismo.Geometry( fname )
        new{T,typeof(geom),typeof(map)}(scale,geom,map)
    end
end
function WaterLily.measure(body::GismoBody,x,t::T) where T
    # x = Float64[body.map(x,t)...,0]./body.scale # force 2D
    # get closest point on geometry
    d,uv = Gismo.closest(body.geom,x)
    # uv = Gismo.asMatrix(uv)
    # get normal at closest point
    n = asVector(Gismo.normal(body.geom,uv))
    # n = SA[n[1],n[2]] # force 2D
    # get velocity at closest point
    v = zero(n)
    return d,n/√(n'*n),v
end
function WaterLily.sdf(body::GismoBody,x,t)
    # x = Float64[body.map(x,t)...,0]./body.scale # force 2D
    # get closest point on geometry
    Gismo.closest(body.geom,x)[1]
end

# parameters
L = 2^4
U = 1
Re = 250.0

# generate a body
center = SA[2L,L]
body = GismoBody("surfaces/cylinder.xml", (x,t)->x.-center, 0.5)

# measure the distance, normal and velocity at a point
x = [0.5,0.5,0.0]
d,n,v = measure(body,x,0.0)
sdf(body,x,0.0)
@btime sdf($body,$x,$0)

d,uv = Gismo.closest(body.geom,[0.5,0.5,0.5])
# uv = Gismo.asMatrix(uv)
n = Gismo.normal(body.geom,uv)
# @btime asSVector($n)
# a = zeros(3)


@btime Gismo.closest($body.geom,$x)
@btime asVector($n)

# make a simulation
# sim = Simulation((4L,2L),(U,0),L;ν=U*L/Re,body,T=Float64)
# measure_sdf!(sim.flow.σ,body,0.0)
