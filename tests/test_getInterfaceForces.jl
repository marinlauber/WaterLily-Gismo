using WaterLily
using ParametricBodies
using StaticArrays
using Plots
using ForwardDiff
include("../../WaterLily/examples/TwoD_plots.jl")
include("../src/Interface.jl")

function plot_integration_points!(body::ParametricBody;N=32)
    uv,q = ParametricBodies._gausslegendre(N,Float64)
    scale=0.5; uv=scale*(uv.+1); q=scale*q 
    for (s,w_) in zip(uv,q)
        x = sim.body.surf(s,0)
        nᵢ = ParametricBodies.norm_dir(sim.body.surf,s,0)
        nᵢ /= √(nᵢ'*nᵢ)
        scatter!([x[1]+nᵢ[1]],[x[2]+nᵢ[2]],color=:black,legend=:none)
    end
end

using LinearAlgebra: norm
function getInterfacedx(flow::Flow{T},body::ParaBodies,quadPoints,weights) where T
    ndx = zeros(length(body.bodies)*length(weights))
    for (i,b) in enumerate(body.bodies)
        I = Index(quadPoints,i)
        ndx[I] .= norm.([ForwardDiff.derivative(uv->b.surf(uv,0.0),uv)*w for (uv,w) in zip(quadPoints[i],weights)])
    end
    return ndx
end

function nds(curve::Function,t,lims::NTuple{2,T};N=64) where T
    # integrate NURBS curve to compute integral
    uv_, w_ = ParametricBodies._gausslegendre(N,T)
    # map onto the (uv) interval, need a weight scalling
    scale=(lims[2]-lims[1])/2; uv_=scale*(uv_.+1); w_=scale*w_ 
    [norm(ForwardDiff.derivative(uv->curve(uv,t),uv))*w for (uv,w) in zip(uv_,w_)]
end

function param_force(f::Function,curve::Function,t,lims::NTuple{2,T};N=64) where T
    # integrate NURBS curve to compute integral
    uv_, w_ = _gausslegendre(N,T)
    # map onto the (uv) interval, need a weight scalling
    scale=(lims[2]-lims[1])/2; uv_=scale*(uv_.+1); w_=scale*w_ 
    [f(uv) for uv in uv_]
end

# parameters
L=2^6
Re=80
U=1
center = SA[3L,3.01L]
radius = L/2

# NURBS points, weights and knot vector for a circle
cps1 = SA[1 1 0;0 1 1]
cps2 = SA[0 -1 -1; 1  1  0]
cps3 = SA[-1 -1  0; 0 -1 -1]
cps4 = SA[ 0  1 1;-1 -1 0]

weights = SA[1.,√2/2,1.]
knots =   SA[0,0,0,1/2,1,1,1]

ParametricBodies.notC¹(l::ParametricBodies.NurbsLocator,uv) = false

# make a nurbs curve
c1 = DynamicBody(NurbsCurve(MMatrix(cps1.*radius.+center),knots,weights),(0,1))
c2 = DynamicBody(NurbsCurve(MMatrix(cps2.*radius.+center),knots,weights),(0,1))
c3 = DynamicBody(NurbsCurve(MMatrix(cps3.*radius.+center),knots,weights),(0,1))
c4 = DynamicBody(NurbsCurve(MMatrix(cps4.*radius.+center),knots,weights),(0,1))

complete_circle = DynamicBody(NurbsCurve(MMatrix(SA[1 1 0 -1 -1 -1  0  1 1
                                                    0 1 1  1  0 -1 -1 -1 0]).*radius .+ center,
                                                    [0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1],
                                                    [1.,√2/2,1.,√2/2,1.,√2/2,1.,√2/2,1.]),(0,1))

# make a body and a simulation
Body = ParametricBody([c1,c2,c3,c4])
sim = Simulation((8L,6L),(U,0),L;U,ν=U*L/Re,body=Body,T=Float64,mem=Array)

measure_sdf!(sim.flow.σ,sim.body,0.0)
flood(sim.flow.σ)

# perimenter of the circle
# len = 0
# for b ∈ sim.body.bodies
    # len += ParametricBodies.integrate(b.surf,b.locate.lims)
# end
# println("Exact: $(L/2) , Numerical: ", round(len/2π,digits=4))

# uniform pressure field
WaterLily.apply!(x->x[2],sim.flow.p); h = zeros(2)
for b ∈ sim.body.bodies
    h .+= ParametricBodies.integrate(s->ParametricBodies._pforce(b.surf,sim.flow.p,
                                     s,0,Val{false}()),b.surf,0,b.locate.lims;N=16)/(π*(radius)^2)
end
println(h)
flood(sim.flow.p,shift=[-.5,-.5])
for b ∈ sim.body.bodies
    plot!(b.surf)
end
plot!()

using ParametricBodies: _gausslegendre,_pforce
# some bais checks
let
    quadPoints,weights = ParametricBodies._gausslegendre(32,Float64)
    quadPoints =  (quadPoints.+1)/2; weights ./= 2
    quadPoints = [quadPoints for _ in 1:length(sim.body.bodies)]
    ndx = getInterfacedx(sim.flow,sim.body,quadPoints,weights)
    forces = zeros(2,length(ndx))
    getInterfaceForces!(forces,sim.flow,sim.body,quadPoints,[1,1,1,1])
    ∑f = sum(forces.*ndx',dims=2)/(π*(radius)^2)
end
# @assert all(nds(sim.body.bodies[1].surf,0,sim.body.bodies[1].locate.lims;N=16) .≈ ndx[1:16])
# @assert all(nds(sim.body.bodies[1].surf,0,sim.body.bodies[1].locate.lims;N=16) .≈ ndx[17:32])
# @assert all(nds(sim.body.bodies[1].surf,0,sim.body.bodies[1].locate.lims;N=16) .≈ ndx[33:48])
# @assert all(nds(sim.body.bodies[1].surf,0,sim.body.bodies[1].locate.lims;N=16) .≈ ndx[49:64])

function pflow_circle(x,center,R)
    x = x .- center .+ 1.5 # whis is that here?
    r,θ = √sum(abs2,x),atan(x[2],x[1])
    r>R-10 ? 2*R^2/r^2*cos(2θ) - R^4/r^4 : 0.0
end
apply!(x->pflow_circle(x,center,radius),sim.flow.p)
flood(sim.flow.p,clims=(-3,1))

my_norm(a) = .√sum(abs2,a;dims=1)[:]

θ = 0:0.01:2π
plot(θ, ones(length(θ)), color=:black, proj=:polar, lw=2.0, label="Circle",dpi=1200)
plot!(θ,abs.(1.0.-4.0.*sin.(θ).^2).+0.5,color=:blue, lw=1.5,label="Analytical") # pressure force norm
for N ∈ [64,128,256]
    quadPoints,weights = ParametricBodies._gausslegendre(N,Float64)
    quadPoints =  (quadPoints.+1)/2; weights ./= 2
    quadPoints = [quadPoints for _ in 1:length(sim.body.bodies)]
    ndx = getInterfacedx(sim.flow,sim.body,quadPoints,weights)
    forces = zeros(2,length(ndx))
    getInterfaceForces!(forces,sim.flow,sim.body,quadPoints,[1,1,1,1])
    θ = reduce(vcat,[(π/2.0.*(q.+i.-1)) for (i,q) in enumerate(quadPoints)])
    plot!(θ, my_norm(forces).+0.5,label="getInterface N = $N")
    # uv = (_gausslegendre(4*N,Float64)[1].+1)./2.0
    # θ2 = uv*2π
    # get the uc from the previsou theta
    uv = θ./(2π)
    pf = reduce(hcat,[-1.0*_pforce(complete_circle.surf,sim.flow.p,s,0.0,Val{false}()) for s ∈ uv])
    plot!(θ, my_norm(pf).+0.5, label="Full circle N = $N")
end
plot!()
savefig("test_forces_256.png")

# intialize
t₀ = sim_time(sim)
duration = 10.
tstep = 0.2
vforces = []; pforces = []; preCICE = []

# reset pressure
sim.flow.p.=0.0;

# for thesting the integration
quadPoints,weights = ParametricBodies._gausslegendre(16,Float64)
quadPoints =  (quadPoints.+1)/2; weights ./= 2
quadPoints = [quadPoints for _ in 1:length(sim.body.bodies)]
ndx = getInterfacedx(sim.flow,sim.body,quadPoints,weights)
forces = zeros(2,length(ndx))
# run
anim = @animate for tᵢ in range(t₀,t₀+duration;step=tstep)

    # update until time tᵢ in the background
    t = sum(sim.flow.Δt[1:end-1])
    while t < tᵢ*sim.L/sim.U
        mom_step!(sim.flow,sim.pois) # evolve Flow
        t += sim.flow.Δt[end]

        @inside sim.flow.p[I] = WaterLily.μ₀(sdf(sim.body,loc(0,I),0.0),1)*sim.flow.p[I] # needed because ∮nds gets inside points
        f = zeros(2); s = zeros(2);
        for b ∈ sim.body.bodies
            f .+= ParametricBodies.∮nds(sim.flow.p,b;N=16)
            s .+= ParametricBodies.∮τnds(sim.flow.u,b;N=16)
        end
        getInterfaceForces!(forces,sim.flow,sim.body,quadPoints,[1,1,1,1])
        push!(pforces,f[1]); push!(vforces,s[1]); push!(preCICE,sum(forces.*ndx',dims=2)[1])
    end

    # flood plot
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u) * sim.L / sim.U
    contourf(clamp.(sim.flow.σ,-10,10)'|>Array,dpi=300,
             color=palette(:RdBu_11), clims=(-10,10), linewidth=0,
             aspect_ratio=:equal, legend=false, border=:none)
    for b ∈ sim.body.bodies
        plot!(b.surf; add_cp=true)
    end
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
# save gif
gif(anim, "test_forces.gif", fps=24)

# make fig
times = cumsum(sim.flow.Δt)[4:end-1]./sim.L
plt = plot(label="Pressure force",xlabel="Convective Time",ylabel="Force")
plot!(plt,times,-2pforces[4:end]./sim.L,label="Pressure")
plot!(plt,times, 2preCICE[4:end]./sim.L,label="Pressure PreCICE")
plot!(plt,times,sim.flow.ν.*2vforces[4:end]/sim.L,label="Viscous")
ylims!(0,2) # xlims!(0,100);
savefig(plt,"test_forces_circle.png")