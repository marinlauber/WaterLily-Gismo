using WaterLily
using ParametricBodies
using StaticArrays
using Plots
include("../../WaterLily/examples/TwoD_plots.jl")

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

Index(Qs,i) = sum(length.(Qs[1:i-1]))+1:sum(length.(Qs[1:i]))
function getInterfaceForces!(forces,flow::Flow{T},body::ParaBodies,quadPoints) where T
    for (i,b) in enumerate(body.bodies)
        I = Index(quadPoints,i)
        forces[:,I] .= reduce(hcat,[-1.0*ParametricBodies._pforce(b.surf,flow.p,s,zero(T),Val{false}()) for s ∈ quadPoints[i]])
    end
end

# parameters
L=2^5
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


# reset pressure
sim.flow.p.=0.0;

# intialize
t₀ = sim_time(sim)
duration = 100.
tstep = 0.2
vforces = []; pforces = []
# run
anim = @animate for tᵢ in range(t₀,t₀+duration;step=tstep)

    # update until time tᵢ in the background
    t = sum(sim.flow.Δt[1:end-1])
    while t < tᵢ*sim.L/sim.U
        mom_step!(sim.flow,sim.pois) # evolve Flow
        t += sim.flow.Δt[end]

        @inside sim.flow.p[I] = WaterLily.μ₀(sdf(sim.body,loc(0,I),0.0),1)*sim.flow.p[I] # fix pressure
        f = zeros(2); s = zeros(2)
        for b ∈ sim.body.bodies
            f .+= ParametricBodies.∮nds(sim.flow.p,b;N=16)
            s .+= ParametricBodies.∮τnds(sim.flow.u,b;N=16)
        end
        push!(pforces,f); push!(vforces,s)
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
vforces = vforces isa Matrix ? vforces : hcat(vforces...)
pforces = pforces isa Matrix ? pforces : hcat(pforces...)
# save gif
gif(anim, "test_forces.gif", fps=24)

# make fig
times = cumsum(sim.flow.Δt)[4:end-1]./sim.L
plt = plot(label="Pressure force",xlabel="Convective Time",ylabel="Force")
plot!(plt,times,-2pforces[1,4:end]./sim.L,label="Pressure")
plot!(plt,times,sim.flow.ν.*2vforces[1,4:end]/sim.L,label="Viscous")
xlims!(0,100); ylims!(0,2)