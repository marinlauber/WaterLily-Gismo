using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
using Plots
using WriteVTK
using FileIO,JLD2
include("../../src/Interface.jl")

# overwrite the momentum function so that we get the correct BC
@fastmath function WaterLily.mom_step!(a::Flow,b::AbstractPoisson)
    a.u⁰ .= a.u; WaterLily.scale_u!(a,0)
    # predictor u → u'
    WaterLily.conv_diff!(a.f,a.u⁰,a.σ,ν=a.ν,perdir=a.perdir)
    WaterLily.BDIM!(a); BC_!(a.u,a.U)
    a.exitBC && WaterLily.exitBC!(a.u,a.u⁰,a.U,a.Δt[end]) # convective exit
    WaterLily.project!(a,b); BC_!(a.u,a.U)
    # corrector u → u¹
    WaterLily.conv_diff!(a.f,a.u,a.σ,ν=a.ν,perdir=a.perdir)
    WaterLily.BDIM!(a); WaterLily.scale_u!(a,0.5); BC_!(a.u,a.U)
    WaterLily.project!(a,b,0.5); BC_!(a.u,a.U)
    push!(a.Δt,WaterLily.CFL(a))
end

# BC function for no slip on noth and south face
function BC_!(a,A;saveexit=true)
    N,n = WaterLily.size_u(a)
    for j ∈ 1:n, i ∈ 1:n # face then component
        if j==2 && i==1 # Dirichlet cannot be imposed, must interpolate (u[i,j+1,1]+u[i,j,1])/2 = 0
            @WaterLily.loop a[I,i] = -a[I-δ(j,I),i] over I ∈ WaterLily.slice(N,N[j],j)
            @WaterLily.loop a[I,i] = -a[I+δ(j,I),i] over I ∈ WaterLily.slice(N,1,j)
        elseif i==j # Normal direction, homoheneous Dirichlet
            @WaterLily.loop a[I,i] = A[i] over I ∈ WaterLily.slice(N,1,j)
            @WaterLily.loop a[I,i] = A[i] over I ∈ WaterLily.slice(N,2,j)
            (!saveexit || i>1) && (@WaterLily.loop a[I,i] = A[i] over I ∈ WaterLily.slice(N,N[j],j))
        else  # Tangential directions, interpolate ghost cell to homogeneous Dirichlet
            @WaterLily.loop a[I,i] = a[I+δ(j,I),i] over I ∈ WaterLily.slice(N,1,j)
            @WaterLily.loop a[I,i] = a[I-δ(j,I),i] over I ∈ WaterLily.slice(N,N[j],j)
        end
    end
end

# make a writer with some attributes
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body); 
                        a.flow.σ |> Array;)
vorticity(a::Simulation) = (@inside a.flow.σ[I] = 
                            WaterLily.curl(3,I,a.flow.u)*a.L/a.U;
                            a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;
mu0(a::Simulation) = a.flow.μ₀ |> Array;

custom_attrib = Dict(
    "u" => velocity,
    "p" => pressure,
    "d" => _body,
    "ω" => vorticity,
    "v" => _vbody,
    "μ₀" => mu0
)# this maps what to write to the name in the file

# writer for the sim
wr = vtkWriter("WaterLily-Gismo"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop
    
    # Simulation parameters
    L,Re,U,Lref,Uref = 2^6,250,1,1.0,10.0
    center = SA[2.95L,0]
    
    # coupling interface
    interface, body = initialize!(Uref,L,center;dir=[1,-1,-1])

    # this makes sur the spline extend to infinity
    ParametricBodies.notC¹(l::ParametricBodies.NurbsLocator,uv) = false

    # construct the simulation
    sim = Simulation((6L,4L),(U,0),L;U,ν=U*L/Re,body,T=Float64)
    sim.flow.Δt[end] = 0.4
    store = Store(sim) # allows checkpointing

    # simulations time
    iter,every = 0,5 # for outputing VTK file

    # result storage
    results = []

    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(interface, sim, store)

        # measure the participant
        update!(interface, sim; center)

        # update the this participant
        step!(sim.flow, sim.pois, sim.body, interface)
        interface.forces .*= interface.U^2 # scale

        # interface.forces .= 0.0 # scale
        # interface.forces[1,:] .= interface.U^2 

        # write data to the other participant
        writeData!(interface, sim, store)
        
        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            mod(iter,every)==0 && write!(wr, sim)
            iter += 1
            push!(results,[interface.dt[end]*sim.L/sim.U,sim.body.bodies[2].surf.pnts[1,1]])
            println("tU/L=",round(sim_time(sim),digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
        end
    end
    close(wr)
    save("results.jld2","results",results)
    @show sim.flow.Δt
    @show sim.pois.n
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")

# using FileIO, JLD2,Plots
a = load("/home/marin/Workspace/WaterLily-Gismo/examples/perpendicular-flap/results.jld2")["results"]
pos = getindex.(a,2); time = cumsum(getindex.(a,1))
plot(time./2^6,pos./2^6.0.-3.0, xlabel="Time", ylabel="Y position", legend=:none)