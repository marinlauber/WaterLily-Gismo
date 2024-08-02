using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
using WriteVTK
using FileIO,JLD2
include("../../src/Interface.jl")

# BC function for no slip on north and south face
function WaterLily.BC!(a,A,saveexit=false,perdir=())
    N,n = WaterLily.size_u(a)
    for j ∈ 1:n, i ∈ 1:n # face then component
        if j==2 && i==1 # Dirichlet cannot be imposed, must interpolate (u[i,j+1,1]+u[i,j,1])/2 = 0
            @WaterLily.loop a[I,i] = -a[I-δ(j,I),i] over I ∈ WaterLily.slice(N,N[j],j)
            @WaterLily.loop a[I,i] = -a[I+δ(j,I),i] over I ∈ WaterLily.slice(N,1,j)
            if all(A.≈0) # μ₀ case on the boundary
                @WaterLily.loop a[I,i] = A[i] over I ∈ WaterLily.slice(N,N[j],j)
                @WaterLily.loop a[I,i] = A[i] over I ∈ WaterLily.slice(N,1,j)
            end
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
wr = vtkWriter("Perp-Flap"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop
    
    # Simulation parameters
    L,Re,U,Lref,Uref = 2^5,100,1,1.0,10.0
    center = SA[2.95L,0]

    # coupling interface
    interface, body = initialize!(Uref,L,center;dir=[1,-1,-1])

    # slow ramp up of the velocity
    a = 2.0
    Ut(i,t::T) where T = i==1 ? convert(T,a*t/L+(1.0+tanh(31.4*(t/L-1.0/a)))/2.0*(1-a*t/L)) : zero(T) # velocity BC

    # construct the simulation
    sim = Simulation((8L,4L),Ut,L;U,ν=U*L/Re,body,T=Float64,exitBC=false)
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
        length(sim.flow.Δt)==2 && (interface.forces .= 0.0) # first time step
        interface.forces .*= interface.U^2 # scale

        # write data to the other participant
        writeData!(interface, sim, store)
        
        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            mod(iter,every)==0 && write!(wr, sim)
            iter += 1
            push!(results,[sum(@view(sim.flow.Δt[1:end-1]))*interface.U/sim.L,sim.body.bodies[2].curve.pnts[1,1]])
            println("tU/L=",round(sim_time(sim),digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
            println("time $(sum(@view(sim.flow.Δt[1:end-1]))*interface.U/sim.L)")
        end
    end
    close(wr)
    save("results.jld2","results",results)
    @show sim.flow.Δt
    @show sim.pois.n
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")

# using FileIO, JLD2, Plots
# a = load("/home/marin/Workspace/WaterLily-Gismo/examples/perpendicular-flap/results.jld2")["results"]
# pos = getindex.(a,2); time = getindex.(a,1)
# plot(time./2^5, pos./2^5.0.-3.0, xlabel="Time", ylabel="Y position", legend=:none)