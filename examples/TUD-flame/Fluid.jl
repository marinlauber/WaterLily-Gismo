using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
using WriteVTK
using FileIO,JLD2
include("../../src/Interface.jl")

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
wr = vtkWriter("TUD-flame"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop
    
    # Simulation parameters
    D,Re,U,Uref = 2^4,250,1,2.0
    L = D*3.5 # diameter of the cylinder compared to the flap
    # center = SA[2.4D,1.85D]
    center = SA[3.4D,3.85D]
  
    # coupling interface
    interface, body = initialize!(Uref,L,center;dir=[1,-1,-1,1])

    # construct the simulation
    # slow ramp up of the velocity
    a = 4.0
    Ut(i,t::T) where T = i==1 ? convert(T,a*t/L+(1.0+tanh(31.4*(t/L-1.0/a)))/2.0*(1-a*t/L)) : zero(T) # velocity BC
    # sim = Simulation((25D,4D),Ut,L;U,ν=U*D/Re,body,T=Float64,uλ=(i,x)->uBC(i,x,4D,1.0),exitBC=true)
    sim = Simulation((16D,8D),Ut,L;U,ν=U*D/Re,body,T=Float64)
    store = Store(sim) # allows checkpointing

    # simulations time
    iter,every = 0,20 # for outputing VTK file
    results = []

    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(interface, sim, store)

        # measure the participant
        update!(interface, sim; center)

        # update the this participant
        step!(sim.flow, sim.pois, sim.body, interface)
        interface.forces .*= interface.U^2 # scale
        sim_time(sim)<2 && (interface.forces .= 0.0) # initial condition

        # write data to the other participant
        writeData!(interface, sim, store)
        
        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            mod(iter,every)==0 && write!(wr, sim)
            push!(results,[sum(@view(sim.flow.Δt[1:end-1]))*interface.U/sim.L,sim.body.bodies[1].curve.pnts[2,end]])
            iter += 1
        end
    end
    close(wr)
    @show sim.pois.n
    @show sim.flow.Δt
    save("results.jld2","results",results)
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")

# using FileIO,DataFrames,Plots
# data = DataFrame(load("examples/Turek-Hron/data/turek_hron_data.csv"),
#                  ["time","disp"])
# plot(data.time,data.disp,xlims=(3,6.5),ylims=(-0.04,0.04),ls=:dash,
#      ylabel="y-displacement [m]",xlabel="time [s]",lw=1,
#      title="y-displaceent of the flap tip",label=:none)
# a = load("examples/Turek-Hron/results.jld2")["results"]
# pos = getindex.(a,2); time = getindex.(a,1)
# plot!(time, pos, label="WaterLily-Gismo")
# savefig("examples/Turek-Hron/data/turek_hron.png")