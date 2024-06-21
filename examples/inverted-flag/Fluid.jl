using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
using Plots
using WriteVTK
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
wr = vtkWriter("WaterLily-Gismo"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop
    
    # Simulation parameters
    L,Re,U,Uref = 2^6,250,1,10.0
    center = SA[2L,2L]
    
    # coupling interface
    interface, body = initialize!(Uref,L,center;dir=[1,-1,-1,1])

    # this makes sur the spline extend to infinity
    ParametricBodies.notC¹(l::ParametricBodies.NurbsLocator,uv) = false

    # construct the simulation
    sim = Simulation((6L,4L),(U,0),L;U,ν=U*L/Re,body,T=Float64)
    store = Store(sim) # allows checkpointing

    # simulations time
    iter,every = 0,5 # for outputing VTK file

    # result storage
    results = []

    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(interface, sim, store)
        @show interface.deformation

        # measure the participant
        update!(interface, sim; center)

        # update the this participant
        step!(sim.flow, sim.pois, sim.body, interface)
        interface.forces .= 0.0 # scale
        # interface.forces .*= interface.U^2/interface.L # scale
        WaterLily.time(sim)<0.1sim.L && (interface.forces[2,Index(interface.quadPoint,3)] .= -4sim.L)
        
        # write data to the other participant
        writeData!(interface, sim, store)
        
        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            mod(iter,every)==0 && write!(wr, sim)
            iter += 1
            push!(results, sim.body.bodies[1].surf.pnts[:,end])
        end
    end
    close(wr)
    @show results
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")