using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
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
    L,Re,U,Uref = 2^5,250,1,1.0
    center = SA[4L,4L]
    
    # coupling interface
    interface, body = initialize!(Uref,L,center;dir=[1,-1,-1,1])

    # slow ramp up of the velocity
    a = 1.0
    Ut(i,t::T) where T = i==1 ? convert(T,a*t/L+(1.0+tanh(31.4*(t/L-1.0/a)))/2.0*(1-a*t/L)) : zero(T) # velocity BC

    # construct the simulation
    sim = Simulation((12L,8L),Ut,L;U,ν=U*L/Re,body,T=Float64)
    store = Store(sim) # allows checkpointing

    # simulations time
    iter,every = 0,25 # for outputing VTK file

    while PreCICE.isCouplingOngoing()

        # read the data from the other participant
        readData!(interface, sim, store)

        # update the participant
        update!(interface, sim; center)

        # step this participant
        step!(sim.flow, sim.pois, sim.body, interface)
        sim_time(sim)<1.5 && (interface.forces[2,Index(interface.quadPoint,3)] .-= 1)
        
        # write data to the other participant
        writeData!(interface, sim, store)
        
        # if we have finshed this time step, save stiff
        if PreCICE.isTimeWindowComplete()
            mod(iter,every)==0 && write!(wr, sim)
            iter += 1
        end
    end
    close(wr)
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")