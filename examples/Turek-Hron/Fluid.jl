using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
using WriteVTK
using FileIO,JLD2
include("../../src/Interface.jl")

# velocity profile of Turek Hron
function uBC(i,xy,N,U1=1)
    x,y = @. xy .- 2
    i!=1 && return 0.0
    ((y < 0) && (y > N-1)) && return 0.0 # correct behaviour on ghost cells
    return 1.5*U1*y/(N-1)*(1.0-y/(N-1))/(0.5)^2
end

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
            if all(A.≈0)
                for s ∈ (1,2)
                    @WaterLily.loop a[I,i] = A[i] over I ∈ WaterLily.slice(N,s,j)
                end
                (!saveexit || i>1) && (@WaterLily.loop a[I,i] = A[i] over I ∈ WaterLily.slice(N,N[j],j)) # overwrite exit
            else
                for s ∈ (1,2)
                    @WaterLily.loop a[I,i] = uBC(i,loc(i,I),N[2],A[1]) over I ∈ WaterLily.slice(N,s,j)
                end
                (!saveexit || i>1) && (@WaterLily.loop a[I,i] = uBC(i,loc(i,I),N[2],A[1]) over I ∈ WaterLily.slice(N,N[j],j))
            end
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
wr = vtkWriter("Turek-Hron"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop
    
    # Simulation parameters
    D,Re,U,Uref = 2^4,250,1,2.0
    L = D*3.5 # diameter of the cylinder compared to the flap
    # center = SA[2.4D,1.85D]
    center = SA[3.4D,3.85D]
    
    # circle for Turek Hron
    # circle = AutoBody((x,t)->√sum(abs2,x.-SA[2D,1.95D])-D/2)
    circle = AutoBody((x,t)->√sum(abs2,x.-SA[3D,3.95D])-D/2)
    
    # coupling interface
    interface, body = initialize!(Uref,L,center;dir=[1,-1,-1,1]) #,curves=[circle])

    # construct the simulation
    # slow ramp up of the velocity
    a = 4.0
    Ut(i,t::T) where T = i==1 ? convert(T,a*t/L+(1.0+tanh(31.4*(t/L-1.0/a)))/2.0*(1-a*t/L)) : zero(T) # velocity BC
    # sim = Simulation((25D,4D),Ut,L;U,ν=U*D/Re,body,T=Float64,uλ=(i,x)->uBC(i,x,4D,1.0))
    sim = Simulation((16D,8D),Ut,L;U,ν=U*D/Re,body,T=Float64)
    store = Store(sim) # allows checkpointing
    write!(wr, sim) # write the initial state
    @show sim.body.bodies[1].curve.pnts
    @show sim.body.bodies[2].curve.pnts
    @show sim.body.bodies[3].curve.pnts
    @show sim.body.bodies[4].curve.pnts
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
        sim_time(sim)<2 && (interface.forces .= 0.0;
                            interface.forces[2,:] .-= 0.01sin(π*sim_time(sim))) # initial condition

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
let
    using FileIO,DataFrames,Plots
    data = DataFrame(load("Workspace/WaterLily-Gismo/examples/Turek-Hron/data/turek_hron_data.csv"),
                    ["time","disp"])
    plot(data.time,data.disp,xlims=(4,6.5),ylims=(-0.05,0.05),ls=:dash,
        ylabel="y-displacement [m]",xlabel="time [s]",lw=1,
        title="y-displacement of the flap tip",label="reference")
    a = load("Workspace/WaterLily-Gismo/examples/Turek-Hron/results.jld2")["results"]
    pos = getindex.(a,2); time = getindex.(a,1)
    plot!(time/16.0.+1.5, (pos.-sum(pos[end-20000:end])/20000)./16/3.5/2, label="WaterLily-Gismo")
    savefig("turek_hron.png")
end