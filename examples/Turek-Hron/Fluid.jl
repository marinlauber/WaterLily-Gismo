using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
using Plots
using WriteVTK
include("../../src/Interface.jl")


# velocity profile of Turek Hron
# function uλ(i,xy)
#     x,y = @. xy .- 2
#     i!=1 && return 0.0
#     ((y < 0) && (y > 4D-1)) && return 0.0 # correct behaviour on ghost cells
#     return 1.5*U*y/(4D-1)*(1.0-y/(4D-1))/(0.5)^2
# end

# # overwrite the momentum function so that we get the correct BC
# @fastmath function WaterLily.mom_step!(a::Flow,b::AbstractPoisson)
#     a.u⁰ .= a.u; WaterLily.scale(a,0)
#     # predictor u → u'
#     WaterLily.conv_diff!(a.f,a.u⁰,a.σ,ν=a.ν)
#     WaterLily.BDIM!(a); BC_TH!(a.u,a.U)
#     WaterLily.project!(a,b); BC_TH!(a.u,a.U)
#     # corrector u → u¹
#     WaterLily.conv_diff!(a.f,a.u,a.σ,ν=a.ν)
#     WaterLily.BDIM!(a); BC_TH!(a.u,a.U)
#     WaterLily.project!(a,b); WaterLily.scale(a,0.5); BC_TH!(a.u,a.U)
#     push!(a.Δt,WaterLily.CFL(a))
# end

# # BC function using the profile
# function BC_TH!(a,A,f=1)
#     N,n = WaterLily.size_u(a)
#     for j ∈ 1:n, i ∈ 1:n
#         if i==j # Normal direction, impose profile on inlet and outlet
#             for s ∈ (1,2,N[j])
#                 @WaterLily.loop a[I,i] = f*uλ(i,loc(i,I)) over I ∈ WaterLily.slice(N,s,j)
#             end
#         else  # Tangential directions, interpolate ghost cell to homogeneous Dirichlet
#             @WaterLily.loop a[I,i] = -a[I+δ(j,I),i] over I ∈ WaterLily.slice(N,1,j)
#             @WaterLily.loop a[I,i] = -a[I-δ(j,I),i] over I ∈ WaterLily.slice(N,N[j],j)
#         end
#     end
# end

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
    D = L÷2
    center = SA[4D,2D]
    
    # circle for Turek Hron
    cps = SA[1 1 0 -1 -1 -1  0  1 1.
             0 1 1  1  0 -1 -1 -1 0]*D/2 .+ center .-SA[D/2,0]
    weights = [1.,√2/2,1.,√2/2,1.,√2/2,1.,√2/2,1.]
    knots =   [0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1]

    # make a nurbs circle
    circle = DynamicBody(NurbsCurve(MMatrix(cps),knots,weights),(0,1))

    # coupling interface
    interface, body = initialize!(Uref,L,center;dir=[1,-1,-1,1],curves=[circle])

    # this makes sur the spline extend to infinity
    ParametricBodies.notC¹(l::ParametricBodies.NurbsLocator,uv) = false

    # construct the simulation
    sim = Simulation((11D,4D),(U,0),L;U,ν=U*L/Re,body,T=Float64)
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