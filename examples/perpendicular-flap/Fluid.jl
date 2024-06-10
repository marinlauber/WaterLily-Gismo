using PreCICE
using WaterLily
using ParametricBodies
using StaticArrays
using Plots
using WriteVTK

# structure to store fluid state
struct Store
    uˢ::AbstractArray
    pˢ::AbstractArray
    xˢ::AbstractArray
    ẋˢ::AbstractArray
    function Store(sim::AbstractSimulation)
        xs = [copy(b.surf.pnts) for b in sim.body.bodies]
        vs = [copy(b.velocity.pnts) for b in sim.body.bodies]
        new(copy(sim.flow.u),copy(sim.flow.p),xs,vs)
    end
end
function store!(s::Store,sim::AbstractSimulation)
    s.uˢ .= sim.flow.u; s.pˢ .= sim.flow.p;
    for i ∈ 1:length(sim.body.bodies)
        s.xˢ[i] .= sim.body.bodies[i].surf.pnts
        s.ẋˢ[i] .= sim.body.bodies[i].velocity.pnts
    end
end
function revert!(s::Store,sim::AbstractSimulation)
    sim.flow.u .= s.uˢ; sim.flow.p .= s.pˢ; pop!(sim.flow.Δt)
    for i ∈ 1:length(sim.body.bodies)
        sim.body.bodies[i].surf.pnts     .= s.xˢ[i]
        sim.body.bodies[i].velocity.pnts .= s.ẋˢ[i]
    end
end

# unpack subarray of inceasing values of an hcat
function unpack(a)
    tmp=[a[1]]; ks=Vector{Number}[]
    for i ∈ 2:length(a)
        if a[i]>=a[i-1]
            push!(tmp,a[i])
        else
            push!(ks,tmp); tmp=[a[i]]
        end
    end
    push!(ks,tmp)
    return ks
end
function knotVectorUnpack(knots)
    knots = reshape(knots,reverse(size(knots)))[1,:]
    unpack(knots)
end
function getControlPoints(points, knots)
    points = reshape(points,reverse(size(points)))
    ncp = [length(knot)-count(i->i==0,knot) for knot in knots]
    @assert sum(ncp) == size(points,2) "Number of control points does not match the number of points"
    return [points[:,1+sum(ncp[1:i])-ncp[1]:sum(ncp[1:i])] for i in 1:length(ncp)]
end
function quadPointUnpack(quadPoints)
    quadPoints = reshape(quadPoints,reverse(size(quadPoints)))
    quadPoints = [filter(!isone,filter(!iszero,quadPoints[:,i]))[1] for i in 1:size(quadPoints,2)]
    unpack(quadPoints)
end
function getDeformation(points,knots)
    ncp = [length(knot)-count(i->i==0,knot) for knot in knots]
    return [points[:,1+sum(ncp[1:i])-ncp[1]:sum(ncp[1:i])] for i in 1:length(ncp)]
end
using ParametricBodies: _pforce

Index(Qs,i) = sum(length.(Qs[1:i-1]))+1:sum(length.(Qs[1:i]))
function getInterfaceForces!(forces,flow::Flow,body::ParaBodies,quadPoints)
    for (i,b) in enumerate(body.bodies)
        I = Index(quadPoints,i)
        forces[:,I] .= reduce(hcat,[-1.0*_pforce(b.surf,flow.p,s,WaterLily.time(flow),Val{false}()) for s ∈ quadPoints[i]])
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


# keyword aguments
if size(ARGS, 1) < 1
    @error(
        "ERROR: pass config path, example: julia solverdummy.jl ./precice-config.xml",
    )
    exit(1)
end

configFileName = ARGS[1]

# coupling
PreCICE.createParticipant("Fluid", configFileName, 0, 1)
println("I am participant WaterLily")


# writer for the sim
wr = vtkWriter("WaterLily-Gismo"; attrib=custom_attrib)

let # setting local scope for dt outside of the while loop
    
    PreCICE.initialize()

    # get the mesh verticies from the fluid solver
    (_, knots) = getMeshVertexIDsAndCoordinates("KnotMesh")
    knots = knotVectorUnpack(knots)
   
    # get the mesh verticies from the structural solver
    (ControlPointsID, ControlPoints) = getMeshVertexIDsAndCoordinates("ControlPointMesh")
    ControlPointsID = Array{Int32}(ControlPointsID)
    ControlPoints = getControlPoints(ControlPoints, knots)
    
    (quadPointID, quadPoint) = getMeshVertexIDsAndCoordinates("ForceMesh")
    forces = zeros(reverse(size(quadPoint))...)
    quadPointID = Array{Int32}(quadPointID)
    quadPoint = quadPointUnpack(quadPoint)
    
    # Simulation parameters
    L,Re,U,ϵ = 2^6,250,1,0.5

    # make the NURBS curve of the interface'
    center = SA[2L,0] # move to the correct position in the domain
    east = DynamicBody(NurbsCurve(MMatrix{size(ControlPoints[1])...}(ControlPoints[1])*L.+center,
                          knots[1],ones(size(ControlPoints[1],2))),(0,1))
    north= DynamicBody(NurbsCurve(MMatrix{size(ControlPoints[2])...}(reverse(ControlPoints[2],dims=2))*L.+center,
                          knots[2],ones(size(ControlPoints[2],2))),(0,1))
    west = DynamicBody(NurbsCurve(MMatrix{size(ControlPoints[3])...}(reverse(ControlPoints[3],dims=2))*L.+center,
                          knots[3],ones(size(ControlPoints[3],2))),(0,1))

    # make a combined body, carefull with winding direction of each curve
    body = ParametricBody([east,north,west])

    # this makes sur the spline extend to infinity
    ParametricBodies.notC¹(l::ParametricBodies.HashedLocator,uv) = false
    ParametricBodies.notC¹(l::ParametricBodies.NurbsLocator,uv) = false

    # construct the simulation
    sim = Simulation((6L,4L),(U,0),L;U,ν=U*L/Re,body,T=Float64)
    store = Store(sim) # allows checkpointing

    # simulations time
    t₀ = 0.0
    dt = dt_precice = PreCICE.getMaxTimeStepSize()
    iter,every = 0,20 # for outputing VTK file

    while PreCICE.isCouplingOngoing()
        
        # set time step
        dt_precice = PreCICE.getMaxTimeStepSize()
        dt = dt_precice # fix the time step to that of the precici-config file
        sim.flow.Δt[end] = dt*sim.L
           
        if PreCICE.requiresWritingCheckpoint()
            store!(store,sim)
        end

        # Read control point displacements
        readData = transpose(PreCICE.readData("ControlPointMesh", "ControlPointData", ControlPointsID, dt))
        deformation = getDeformation(readData,knots) # repack correctly
        for i in 1:length(sim.body.bodies)
            new = reverse(ControlPoints[i].+deformation[i],dims=2) # winding correctly
            i==1 && (new = ControlPoints[i].+deformation[i])       # this one is not reversed
            ParametricBodies.update!(sim.body.bodies[i],new.*sim.L.+center,dt*sim.L)
        end
 
        # solver update
        measure!(sim); mom_step!(sim.flow,sim.pois)
        @inside sim.flow.p[I] = WaterLily.μ₀(sdf(sim.body,loc(0,I),0.0),sim.ϵ)*sim.flow.p[I] # fix pressure
        getInterfaceForces!(forces,sim.flow,sim.body,quadPoint)
        forces = clamp(forces./sim.L,-sim.L,sim.L)
        
        # write the force at the integration points
        PreCICE.writeData("ForceMesh", "ForceData", quadPointID, permutedims(forces)) #reshape(forces,reverse(size(forces))))

        # do the coupling    
        dt_precice = PreCICE.advance(dt)
        
        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint()
            revert!(store,sim)
        else
            t₀ += dt*sim.L # this is not really usefull
            mod(iter,every)==0 && write!(wr, sim)
            iter += 1
        end
    end
    close(wr)
end
PreCICE.finalize()
println("WaterLily: Closing Julia solver...")