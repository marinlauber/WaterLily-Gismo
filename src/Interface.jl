using ParametricBodies
using WaterLily
using StaticArrays

# structure to store fluid state
struct Store
    uˢ::AbstractArray
    pˢ::AbstractArray
    b::Vector{Bodies}
    function Store(sim::Simulation)
        new(copy(sim.flow.u),copy(sim.flow.p),[copy(sim.body)])
    end
end
function store!(s::Store,sim::Simulation)
    s.uˢ .= sim.flow.u; s.pˢ .= sim.flow.p
    s.b[1] = copy(sim.body)
end
function revert!(s::Store,sim::Simulation)
    sim.flow.u .= s.uˢ; sim.flow.p .= s.pˢ; pop!(sim.flow.Δt)
    pop!(sim.pois.n); pop!(sim.pois.n) # pop predictor and corrector
    sim.body = s.b[1] # nice and simple
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
using ParametricBodies: _pforce, _vforce
"""
    getInterfaceForces!(forces,flow::Flow,body::Bodies,quadPoints)

Compute the interface forces at the quadrature points `quadPoints` for each body in `body.bodies`.
"""
Index(Qs,i) = sum(length.(Qs[1:i-1]))+1:sum(length.(Qs[1:i]))
function getInterfaceForces!(forces,flow::Flow{T},body::Bodies,quadPoints,dir) where T
    for (i,b) in enumerate(body.bodies[1:length(quadPoints)]) # only select the active curves
        I = Index(quadPoints,i)
        fi = reduce(hcat,[-1.0*_pforce(b.curve,flow.p,s,zero(T),Val{false}()) for s ∈ quadPoints[i]])
        dir[i] != 1 && (fi = reverse(fi;dims=2))
        forces[:,I] .= fi
    end
end

struct CouplingInterface
    U::Float64
    # L::Float64
    ControlPointsID::AbstractArray
    ControlPoints::AbstractArray
    quadPointID::AbstractArray
    quadPoint::AbstractArray
    forces::AbstractArray
    deformation::AbstractArray
    knots::AbstractArray
    dt::Vector{Float64}
    dir::Vector{Int16}
    N::Int16
end

using PreCICE
function initialize!(U,L,center;KnotMesh="KnotMesh",ControlPointMesh="ControlPointMesh",
                     ForceMesh="ForceMesh",dir=nothing,curves=nothing)
    
    # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # coupling
    PreCICE.createParticipant("Fluid", configFileName, 0, 1)
    println("I am participant WaterLily")

    # initilise PreCICE
    PreCICE.initialize()

    # get the mesh verticies from the fluid solver
    (_, knots) = getMeshVertexIDsAndCoordinates(KnotMesh)
    knots = knotVectorUnpack(knots)
   
    # get the mesh verticies from the structural solver
    (ControlPointsID, ControlPoints) = getMeshVertexIDsAndCoordinates(ControlPointMesh)
    ControlPointsID = Array{Int32}(ControlPointsID)
    ControlPoints = getControlPoints(ControlPoints, knots)
    deformation = copy(ControlPoints)
    isnothing(dir) ? (direction = ones(length(ControlPoints))) : direction = dir
    
    # get the quad points in parameter space
    (quadPointID, quadPoint) = getMeshVertexIDsAndCoordinates(ForceMesh)
    forces = zeros(reverse(size(quadPoint))...)
    quadPointID = Array{Int32}(quadPointID)
    quadPoint = quadPointUnpack(quadPoint)

    dt = PreCICE.getMaxTimeStepSize()
    
    # construct the interface curves
    bodies = AbstractBody[]; ops = Function[]
    for (i,(cps,knot)) in enumerate(zip(ControlPoints,knots))
        direction[i] != 1 && (cps = reverse(cps;dims=2))
        cps = SMatrix{2,size(cps,2)}(cps)
        knot = SVector{length(knot)}(knot)
        weights = SA[ones(size(cps,2))...]
        push!(bodies,DynamicNurbsBody(NurbsCurve(cps*L.+center,knot,weights)))
        push!(ops, ∩) # always interset with the next curve
    end

    # return coupling interface
    interface = CouplingInterface(U, ControlPointsID, ControlPoints, quadPointID, quadPoint, forces,
                                  deformation, knots, [dt], direction, length(bodies))
    
    # add some passive curves if we want
    !isnothing(curves) && for crv in curves
        push!(bodies,crv); push!(ops, +) # always union with the next curve
        println("Adding a curve to the stack...")
    end

    # return the interface and the body
    return interface, Bodies(bodies, ops)
end

function readData!(interface::CouplingInterface,sim::Simulation,store::Store)

    # set time step
    dt_precice = PreCICE.getMaxTimeStepSize()
    interface.dt[end] = dt_precice# min(sim.flow.Δt[end]*sim.U/sim.L, dt_precice) # min physical time step
    sim.flow.Δt[end] = interface.dt[end]*sim.L/sim.U # numerical time step

    if PreCICE.requiresWritingCheckpoint()
        store!(store,sim)
    end

    # Read control point displacements
    readData = transpose(PreCICE.readData("ControlPointMesh", "ControlPointData", 
                                          interface.ControlPointsID, interface.dt[end]))
    interface.deformation .= getDeformation(readData,interface.knots) # repack correctly
end

function update!(interface::CouplingInterface,sim::Simulation;center=0)
    # update the geom as this has not been done yet
    for i in 1:interface.N
        new = interface.ControlPoints[i].+interface.deformation[i]
        interface.dir[i] != 1 && (new = reverse(new;dims=2))
        new = SMatrix{size(new)...}(new.*sim.L.+center)
        # time step is the (numerical) time between data exchange
        sim.body.bodies[i] = ParametricBodies.update!(sim.body.bodies[i],new,sim.flow.Δt[end])
    end
    # solver update
    WaterLily.measure!(sim)
end

using ParametricBodies: _pforce, _vforce
function step!(flow::Flow,pois::AbstractPoisson,body,interface::CouplingInterface)
    mom_step!(flow,pois)
    @inside flow.p[I] = WaterLily.μ₀(sdf(body,loc(0,I),0.0),1)*flow.p[I] # fix pressure
    getInterfaceForces!(interface.forces,flow,body,interface.quadPoint,interface.dir)
end

function writeData!(interface::CouplingInterface,sim::Simulation,store::Store)
    # write the force at the integration points
    PreCICE.writeData("ForceMesh", "ForceData", interface.quadPointID, permutedims(interface.forces))

    # do the coupling    
    interface.dt[end] = PreCICE.advance(interface.dt[end])
    
    # read checkpoint if required or move on
    if PreCICE.requiresReadingCheckpoint()
        revert!(store,sim)
    end
end