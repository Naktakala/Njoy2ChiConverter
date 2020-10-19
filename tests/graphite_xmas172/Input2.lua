--########################################## Mesh
chiMeshHandlerCreate()
N = 100
L = 10.0
ds = L/N
nodes={}
for i=0,N do 
  table.insert(nodes,i*ds)
end 
chiMeshCreateUnpartitioned1DOrthoMesh(nodes)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();

vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--########################################## Cross-section
xs_grph = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_grph,CHI_XSFILE,"../../output/ENDF-B-VII-1/xmas172/Cnat_graphite.csx")
xs = chiPhysicsTransportXSGet(xs_grph)
num_groups = xs["G"]
for k,v in pairs(xs) do 
  print(k,v)
end
print("number of groups: ",num_groups)

--########################################## Materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Graphite");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              EXISTING,xs_grph)
src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[21] = 1.0
chiPhysicsMaterialSetProperty(materials[1],
                              ISOTROPIC_MG_SOURCE,
                              FROM_ARRAY,src)

--########################################## Solver
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)
--
--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end
--========== ProdQuad
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE,24)
pquad1 = chiCreateProductQuadrature(GAUSS_LEGENDRE,48)
--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,49)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad1)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-8)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
-- chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
-- chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--========== Groupset def
gs1 = chiLBSCreateGroupset(phys1)
cur_gs = gs1
chiLBSGroupsetAddGroups(phys1,cur_gs,50,80)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad1)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
-- chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
-- chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

gs1 = chiLBSCreateGroupset(phys1)
cur_gs = gs1
chiLBSGroupsetAddGroups(phys1,cur_gs,81,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,7)
--
chiLBSInitialize(phys1)
chiLBSExecute(phys1)
--
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

for g=1,num_groups do 
  ffi = chiFFInterpolationCreate(VOLUME)
  curffi = ffi

  function Adder(ff_value, mat_id)
    ret_val = ff_value;   --Or some computation
    return ret_val
  end

  chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM_LUA,"Adder")
  chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[g])
  chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
  chiFFInterpolationInitialize(curffi)
  chiFFInterpolationExecute(curffi)
  print(g-1,chiFFInterpolationGetValue(curffi))
end

