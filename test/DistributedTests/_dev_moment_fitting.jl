using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test
using Gridap.ReferenceFEs

using GridapEmbedded.CSG

# using GridapEmbedded: AgFEMSpace

function circle_geometry(ranks,parts,cells)
  L = 1
  p0 = Point(0.0,0.0)
  pmin = p0-L/2
  pmax = p0+L/2
  R = 0.35
  geo = disk(R,x0=p0)
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,cells)
  bgmodel,geo
end

function remotes_geometry(ranks,parts,cells)
  x0 = Point(0.05,0.05)
  d1 = VectorValue(0.9,0.0)
  d2 = VectorValue(0.0,0.1)
  geo1 = quadrilateral(;x0=x0,d1=d1,d2=d2)

  x0 = Point(0.15,0.1)
  d1 = VectorValue(0.25,0.0)
  d2 = VectorValue(0.0,0.6)
  geo2 = quadrilateral(;x0=x0,d1=d1,d2=d2)
  geo = union(geo1,geo2)

  domain = (0, 1, 0, 1)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells)
  bgmodel,geo
end



  distribute = Vector
  threshold=1
  parts = (2,2)
  cells=(6,6)
  geometry=:circle

  ranks = distribute(LinearIndices((prod(parts),)))

  u(x) = x[1] - x[2]
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  geometries = Dict(
    :circle => circle_geometry,
    :remotes => remotes_geometry,
  )

  bgmodel,geo = geometries[geometry](ranks,parts,cells)

  writevtk(bgmodel.models[4],"model0")

  D = 2
  cell_meas = map(get_cell_measure∘Triangulation,local_views(bgmodel))
  meas = map(first,cell_meas) |> PartitionedArrays.getany
  h = meas^(1/D)

  cutgeo = cut(bgmodel,geo)
  cutgeo_facets = cut_facets(bgmodel,geo)

  strategy = AggregateCutCellsByThreshold(threshold)
  bgmodel,cutgeo,aggregates = aggregate(strategy,cutgeo)

  Ω_bg = Triangulation(bgmodel)
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  order = 1
  degree = 2*order
  dΩ = Measure(Ω,degree)
  quad = Quadrature(momentfitted,cutgeo,degree;in_or_out=IN)
  dΩ_act = Measure(Ω_act,quad)

  @show sum(∫(1)dΩ)
  @show sum(∫(1)dΩ_act)

  dΩ_act.measures[1].quad.cell_weight[1]
  dΩ_act.measures[1].quad.cell_point[1]
  sum(∫(1)dΩ_act.measures[1])
  sum(∫(1)dΩ.measures[1])

writevtk(Ω_bg,"trian")


  writevtk(Ω_act.trians[1],"act_1")
  writevtk(Γ,"bnd")

  p = 1
  fids = get_face_gids(bgmodel,1).partition[p]

  cf = cutgeo_facets.discretizations[p]
  cf = GridapEmbedded.Distributed.disable_ghost_facets(cf,fids)
  t = SkeletonTriangulation(cf,IN)
  writevtk(t,"sk_1")
  t = BoundaryTriangulation(cf,PHYSICAL_IN)
  writevtk(t,"bt_1")


  #TODO:
  #   - debug moment fitting integration when not all cells are cut
  #   - revise moment fitting for non-cut cells
  #   - create moment fitting tests
  #   - remove this file
