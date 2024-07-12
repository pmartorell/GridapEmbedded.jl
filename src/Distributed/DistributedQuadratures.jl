function ReferenceFEs.Quadrature(
  ::MomentFitted,
  cut::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  degree::Integer;
  kwargs...)

  cutfacets = cut_facets(cut,geo)
  Quadrature(momentfitted,cut,cutfacets,geo,degree;kwargs...)
end

function ReferenceFEs.Quadrature(
  ::MomentFitted,
  cut::DistributedEmbeddedDiscretization,
  degree::Integer;
  kwargs...)

  geo = get_geometry(cut)
  Quadrature(momentfitted,cut,geo,degree;kwargs...)
end

function CellData.Measure(
  trian::DistributedTriangulation{D},
  quad::Tuple{MomentFitted,Vararg};
  kwargs...) where D
  name, _args, _kwargs = quad
  cut,cutfacets,geo,_args... = _args
  model = get_background_model(cut)
  facet_gids = get_face_gids(model,D-1)
  measures = map(
    local_views(trian),
    local_views(cut),
    local_views(cutfacets),
    partition(facet_gids)) do trian,cut,cutfacets,facet_ids
      cutfacets = disable_ghost_facets(cutfacets,facet_ids)
      quad = name, (cut,cutfacets,geo,_args...), _kwargs
      Measure(trian,quad;kwargs...)
  end
  measures = move_ghost_contributions(measures,trian,cut,geo)
  DistributedMeasure(measures,trian)
end

function disable_facets!(cut::EmbeddedFacetDiscretization,mask)
  inactive = findall(lazy_map(!,mask))
  map(cut.ls_to_facet_to_inoutcut::Vector{Vector{Int8}}) do facet_to_inoutcut
    facet_to_inoutcut[inactive] .= typemax(Int8)
  end
end

function disable_facets(cut::EmbeddedFacetDiscretization,mask)
  cut = deepcopy(cut)
  disable_facets!(cut,mask)
  cut
end

function disable_ghost_facets(cut::EmbeddedFacetDiscretization,facet_ids)
  facet_mask = lazy_map(!=(0),local_to_own(facet_ids))
  disable_facets(cut,facet_mask)
end


function move_ghost_contributions(
  measures::AbstractArray{<:Measure},
  trian::DistributedTriangulation,
  cut::DistributedEmbeddedDiscretization)

  geo = get_geometry(cut)
  move_ghost_contributions(measures,trian,cut,geo)
end

function move_ghost_contributions(
  measures::AbstractArray{<:Measure},
  trian::DistributedTriangulation,
  cut::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry)

  bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cut,geo)
  move_ghost_contributions(measures,trian,bgcell_to_inoutcut)
end

function move_ghost_contributions(
  measures::AbstractArray{<:Measure},
  trian::DistributedTriangulation{D},
  bgcell_to_inoutcut::AbstractArray) where D

  # Setup graph
  model = get_background_model(trian)
  gids = get_cell_gids(model)
  graph = assembly_graph(partition(gids))

  # Setup snd ghost cells
  tcell_to_mcell = map(local_views(trian)) do trian
    get_glue(trian,Val{D}()).tface_to_mface
  end
  mcell_to_tcell = map(local_views(trian)) do trian
    get_glue(trian,Val{D}()).mface_to_tface
  end
  tcell_ghost = map(
    tcell_to_mcell,
    local_views(gids),
    graph.snd) do tcell_to_mcell,ids,snd

    tc_to_o = lazy_map(Reindex(local_to_owner(ids)),tcell_to_mcell)
    map(snd) do pid
      map((p)-> p == (pid),tc_to_o) |> findall
    end
  end
  tcells_snd = map(
    tcell_to_mcell,
    bgcell_to_inoutcut,
    local_views(gids),
    graph.snd) do tcell_to_mcell,bgcell_to_inoutcut,ids,snd

    tc_to_o = lazy_map(Reindex(local_to_owner(ids)),tcell_to_mcell)
    tc_to_ioc = lazy_map(Reindex(bgcell_to_inoutcut),tcell_to_mcell)
    map(snd) do pid
      map((p,ioc)-> p == (pid) && ioc == CUT,tc_to_o,tc_to_ioc) |> findall
    end
  end
  @show tcells_snd
  @show tcell_ghost
  gcells_snd = map(
    tcell_to_mcell,
    tcells_snd,
    local_views(gids)) do tcell_to_mcell,tcells,ids

    tc_to_g = lazy_map(Reindex(local_to_global(ids)),tcell_to_mcell)
    map(tcells) do tcells
      map(Reindex(tc_to_g),tcells)
    end
  end
  cell_weight_snd = map(tcells_snd,measures) do tcells,measure
    cell_weight = measure.quad.cell_weight
    map(Reindex(cell_weight),tcells)
  end

  # Exchange
  gcells_rcv = exchange(gcells_snd,graph) |> fetch
  cell_weight_rcv = exchange(cell_weight_snd,graph) |> fetch

  # Add contributions
  map(
    measures,
    partition(gids),
    mcell_to_tcell,
    gcells_rcv,
    cell_weight_rcv) do measure,ids,mcell_to_tcell,gcells,cell_weight

    @assert measure.quad.cell_weight.ptrs == 1:length(measure.quad.cell_weight.ptrs)
    map(gcells,cell_weight) do gcells,cell_weight
      lcells = map(Reindex(global_to_local(ids)),gcells)
      tcells = map(Reindex(mcell_to_tcell),lcells)
      measure.quad.cell_weight.values[tcells] .+= cell_weight
    end
  end

  # Remove ghost cells contributions
  @show tcell_ghost
  map(measures,tcell_ghost) do measure,tcells
    @assert measure.quad.cell_weight.ptrs == 1:length(measure.quad.cell_weight.ptrs)
    map(tcells) do tcells
      for tcell in tcells
        measure.quad.cell_weight.values[tcell] .= 0
      end
    end
  end

  measures
end

# TODO: Look how cut cells contribute
