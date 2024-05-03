
function CellData.Measure(
  t::DistributedTriangulation,
  quad::Tuple{MomentFitted,Vararg};
  kwargs...)
  @notimplemented
  name, _args, _kwargs = quad
  cut,cutfacets,_args... = _args
  measures = map(
    local_views(t),
    local_views(cut),
    local_views(cutfacets)) do trian,cut,cutfacets
      # TODO: deactivate ghost facets, e.g., set to CUT or typemax(Int8)
      cutfacets = disable_ghost_facets(cutfacets)
      quad = name, (cut,cutfacets,_args...), _kwargs
      Measure(trian,quad;kwargs...)
  end
  # add contributions from ghost cells ( measures )
  # i.e., sum the weights to the own cell (need a reverse map)
  # remove ghost cells ( measures )
  DistributedMeasure(measures)
end

function disable_facets(cut::EmbeddedFacetDiscretization,active_cells)
  @notimplemented
  map(cut.ls_to_facet_to_inoutcut::Vector{Vector{Int8}}) do facet_to_inoutcut
  end
end

function disable_ghost_facets(cut::DistributedEmbeddedDiscretization)
  @notimplemented
end

function add_ghost_contributions(
  measures::AbstractArray{<:Measure},
  gids::PRange)
  @notimplemented
end
