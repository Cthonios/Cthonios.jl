# starting with 2D initially
# for now only works with tris most likely
function edges(domain::Domain)
  es = Matrix{Int}[]
  for section in domain.sections
    e_edges = section.fspace.ref_fe.edge_nodes
    for e in 1:num_elements(section.fspace)
      conn = connectivity(section.fspace, e)
      edges_temp = conn[e_edges]
      push!(es, edges_temp)
    end
  end
  return es
end

function edge_coordinates(coords::NodalField, edges)
  return map(x -> coords[:, x], edges)  
end

function edge_normals(edge_coords)
  # loop over elements
  normals = zeros(2, size(edge_coords[1], 3), length(edge_coords))
  for e in axes(edge_coords, 1)
    for edge in axes(edge_coords[e], 3)
      @views p1 = edge_coords[e][:, 1, edge]
      @views p2 = edge_coords[e][:, 2, edge]
      tangent = p2 .- p1
      normals[1, edge, e] = tangent[2]
      normals[2, edge, e] = -tangent[1]
      @views normals[:, edge, e] = normals[:, edge, e] / norm(normals[:, edge, e])
    end
  end
  normals
end

function edge_vectors(edge_coords)
  # normals = 
end