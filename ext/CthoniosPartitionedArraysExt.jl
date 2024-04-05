module CthoniosPartitionedArraysExt

using Cthonios
using Exodus
using PartitionedArrays

struct QuasiStaticDomain
end

function exodus_pad(n_procs::Int32)

  if n_procs < 10
    pad_size = 1
  elseif n_procs < 100
    pad_size = 2
  elseif n_procs < 1000
    pad_size = 3
  elseif n_procs < 10000
    pad_size = 4
  else
    throw(ErrorException("Holy crap that's a big mesh. We need to check if we support that!"))
  end
  return pad_size
end

function QuasiStaticDomain(input_settings, rank)

  # first read nemesis file
  mesh_file = input_settings[:mesh][Symbol("file name")]
  # @info "Reading from $mesh_file"

  nemesis_file = mesh_file * ".nem"

  @info "Nemesis file = $nemesis_file"

  nem = ExodusDatabase(nemesis_file, "r")
  num_proc, num_proc_in_f, _ = Exodus.read_init_info(nem)


  init_global = Exodus.InitializationGlobal(nem)

  n_nodes_global = Exodus.num_nodes(init_global) |> Int64

  # now read the mesh from this proc
  mesh_file_rank = mesh_file * ".$num_proc" * ".$(lpad(rank - 1, exodus_pad(num_proc), '0'))"

  exo = ExodusDatabase(mesh_file_rank, "r")

  # get load balance stuff
  @show lb_params = Exodus.LoadBalanceParameters(exo, rank - 1)
  cmap_params = Exodus.CommunicationMapParameters(exo, lb_params, rank - 1)
  # @show cmap_params

  # @show Exodus.read_map(exo)
  @show Exodus.ProcessorNodeMaps(exo, rank)

  close(exo)
  close(nem)
end

struct ForwardProblem
end

function ForwardProblem(_, input_settings::Dict, common::Cthonios.CthoniosCommon, rank)
  domain_settings = input_settings[:domain]  
  domain = QuasiStaticDomain(domain_settings, rank)
end

function Cthonios.cthonios_main_mpi(
  input_file::String, n_procs::Int, distribute::F;
  verbose::Bool=true
) where F <: Function

  # set up ranks
  ranks = distribute(LinearIndices((n_procs,)))

  # log file set up
  # TODO more to do here. Need a parallel logger
  log_file_name = splitext(input_file)[1] * ".log"
  common = Cthonios.CthoniosCommon(log_file_name, Cthonios.CthoniosBackend(Cthonios.NoKABackend()))

  map(ranks) do rank
    println("I am proc $rank of $n_procs.")

    # read inputs file
    input_settings = Cthonios.parse_input_file(input_file)

    # TODO make this not hardcoded
    problem_temp = input_settings[:problems][1]
    
    prob = ForwardProblem(nothing, problem_temp, common, rank)
  end
end

end # module