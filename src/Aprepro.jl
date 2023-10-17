function aprepro(file_name::String)
  file_out = file_name * ".aprepro"

  aprepro_out = @capture_out @capture_err aprepro_exe() do exe
    run(`$exe --comment=\# $file_name $file_out`, wait=true)
  end

  @info aprepro_out

  return file_out
end 