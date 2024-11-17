@testset ExtendedTestSet "TimeSteppers" begin
  time = QuasiStatic(0.0, 1.0, 0.1)
  @show time
  @test time.start_time[1] ≈ 0.0
  @test time.end_time[1] ≈ 1.0
  @test time.current_time[1] ≈ 0.0
  @test time.current_time_step[1] == 1
  @test time.Δt[1] ≈ 0.1

  n = 1
  while time.current_time[1] <= time.end_time[1]
    Cthonios.step!(time)
    n = n + 1
    @test time.start_time[1] ≈ 0.0
    @test time.end_time[1] ≈ 1.0
    @test time.current_time[1] ≈ (n - 1) * time.Δt[1]
    @test time.current_time_step[1] == n
    @test time.Δt[1] ≈ 0.1
  end
end
