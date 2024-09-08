@testset ExtendedTestSet "TimeSteppers" begin
  time = ConstantTimeStepper(0.0, 1.0, 0.1)
  @show time
  @test time.start_time ≈ 0.0
  @test time.end_time ≈ 1.0
  @test time.current_time ≈ 0.0
  @test time.current_time_step == 1
  @test time.Δt ≈ 0.1

  n = 1
  while time.current_time <= time.end_time
    Cthonios.step!(time)
    n = n + 1
    @test time.start_time ≈ 0.0
    @test time.end_time ≈ 1.0
    @test time.current_time ≈ (n - 1) * time.Δt
    @test time.current_time_step == n
    @test time.Δt ≈ 0.1
  end
end
