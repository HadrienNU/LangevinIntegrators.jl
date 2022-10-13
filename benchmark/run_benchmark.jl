using PkgBenchmark
using LangevinIntegrators
using Dates

# res=benchmarkpkg(LangevinIntegrators; script="benchmark/benchmarks_kernel.jl")
res=benchmarkpkg(LangevinIntegrators)

export_markdown("results_$(Dates.format(now(),"ddmmyyyy")).md", res)


rm("colvar")
rm("plumed.log")
