using CSV
using MetricForestCompletion

base_name = joinpath("~", "Projects", "clustering-data-v1", "mnist", "digits")
base_name = expanduser(base_name)
data = CSV.read(base_name * ".data.gz", CSV.Tables.matrix; header=false)
labels = CSV.read(base_name * ".labels0.gz", CSV.Tables.matrix; header=false)