exr_GTEx[1,:]
exr_GTEx = DataFrame([[names(ex_GTEx)]; collect.(eachrow(ex_GTEx))], [:column; Symbol.(axes(ex_GTEx,2))])
exr_GTEx[1,2:length(exr_GTEx[1,:])]
col1 = Array(exr_GTEx[1,2:length(exr_GTEx[1,:])])
append(col1, "row_num")

cat("row_num",col1,dims=1)
cat("row_num",exr_GTEx[1,2:length(exr_GTEx[1,:])],dims=1)


names(exr_GTEx[2:length(exr_GTEx[1,:])])
names!(exr_GTEx[2:length(exr_GTEx[1,:])],exr_GTEx[1,2:length(exr_GTEx[1,:])])

# names(exr_GTEx[2:length(exr_GTEx[1,:])],exr_GTEx[1,2:length(exr_GTEx[1,:])])
# `getindex(df::DataFrame, col_inds::Union{AbstractVector, Regex, Not})` is deprecated, use `df[:, col_inds]` instead.
# │   caller = top-level scope at REPL[44]:1


###
names!(df, [Symbol("Col$i") for i in 1:size(df,2)])
