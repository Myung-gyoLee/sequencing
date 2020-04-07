exr_GTEx[1,:]
exr_GTEx = DataFrame([[names(ex_GTEx)]; collect.(eachrow(ex_GTEx))], [:column; Symbol.(axes(ex_GTEx,2))])
exr_GTEx[1,2:length(exr_GTEx[1,:])]
Array(exr_GTEx[1,2:length(exr_GTEx[1,:])])

###
names!(df, [Symbol("Col$i") for i in 1:size(df,2)])
