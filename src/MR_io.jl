#### Functions to make IO easier

### NETCDF stuff
using NCDatasets
using Chain

function ncreadall(file::AbstractString)
    #ncd = NCDataset(file)
    #return [Symbol(k)=>ncd[k][:] for k in keys(ncd)] |> NamedTuple

    @chain file begin
        NCDataset
        [Symbol(k)=>_[k][:] for k in keys(_)]
        NamedTuple
    end
end

