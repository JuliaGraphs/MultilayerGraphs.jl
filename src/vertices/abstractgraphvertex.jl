abstract type AbstractGraphVertex{G <: AbstractGraph} end

struct LayerVertex{L} <: AbstractGraphVertex{G} end
