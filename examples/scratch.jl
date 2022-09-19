# ψ0 is my wavefunction
# n is some index in the MPS 
T = ψ0[n] # This is an ITensor
T̃ = T.tensor # This is a block-sparse tensor 
T̃_size = prod(size(T̃.storage)) # size of the block-sparse tensor
T_size = prod(size(T)) # size of the original tensor
@show T̃_size/T_size
