"""
Parameters.py << Parameters.cpp
2024/8/1, T. Masuda
Amagasa Laboratory, University of Tsukuba
"""
# //
# // Created by hugo on 19-2-19.
# //
#
# #ifndef DKTREE_PARAMETERS
# #define DKTREE_PARAMETERS


class Parameters:
    pass
    # /// The three main parameters for the TTree and LTree representation
    k: int = 2
    # static const unsigned int k = 2; // The `k` in the k2-tree
    BLOCK_SIZE: int = k * k  # 4
    # static const unsigned int BLOCK_SIZE =  k * k; // The number of bits in one block of the bit vector
    B: int = 256 - BLOCK_SIZE  # 252
    # static const unsigned int B = 256 - BLOCK_SIZE; // The maximum size (in bits) of a leaf bitvector

    # /// The maximum/minimum number of children/blocks an internal node/leaf node
    # /// is allowed to have, as per the rules of the B+tree
    nodeSizeMax: int = 3
    # static const unsigned int nodeSizeMax = 3;
    nodeSizeMin: int = int((nodeSizeMax + 1) / 2)  # 2
    # static const unsigned int nodeSizeMin = (nodeSizeMax + 1) / 2;
    leafSizeMax: int = int(B / BLOCK_SIZE)  # 252 /4 = 63
    # static const unsigned int leafSizeMax = B / BLOCK_SIZE;
    leafSizeMin: int = int((leafSizeMax + 1) / 2)  # (63 + 1) / 2 = 32
    # static const unsigned int leafSizeMin = (leafSizeMax + 1) / 2;

# #endif // DKTREE_PARAMETERS
