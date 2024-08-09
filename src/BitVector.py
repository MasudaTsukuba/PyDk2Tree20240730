"""
BitVector.py << BitVector.cpp
2024/8/1, T. Masuda
Amagasa Laboratory, University of Tsukuba
"""
# //
# // Created by hugo on 2-1-19.
# //

# #ifndef DK2TREE_BIT_VECTOR_H
# #define DK2TREE_BIT_VECTOR_H

# #include <vector>
# #include <cstdint>
# #include <cassert>
# #include <cstdio>
# #include "parameters.cpp"
from __future__ import annotations
import numpy as np
from src.Parameters import Parameters

# using std::vector;
# typedef uint64_t u64;
# typedef uint8_t u8;

MAX_BIT: np.uint64 = np.left_shift(np.uint64(1), np.uint64(63))
# #define MAX_BIT (((u64) 1) << 63)

# // This static table contains the number of 1-bits in each 8-bit integer
# // It is used to allow very fast computation of one counts in larger integers
ONE_BITS = [
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
]
# static const u8 ONE_BITS[256] = {
#         0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
#         1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
#         1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
#         2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
#         1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
#         2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
#         2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
#         3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
#         1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
#         2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
#         2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
#         3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
#         2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
#         3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
#         3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
#         4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
# };


def ones(n: np.uint64):
    pass
    # // Efficiently counts the number of 1-bits in a 64-bit integer using the table
    # // defined above.
    # // u64 ones(u64 n) {
    xff = np.uint64(0xFF)
    return ONE_BITS[n & xff] \
        + ONE_BITS[(n >> np.uint64(8)) & xff] \
        + ONE_BITS[(n >> np.uint64(16)) & xff] \
        + ONE_BITS[(n >> np.uint64(24)) & xff] \
        + ONE_BITS[(n >> np.uint64(32)) & xff] \
        + ONE_BITS[(n >> np.uint64(40)) & xff] \
        + ONE_BITS[(n >> np.uint64(48)) & xff] \
        + ONE_BITS[(n >> np.uint64(56)) & xff]
    # //     return ONE_BITS[n & 0xFF]
    # //            + ONE_BITS[(n >> 8) & 0xFF]
    # //            + ONE_BITS[(n >> 16) & 0xFF]
    # //            + ONE_BITS[(n >> 24) & 0xFF]
    # //            + ONE_BITS[(n >> 32) & 0xFF]
    # //            + ONE_BITS[(n >> 40) & 0xFF]
    # //            + ONE_BITS[(n >> 48) & 0xFF]
    # //            + ONE_BITS[(n >> 56) & 0xFF];
    # // }
# u64 ones(u64 n);  ////


LENGTH: int = int((Parameters.B + Parameters.BLOCK_SIZE + 63) / 64)  # (252 + 4 + 63) / 64 = 319 / 64 = 4
# template<unsigned long LENGTH = (B + BLOCK_SIZE + 63) / 64>


class BitVector:
    pass
    # /**
    #  * A simple bitvector containing the `raw` bits in a vector<bool>, as well as
    #  * a list of the number of ones in each block, to speed up rank operations
    #  */
    # struct BitVector {

    def __init__(self):
        pass
        # /* The size of this bitvector in bits */
        self.bits: int = 0
        # u64 bits;

        # /* The array of integers representing the bitvector */
        # self.data = [0 for i in range(LENGTH)]
        xxx = np.zeros(LENGTH, dtype=np.uint64)
        self.data: np.ndarray[np.uint64] = np.zeros(LENGTH, dtype=np.uint64)
        # u64 data[LENGTH];

        # /*
        #  * For each index `i`, block_counts[i] stores the number of 1-bits in data[i]
        #  */
        self.block_counts: list[np.uint8] = [np.uint8(0) for _ in range(LENGTH)]
        #  u8 block_counts[LENGTH];
        pass

    def __getitem__(self, n):
        pass
        #     /**
        #      * Gives the value of the n-th bit in the bitvector. This is a read-only operator,
        #      * since the block_counts must also be updated when writing
        #      *
        #      * @param n an index with 0 <= n < bv.size()
        #      * @return the value of the n-th bit of the bitvector
        #      */

        #   const bool operator[](unsigned long n) const {
        idx: np.uint64 = np.uint64(int(n / 64))
        #  u64 idx = n / 64;
        mask: np.uint64 = np.right_shift(MAX_BIT, np.uint64(n % 64))
        # mask = MAX_BIT >> (n % 64)
        # u64 mask = MAX_BIT >> (n % 64);
        xxx: np.uint64 = np.bitwise_and(self.data[idx], mask)
        return xxx != np.uint64(0)
        # return (self.data[idx] & mask) != 0
        # return (data[idx] & mask) != 0;
        # }
        pass

    def set(self, n: int, b: bool) -> bool:
        pass
        #   /**
        #    * Sets the n-th bit to value b, and returns true if the value changed
        #    * @param n an index with 0 <= n < size()
        #    * @param b a boolean
        #    * @return true iff the previous value of bit n was unequal to b
        #    */
        # const bool set(unsigned long n, bool b) {
        block: np.uint64 = np.uint64(int(n / 64))  # 0~63 -> block 0, 64~127 -> block 1
        # u64 block = n / 64;
        mask: np.uint64 = np.right_shift(MAX_BIT, np.uint64(n % 64))
        # mask = MAX_BIT >> (n % 64)
        # u64 mask = MAX_BIT >> (n % 64);

        xxx: np.uint64 = np.bitwise_and(self.data[block], mask)
        changed: bool = (xxx != np.uint64(0)) ^ b
        # changed = ((self.data[block] & mask) != 0) ^ b
        # bool changed = ((data[block] & mask) != 0) ^b;
        if changed:
            # if (changed) {
            if b:  # changed to True
                # if (b) {
                self.data[block] = np.bitwise_or(self.data[block], mask)  # change 0 to 1
                # self.data[block] |= mask
                # data[block] |= mask;
                self.block_counts[block] += 1  # Number of 1's in this 64-bit block
                # block_counts[block]++;
            else:
                # } else {
                self.data[block] = np.bitwise_and(self.data[block], np.bitwise_not(mask))  # change to 0
                # self.data[block] &= ~mask
                # data[block] &= ~mask;
                self.block_counts[block] -= 1  # Number of 1's in this 64-bit block
                # block_counts[block]--;
            #  }
        # }
        return changed
        # return changed;
        # }
        pass

    def rank1(self, n: int) -> int:
        pass
        #   /**
        #    * Performs the rank-operation on this bitvector
        #    * @param n an index with 0 <= n <= size()
        #    * @return the number of ones in the bits [0 ... n)
        #    */
        # unsigned long rank1(unsigned long n) {
        # // First split the interval [0, n) up into a whole number of blocks and
        # // a remainder, then count the total number of bits
        end_blocks = n - n % 64
        #       unsigned long end_blocks = n - n % 64;
        nr_blocks = int(end_blocks / 64)
        #       unsigned long nr_blocks = end_blocks / 64;
        return self.count_blocks(0, nr_blocks) + self.count_ones_raw(end_blocks, n)
        #       return countBlocks(0, nr_blocks) + countOnesRaw(end_blocks, n);
        #   }
        pass

    def range_rank1(self, lo, hi):
        pass
        #     /**
        #      * Returns the number of 1-bits in the interval [lo, hi), which is equal to
        #      * rank1(hi) - rank1(lo)
        #      */
        #       unsigned long rangeRank1(unsigned long lo, unsigned long hi) {
        block_lo = int((lo + 64 - 1) / 64)
        block_hi = int(hi / 64)
        #       unsigned long blockLo = (lo + 64 - 1) / 64, blockHi = hi / 64;
        if block_lo > block_hi:
            pass
            #       if (blockLo > blockHi) {
            return self.count_ones_raw(lo, hi)
            #           return countOnesRaw(lo, hi);
            #       }
        block_start = block_lo * 64
        block_end = block_hi * 64
        #       unsigned long blockStart = blockLo * 64, blockEnd = blockHi * 64;
        return self.count_ones_raw(lo, block_start) + self.count_ones_raw(block_lo, block_hi) + \
            self.count_ones_raw(block_end, hi)
        #       return countOnesRaw(lo, blockStart) + countBlocks(blockLo, blockHi) +
        #              countOnesRaw(blockEnd, hi);
        #   }
        pass

    def insert(self, begin: int, size: int):
        pass
        #   /**
        #    * Inserts `size` 0-bits at position `begin`
        #    * @param begin an index with 0 <= begin <= size()
        #    * @param size the number of bits to be inserted
        #    */
        #   void insert(unsigned long begin, unsigned long size) {
        block_start: np.uint64 = np.uint64(int(begin / 64))
        #       u64 block_start = begin / 64;
        block_amount: np.uint64 = np.uint64(int(size / 64))
        #       u64 block_amount = size / 64;
        bit_amount: np.uint64 = np.uint64(size % 64)
        #       u64 bit_amount = size % 64;
        #       // We save the first block, so we can set everything but the part to be
        #       // moved to zero, simplifying the rest
        # first_part_mask = (2 << (63 - begin % 64)) - 1
        xxx: np.uint64 = np.uint64(63) - np.uint64(begin % 64)
        yyy: np.uint64 = np.left_shift(np.uint(2), xxx)
        first_part_mask: np.uint64 = np.uint64(np.int64(yyy) - np.int64(1))  # ###
        #       u64 first_part_mask = (2ULL << (63 - begin % 64)) - 1;
        # first_block_keep = self.data[block_start] & ~first_part_mask
        first_block_keep: np.uint64 = np.bitwise_and(self.data[block_start], np.bitwise_not(first_part_mask))  # ###
        #       u64 first_block_keep = data[block_start] & ~first_part_mask;
        # self.data[block_start] &= first_part_mask
        self.data[block_start] = np.bitwise_and(self.data[block_start], first_part_mask)
        #       data[block_start] &= first_part_mask;
        #       // First, shift by whole number of blocks if applicable
        #       // The `if` is necessary since data would be destroyed otherwise
        if block_amount != 0:
            #       if (block_amount != 0) {
            pass
            for idx in range(LENGTH - 1, block_start + block_amount - 1, -1):
                pass
                #           for (u64 idx = LENGTH - 1;
                #               idx >= block_start + block_amount; idx--) {
                self.data[idx] = self.data[idx - block_amount]
                #               data[idx] = data[idx - block_amount];
                self.data[idx - block_amount] = np.uint64(0)
                #               data[idx - block_amount] = 0;
                #           }
                #       }
        if bit_amount != 0:
            pass
            #       // Then, shift by the remaining number of bits if applicable
            #       // The `if` is necessary since the code would otherwise perform bit-
            #       // shifts by 64 bits, which is undefined behaviour
            #       if (bit_amount != 0) {
            for idx in range(LENGTH - 1, block_start, -1):
                #           for (u64 idx = LENGTH - 1; idx >= block_start + 1; idx--) {
                pass
                self.data[idx] = np.bitwise_or(np.right_shift(self.data[idx], np.uint64(bit_amount)),
                                               np.left_shift(self.data[idx - 1], np.uint64(64 - bit_amount)))
                # self.data[idx] = (self.data[idx] >> bit_amount) | \
                #               (self.data[idx - 1] << (64 - bit_amount))
                #               data[idx] = (data[idx] >> bit_amount) |
                #                           (data[idx - 1] << (64 - bit_amount));
                #           }
            self.data[block_start] = np.right_shift(self.data[block_start], np.uint64(bit_amount))
            # self.data[block_start] >>= bit_amount
            #           data[block_start] >>= bit_amount;
            #       }
            #       // Finally, restore the first block and recompute the block counts
        self.data[block_start] = np.bitwise_or(self.data[block_start], first_block_keep)
        # self.data[block_start] |= first_block_keep
        #       data[block_start] |= first_block_keep;
        self.bits += size
        #       bits += size;
        self.recompute(begin)
        #       recompute(begin);
        #   }
        pass

    def insert_range(self, begin: int, from_: BitVector, lo: int, hi: int):
        pass
        #   /**
        #    * Inserts the range [lo, hi) of the bit vector `from` into this bit vector
        #    *
        #    * @param begin the position in this bit vector to insert into
        #    * @param from the bit vector to insert a range from
        #    * @param lo the start of the range in `from` to insert
        #    * @param hi the end of the range in `from` to insert
        #    */
        #   void
        #   insert(unsigned long begin, const BitVector<LENGTH> &from, unsigned long lo,
        #          unsigned long hi) {
        #       // This can probably be done faster, but this operation usually only
        #       // performed with [lo, hi) being a single k^2 block
        self.insert(begin, hi - lo)
        #       insert(begin, hi - lo);
        for idx in range(hi):
            #       for (u64 idx = 0; idx + lo < hi; idx++) {
            self.set(idx + begin, from_[idx + lo])
            # set(idx + begin, from[idx + lo]);
            #       }
        self.recompute(begin)
        #       recompute(begin);
        #   }
        pass

    def append(self, from_: BitVector, lo: int, hi: int):
        pass
        #   /**
        #    * Appends the range [lo, hi) of the bit vector `from` into this bit vector
        #    *
        #    * @param from the bit vector to append bits from
        #    * @param lo the start of the range of bits to append
        #    * @param hi the end of the range of bits to append
        #    */
        #   void append(const BitVector &from, unsigned long lo, unsigned long hi) {
        self.insert_range(self.bits, from_, lo, hi)
        #       insert(bits, from, lo, hi);
        #   }
        pass

    def erase(self, lo: int, hi: int):
        pass
        #   /**
        #    * Deletes the indicated range of bits starting at the indicated index
        #    * @param lo an index with 0 <= lo <= size()
        #    * @param hi the end of the range of bits to be deleted. Should satisfy lo <= hi <= size()
        #    */
        #   void erase(unsigned long lo, unsigned long hi) {
        amount: np.uint64 = np.uint64(hi - lo)
        # u64 amount = hi - lo;
        block_start: np.uint64 = np.uint64(int(lo / 64))
        #       u64 block_start = lo / 64;
        block_amount: np.uint64 = np.uint64(int(amount / 64))
        #       u64 block_amount = amount / 64;
        bit_amount: np.uint64 = np.uint64(amount % 64)
        #       u64 bit_amount = amount % 64;

        #       // We save the first block, so we can set everything but the part to be
        #       // moved to zero, simplifying the rest
        xxx: np.uint64 = np.left_shift(np.uint64(2), np.uint64((63 - lo % 64)))
        first_part_mask: np.uint64 = np.uint64(np.int64(xxx) - np.int64(1))
        #         first_part_mask = (2 << (63 - lo % 64)) - 1
        #       u64 first_part_mask = (2ULL << (63 - lo % 64)) - 1;
        first_block_keep: np.uint64 = np.bitwise_and(self.data[block_start], np.bitwise_not(first_part_mask))
        #         first_block_keep = self.data[block_start] & ~first_part_mask
        #       u64 first_block_keep = data[block_start] & ~first_part_mask;
        #         self.data[block_start] &= first_part_mask
        self.data[block_start] = np.bitwise_and(self.data[block_start], first_part_mask)
        #       data[block_start] &= first_part_mask;

        #       // First, move everything over by the specified number of blocks
        if block_amount != 0:
            #       if (block_amount != 0) {
            for idx in range(block_start, LENGTH - block_amount):
                # for (u64 idx = block_start; idx + block_amount < LENGTH; idx++) {
                self.data[idx] = self.data[idx + block_amount]
                # data[idx] = data[idx + block_amount];
                self.data[idx + block_amount] = np.uint64(0)
                #               data[idx + block_amount] = 0;
                #           }
                #       }

        #       // Then, shift everything over by the correct bit-amount
        if bit_amount != 0:
            #       if (bit_amount != 0) {
            for idx in range(block_start, LENGTH - 1):
                # for (u64 idx = block_start; idx + 1 < LENGTH; idx++) {
                # self.data[idx] = (self.data[idx] << bit_amount) | (self.data[idx+1] >> (64 - bit_amount))
                self.data[idx] = np.bitwise_or(np.left_shift(self.data[idx], np.uint64(bit_amount)),
                                               np.right_shift(self.data[idx+1], np.uint64(64 - bit_amount)))
                # data[idx] = (data[idx] << bit_amount) |
                #             (data[idx + 1] >> (64 - bit_amount));
                #           }
                #             self.data[LENGTH - 1] <<= bit_amount
            self.data[LENGTH - 1] = np.left_shift(self.data[LENGTH - 1], np.uint64(bit_amount))
            #           data[LENGTH - 1] <<= bit_amount;
            #       }
            #       // Finally, fix the first block and recompute
        self.data[block_start] = np.bitwise_and(self.data[block_start], first_part_mask)
        #         self.data[block_start] &= first_part_mask
        #       data[block_start] &= first_part_mask;
        self.data[block_start] = np.bitwise_or(self.data[block_start], first_block_keep)
        #         self.data[block_start] |= first_block_keep
        #       data[block_start] |= first_block_keep;
        self.bits -= amount
        #       bits -= amount;
        self.recompute(lo)
        #       recompute(lo);
        #   }
        pass

    def size(self) -> int:
        pass
        #     /**
        #      * Gets the size (number of bits) of this bit vector
        #      */
        #       unsigned long size() {
        return self.bits
        #       return bits;
        #   }
        pass

    def init_with_size(self, size):
        pass
        #     /**
        #      * Construct an all-zeros bit vector with the given number of bits
        #      *
        #      * @param size the number of bits the constructed vector will contain
        #      */
        #     explicit BitVector(unsigned long size) :
        # self.__init__()
        self.bits = size
        #       bits(size),
        #       data{0},
        #       block_counts{0} {}
        return self
        pass

    def init_with_range(self, from_, lo, hi):
        pass
        #     /**
        #      * Constructs a bit vector from the range [lo, hi) of another bit vector
        #      *
        #      * @param from the bit vector to copy a range of bits from
        #      * @param lo the start of the range of bits to take
        #      * @param hi the end of the range of bits to take
        #      */
        #   BitVector(const BitVector<LENGTH> &from, unsigned long lo, unsigned long hi) :
        # self.__init__()
        self.bits = from_.bits
        #       bits(from.bits),
        #       data{0},
        #       block_counts{0} {
        for idx in range(LENGTH):
            #       for (u8 idx = 0; idx < LENGTH; idx++) {
            pass
            self.data[idx] = from_.data[idx]
            #           data[idx] = from.data[idx];
            self.block_counts[idx] = from_.block_counts[idx]
            #           block_counts[idx] = from.block_counts[idx];
            #       }
        self.erase(hi, self.bits)
        #       erase(hi, bits);
        self.erase(0, lo)
        #       erase(0, lo);
        return self
        #   }
        pass

    def memory_usage(self) -> int:
        #   unsigned long memoryUsage() {
        pass
        #       // For each 64-bit block, we store a 64-bit integer (containing those
        #       // bits), and an 8-bit integer storing the number of one-bits
        #       // So 72 bits = 9 bytes for each block
        return int(int(self.bits + 63) / 64) * 9
        #       return ((bits + 63) / 64) * 9;
        #   }
        pass

    def recompute(self, start=0):
        pass
        # private:
        #     /**
        #      * Private method to re-compute all the values of block_counts from a certain
        #      * starting point. Used when inserting or deleting bits
        #      * @param start the first bit that may have changed and require updating the counters
        #      */
        #     void recompute(unsigned long start = 0) {
        start = int(start / 64)
        #       start /= 64;
        for block in range(start, LENGTH):
            #       for (u64 block = start; block < LENGTH; block++) {
            pass
            xxx: np.uint64 = np.uint64(self.data[block])
            self.block_counts[block] = ones(xxx)
            #           block_counts[block] = (u8) ones(data[block]);
            #       }
            #   }
        pass

    def count_ones_raw(self, lo: int, hi: int) -> int:
        pass
        #     /**
        #      * Counts the number of 1-bits in the interval [lo, hi). This interval
        #      * should be entirely in one b64-bit lock
        #      *
        #      * @param lo the start of the interval
        #      * @param hi the end of the interval
        #      * @return the number of 1-bits in the interval [lo, hi)
        #      */
        #   unsigned long countOnesRaw(unsigned long lo, unsigned long hi) {
        block: np.uint64 = np.uint64(int(lo / 64))
        #       u64 block = lo / 64;
        lo -= block * 64
        #       lo -= block * 64;
        hi -= block * 64
        #       hi -= block * 64;

        aaa: np.uint64 = np.left_shift(np.uint64(1), np.uint64(hi - lo)) - np.uint64(1)
        mask: np.uint64 = np.left_shift(aaa, np.uint64(64 - hi))
        # mask = ((1 << (hi - lo)) - 1) << (64 - hi)
        #       u64 mask = ((1ULL << (hi - lo)) - 1) << (64 - hi);
        # xxx = self.data[block] & mask
        xxx = np.bitwise_and(self.data[block], mask)
        return ones(xxx)
        #       return ones(data[block] & mask);
        #   }
        pass

    def count_blocks(self, lo: int, hi: int) -> int:
        pass
        #     /**
        #      * Counts the total number of ones in the blocks in interval [lo, hi)
        #      * @param lo the start of the interval
        #      * @param hi the end of the interval
        #      * @return the number of 1-bits in the interval [lo, hi) of blocks
        #      */
        #     unsigned long countBlocks(unsigned long lo, unsigned long hi) {
        tot: int = 0
        #       unsigned long tot = 0;
        for k in range(lo, hi):
            #       for (unsigned long k = lo; k < hi; k++) {
            pass
            tot += self.block_counts[k]
            #           tot += block_counts[k];
            #       }
        return tot
        #       return tot;
        #   }
        pass
        # };  // end of struct BitVector
    pass

# #endif //DK2TREE_BIT_VECTOR_H

# """
# BitVector.cpp
# """
# #include "BitVector.h"

# u64 ones(u64 n) {
#     return ONE_BITS[n & 0xFF]
#            + ONE_BITS[(n >> 8) & 0xFF]
#            + ONE_BITS[(n >> 16) & 0xFF]
#            + ONE_BITS[(n >> 24) & 0xFF]
#            + ONE_BITS[(n >> 32) & 0xFF]
#            + ONE_BITS[(n >> 40) & 0xFF]
#            + ONE_BITS[(n >> 48) & 0xFF]
#            + ONE_BITS[(n >> 56) & 0xFF];
# }
