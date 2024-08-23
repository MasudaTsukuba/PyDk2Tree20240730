"""
TTree.py << TTree.cpp
2024/8/1, T. Masuda
Amagasa Laboratory, University of Tsukuba
"""
# //
# // Created by anneke on 18/12/18.
# //

# #include <iostream>
# #include <utility>
# #include "TTree.h"
from __future__ import annotations
import sys
from src.Parameters import Parameters
from src.BitVector import BitVector

# --------------------------------------------------------------------
# //
# // Created by anneke on 18/12/18.
# //

# #ifndef DK2TREE_TTREE_H
# #define DK2TREE_TTREE_H

# #include "BitVector.h"
# #include <utility>
# #include "parameters.cpp"


class Record:
    pass
    # /// Record type containing the number pf preceding bits and ones, and the
    # /// index of a child node in the parent's `entries` list
    # struct Record {

    def __init__(self, b: int, o: int, i: int):
        pass
        self.b = b
        #       unsigned long b;
        self.o = o
        #       unsigned long o;
        self.i = i
        #       unsigned long i;
        # };
    pass  # end of class Record

# /// The three main structs forming the tree
# /// `TTree` represents one node in the tree, which contains a pointer to its
# /// parent (or nullptr for the root) and the index in the parent,
# /// as well as either an InternalNode containing entries (b, o, P) or a LeafNode
# /// containing a BitVector
# struct InternalNode;
# struct LeafNode;
# struct TTree;


class Nesbo:
    pass
    # struct Nesbo {
    # public:
    #     TTree *node;
    #     unsigned long index;
    #     unsigned long size;
    #     unsigned long bitsBefore;
    #     unsigned long onesBefore;

    def __init__(self, node, index: int, size: int, bitsBefore: int, onesBefore: int):
        pass
        #       Nesbo(TTree *node, unsigned long index, unsigned long size,
        #           unsigned long bitsBefore, unsigned long onesBefore) :
        self.node: TTree = node
        #             node(node),
        self.index: int = index
        #             index(index),
        self.size: int = size
        #             size(size),
        self.bitsBefore: int = bitsBefore
        #             bitsBefore(bitsBefore),
        self.onesBefore: int = onesBefore
        #             onesBefore(onesBefore) {}
    # };
    pass  # end of class Nesbo


class InternalNode:
    pass

    # /**
    #  * A struct for the internal nodes of the tree
    #  * This contains a list of Entries of the form <b, o, P>
    #  */
    class Entry:
        pass
        # struct InternalNode {
        #     struct Entry {
        #         unsigned long b;
        #         unsigned long o;
        #         TTree *P;

        def __init__(self):
            pass
            #         Entry() :
            #                 b(0), o(0), P(nullptr) {}
            self.b = 0
            self.o = 0
            self.p = None
            pass

        def init_with_p(self, p):
            pass
            #         explicit Entry(TTree *);
            # # InternalNode::Entry::Entry(TTree *P) :
            self.b = p.bits()
            # #         b(P->bits()),
            self.o = p.ones()
            # #         o(P->ones()),
            self.p = p
            # #         P(P) {}
            return self
            pass

        def init_with_bop(self, b, o, p):
            pass
            #         /// Construct Entry from the three fields
            #         /// This method is not unused, but is used as an initializer list
            #         Entry(unsigned long b, unsigned long o, TTree *P) :
            self.b = b
            self.o = o
            self.p = p
            #                 b(b), o(o), P(P) {}
            return self
            pass

        def remove(self):
            pass
            #         /**
            #          * This function is called in the destructor of `InternalNode`, to
            #          * delete the child nodes. It is not part of the destructor of `Entry`,
            #          * since that causes problems when the struct is used in methods such as
            #          * `findLeaf`
            #          */
            #         void remove();
            # # void InternalNode::Entry::remove() {
            del self.p
            # #     delete P;
            # # }
            pass
        # };
        pass  # end of class Entry

    def __init__(self):
        #     /// The default constructor creates an empty internal node
        #     InternalNode() :
        #             size(0),
        #             entries{Entry()} {}
        self.size: int = 0
        #     /// The number of children this node has
        #     unsigned long size;
        self.entries: list[InternalNode.Entry] = [InternalNode.Entry() for _ in range(Parameters.nodeSizeMax + 1)]
        #     /// An array of pointers to the child nodes
        #     /// This is one more than the maximum, so that we can split nodes
        #     /// after insertion instead of before
        #     Entry entries[nodeSizeMax + 1];

    # def init_with_arguments(self, left: TTree, right: TTree, parent=None):
    def init_with_arguments(self, left: TTree, right: TTree):
        pass
        #     /**
        #      * Creates a new internal node with the given two children
        #      *
        #      * @param left the first child of this node
        #      * @param right the second child of this node
        #      * @param parent the parent node, which has this as its internal node
        #      *        the left and right TTrees have their parent and indexInParent
        #      *        set correctly as well
        #      */
        #     InternalNode(TTree *left, TTree *right, TTree *parent = nullptr);
        # # InternalNode::InternalNode(TTree *left, TTree *right, TTree *parent) :
        self.size = 2
        # #         size(2),
        # left.parent = parent  # parent is set at new TTree(left, right)  # 2024/8/22
        # #     left->parent = parent;
        left.indexInParent = 0
        # #     left->indexInParent = 0;
        # right.parent = parent
        # #     right->parent = parent;
        right.indexInParent = 1
        # #     right->indexInParent = 1;
        self.entries = [self.Entry().init_with_p(left), self.Entry().init_with_p(right), self.Entry(), self.Entry()]
        # #         entries{Entry(left), Entry(right), Entry()} {
        return self
        # # }
        pass

    def bits(self):
        pass
        #     /**
        #      * When an internal node is dropped, clear the entries it points to
        #      */
        #     ~InternalNode() {
        #         for (auto &entry : entries) {
        #             entry.remove();
        #         }
        #     }

        #     /**
        #      * Returns the total number of bits in this tree, by summing up the
        #      * b-parts of the entries
        #      */
        #     unsigned long bits();
        # # unsigned long InternalNode::bits() {
        total = 0
        # #     unsigned long total = 0;
        for entry in self.entries:
            pass
        # #     for (const auto &entry : entries) {
            total += entry.b
        # #         total += entry.b;
        # #     }
        return total
        # #     return total;
        # # }
        pass

    def ones(self):
        pass
        #     /**
        #      * Returns the total number of ones in this tree, by summing up the
        #      * o-parts of the entries
        #      */
        #     unsigned long ones();
        # # unsigned long InternalNode::ones() {
        total = 0
        # #     unsigned long total = 0;
        for entry in self.entries:
            pass
            # #     for (const auto &entry : entries) {
            total += entry.o
            # #         total += entry.o;
            # #     }
        return total
        # #     return total;
        # # }
        pass

    def pop_first(self):
        pass
        #     /**
        #      * Takes the leftmost entry out of this node and returns it
        #      * @return
        #      */
        #     Entry popFirst();
        # # InternalNode::Entry InternalNode::popFirst() {
        # #     // Take the first entry out, move everything else left
        result = self.entries[0]
        # #     InternalNode::Entry result = this->entries[0];
        self.remove(0)
        # #     this->remove(0);
        return result
        # #     return result;
        # # }
        pass

    def pop_last(self):
        pass
        #     /**
        #      * Takes the rightmost entry out of this node and returns it
        #      * @return
        #      */
        #     Entry popLast();
        # # InternalNode::Entry InternalNode::popLast() {
        # #     // Take the last entry out
        self.size -= 1
        # #     size--;
        result = self.entries[self.size]
        # #     InternalNode::Entry result = entries[size];
        self.entries[self.size] = self.Entry()
        # #     entries[size] = Entry();
        return result
        # #     return result;
        # # }
        pass

    def insert(self, idx: int, entry: InternalNode.Entry):
        pass
        #     /**
        #      * Adds the given entry to this node at the specified position
        #      */
        #     void insert(unsigned long, Entry);
        # # void InternalNode::insert(unsigned long idx, InternalNode::Entry entry) {
        # #     // Move everything from idx onwards right
        i = self.size
        while i > idx:
            pass
            # #     for (unsigned long i = size; i > idx; i--) {
            self.entries[i] = self.entries[i - 1]
            # #         entries[i] = entries[i - 1];
            self.entries[i].p.indexInParent = i
            if i >= 3:
                pass  # debug  # 2024/8/23
            # #         entries[i].P->indexInParent = i;
            i -= 1
            # #     }
        self.entries[idx] = entry
        # #     entries[idx] = entry;
        self.entries[idx].p.indexInParent = idx
        # #     entries[idx].P->indexInParent = idx;
        self.size += 1
        # #     size++;
        # # }
        pass

    def append(self, entry):
        pass
        #     /**
        #      * Adds the given entry to the end of this node
        #      */
        #     void append(Entry);
        # # void InternalNode::append(InternalNode::Entry entry) {
        entry.p.indexInParent = self.size
        if self.size >= 3:
            pass  # debug  # 2024/8/23
        # #     entry.P->indexInParent = size;
        self.entries[self.size] = entry
        # #     entries[size] = entry;
        self.size += 1
        # #     size++;
        # # }
        pass

    def remove(self, idx: int):
        pass
        #     /**
        #      * Removes the entry at the specified position from this node
        #      */
        #     void remove(unsigned long);
        # # void InternalNode::remove(unsigned long idx) {
        self.size -= 1
        # #     size--;
        for i in range(idx, self.size):
            pass
            # #     for (unsigned long i = idx; i < size; i++) {
            self.entries[i] = self.entries[i + 1]
            # #         entries[i] = entries[i + 1];
            self.entries[i].p.indexInParent = i
            if i >= 3:
                pass  # debug  # 2024/8/23
            # #         entries[i].P->indexInParent = i;
            # #     }
        self.entries[self.size] = self.Entry()
        # #     entries[size] = Entry();
        # # }
        pass
    # };
    pass  # end of class InternalNode


class LeafNode:
    pass
    # /** A leaf node, which consists of a bitvector (represented by vector<bool>) */
    # struct LeafNode {

    def __init__(self):
        self.bv = BitVector()
        #       BitVector<> bv;
        pass

    def init_with_size(self, size: int):
        pass
        #     /**
        #      * Constructs a leaf with the given number of bits
        #      */
        #       explicit LeafNode(unsigned long size) :
        self.bv = BitVector().init_with_size(size)
        #       bv(size) {}
        return self
        pass

    def init_with_bv(self, bv):
        pass
        #     /**
        #      * Constructs a leaf node from the given bit vector
        #      * @param bv the bitvector to be moved into this leaf node
        #      */
        #       explicit LeafNode(BitVector<> bv) :
        self.bv = bv
        #       bv(bv) {}
        return self
        pass

    def bits(self):
        pass
        #     /**
        #      * Get the total number of bits stored in this leaf node
        #      * @return the size in bits of this leaf
        #      */
        #     unsigned long bits();
        # # unsigned long LeafNode::bits() {
        return self.bv.size()
        # #     return bv.size();
        # # }
        pass

    def ones(self):
        pass
        #     /**
        #      * Get the number of ones in this leaf node
        #      * @return the number of bits in this leaf that are set to 1
        #      */
        #     unsigned long ones();
        # # unsigned long LeafNode::ones() {
        # #     // Count all ones manually
        return self.bv.rank1(self.bv.size())
        # #     return bv.rank1(bv.size());
        # # }
        pass

    # };
    pass  # end of class LeafNode


class TTree:
    pass
    # /** TTree is the struct representing a single node (leaf or internal) of the TTree */
    # /** A single node is either an internal or leaf node, as indicated by the isLeaf value */
    # struct TTree {

    class Node:
        pass

        def __init__(self):
            pass
            # # TTree::Node::Node() {
            self.internalNode = None
            # #     this->internalNode = nullptr;
            self.leafNode = LeafNode().init_with_size(0)
            # #     this->leafNode = new LeafNode(0);
            self.isLeaf = False
            #     bool isLeaf;
            self.parent = None
            #     TTree *parent = nullptr;
            self.indexInParent = 0
            #     unsigned long indexInParent = 0;
            # # }
            pass

        def init_with_ps(self, p1: TTree, p2: TTree):
            pass
            # # TTree::Node::Node(TTree *P1, TTree *P2) {
            self.leafNode = None
            # #     this->leafNode = nullptr;
            self.internalNode = InternalNode().init_with_arguments(p1, p2)
            # #     this->internalNode = new InternalNode(P1, P2);
            return self
            # # }
            pass

        def init_with_bv(self, bv):
            pass
            # # TTree::Node::Node(BitVector<> bv) {
            self.internalNode = None
            # #     this->internalNode = nullptr;
            self.leafNode = LeafNode().init_with_bv(bv)
            # #     this->leafNode = new LeafNode(bv);
            return self
            # # }
            pass

        def init_with_size(self, size: int):
            pass
            # # TTree::Node::Node(unsigned long size) {
            self.internalNode = None
            # #     this->internalNode = nullptr;
            self.leafNode = LeafNode().init_with_size(size)
            # #     this->leafNode = new LeafNode(size);
            return self
            # # }
            pass

        #     union Node {
        #         InternalNode *internalNode;
        #         LeafNode *leafNode;
        #
        #         Node();
        #
        #         Node(TTree *, TTree *);
        #
        #         explicit Node(BitVector<>);
        #
        #         explicit Node(unsigned long);
        #     } node;
        pass  # end of class Node

    def __init__(self):
        pass
        #     /**
        #      * Constructs an empty leaf node
        #      */
        #     TTree() :
        self.isLeaf = True
        #             isLeaf(true),
        self.node = TTree.Node()
        #             node() {}
        self.parent = None
        self.indexInParent: int = 0  # 2024/8/9
        pass

    def init_with_arguments(self, left: TTree, right: TTree):
        pass
        #     /**
        #      * Constructs a node with the two given `TTree`s as children
        #      * @param left the first child of this node
        #      * @param right the second child of this node
        #      */
        #       TTree(TTree *left, TTree *right) :
        self.isLeaf = False
        #             isLeaf(false),
        left.parent = self
        #         left->parent = this;
        left.indexInParent = 0
        #         left->indexInParent = 0;
        right.parent = self
        #         right->parent = this;
        right.indexInParent = 1
        #         right->indexInParent = 1;
        self.node = TTree.Node().init_with_ps(left, right)
        #             node(left, right) {
        return self
        #     }
        pass

    def init_with_bv(self, bv):
        pass
        #     /**
        #      * Constructs a leaf node with the given bit vector
        #      * @param bv the bit vector to be moved into this leaf node
        #      */
        #     explicit TTree(BitVector<> bv) :
        self.isLeaf = True
        #             isLeaf(true),
        self.node = TTree.Node().init_with_bv(bv)
        #             node(bv) {}
        return self
        pass

    def init_with_size(self, size: int):
        pass
        #     /**
        #      * Constructs an all-zeros leaf node with the specified size
        #      * @param size the size of this leaf node in bits
        #      */
        #     explicit TTree(unsigned long size) :
        self.isLeaf = True
        #             isLeaf(true),
        self.node = TTree.Node().init_with_size(size)
        #             node(size) {}
        return self
        pass

    def __del__(self):
        #     /// The TTree destructor decides which variant of the union to destroy
        #     ~TTree();
        # # TTree::~TTree() {
        if self.isLeaf:
            # #     if (isLeaf) {
            del self.node.leafNode
            # #         delete node.leafNode;
        else:
            # #     } else {
            del self.node.internalNode
            # #         delete node.internalNode;
            # #     }
        # # }
        pass

    def depth(self):
        pass
        #     /**
        #      * Returns the depth of this node, which is the length of the path from this
        #      * node to the root
        #      * @return the depth of `this`
        #      */
        #     unsigned long depth();
        # # unsigned long TTree::depth() {
        if self.parent is None:
            pass
            # #     if (parent == nullptr) {
            return 0
            # #         return 0;
        else:
            pass
            # #     } else {
            return 1 + self.parent.depth()
            # #         return 1 + parent->depth();
            # #     }
            # # }
        pass

    def height(self):
        pass
        #     /**
        #      * Returns the height of this node, which is the length of the longest path from
        #      * this node to any leaf in its subtree.
        #      * @return the height of `this`
        #      */
        #     unsigned long height();
        # # unsigned long TTree::height() {
        if self.isLeaf:
            pass
            # #     if (isLeaf) {
            return 0
            # #         return 0;
        else:
            pass
            # #     } else {
            entries = self.node.internalNode.entries
            # #         auto &entries = node.internalNode->entries;
            max_ = 0
            # #         unsigned long max = 0;
            for entry in entries:
                pass
                # #         for (auto &entry : entries) {
                depth = entry.p.height()
                # #             unsigned long depth = entry.P->height();
                if depth > max_:
                    pass
                    # #             if (depth > max) {
                    max_ = depth
                    # #                 max = depth;
                # #             }
            # #         }
            return max_ + 1
            # #         return max + 1;
            # #     }
        # # }
        pass

    def size(self) -> int:
        pass
        #     /**
        #      * Gets the size of this node, either as the number of children (internal
        #      * node) or in number of k^2 blocks (leaf node)
        #      * @return the size of this node in number of children or bit blocks,
        #      *         depending on the type of the node
        #      */
        #     unsigned long size();
        # # unsigned long TTree::size() {
        if self.isLeaf:
            pass
            # #     if (isLeaf) {
            return int(self.node.leafNode.bits() / Parameters.BLOCK_SIZE)
            # #         return node.leafNode->bits() / BLOCK_SIZE;
        else:
            pass
            # #     } else {
            return self.node.internalNode.size
            # #         return node.internalNode->size;
            # #     }
            # # }
        pass

    def find_child(self, n: int) -> Record:
        pass
        #     /**
        #      * Given an integer n, returns the child node containing the n-th bit in this subtree,
        #      *      as well as the numbers of bits and ones preceding it.
        #      * @param n  an integer that satisfies 0 <= n < (number of bits in this subtree)
        #      * @return  a TTree::InternalNode::Entry with:
        #      *      b = number of bits preceding this node in the TTree
        #      *      o = number of ones preceding this node in the TTree
        #      *      i = the index of the relevant subtree in the `entries` of `this`
        #      */
        #     Record findChild(unsigned long);
        # # Record TTree::findChild(unsigned long n) {
        bitsBefore = 0
        # #     unsigned long bitsBefore = 0;
        onesBefore = 0
        # #     unsigned long onesBefore = 0;
        node = self.node.internalNode
        # #     InternalNode *node = this->node.internalNode;
        # #     unsigned long i;
        i = 0
        sss = node.size
        for ii in range(sss):
            i = ii
            pass
            # #     for (i = 0; i < node->size; i++) {
            entry = node.entries[ii]
            # #         auto &entry = node->entries[i];
            if bitsBefore + entry.b > n:
                pass
                # #         if (bitsBefore + entry.b > n) {
                return Record(bitsBefore, onesBefore, ii)
                # #             return {bitsBefore, onesBefore, i};
                # #         }
            bitsBefore += entry.b
            # #         bitsBefore += entry.b;
            onesBefore += entry.o
            # #         onesBefore += entry.o;
        # #     }
        i += 1
        # #     // If the required bit is one after the last bit in this tree,
        # #     // return the last child anyway
        # #     // This is necessary for appending bits
        if i == node.size and bitsBefore == n:
            pass
            # #     if (i == node->size && bitsBefore == n) {
            entry = node.entries[i - 1]
            # #         auto &entry = node->entries[i - 1];
            return Record(bitsBefore - entry.b, onesBefore - entry.o, i - 1)
            # #         return {bitsBefore - entry.b, onesBefore - entry.o, i - 1};
            # #     }
        # #     // If we reach this point, that means that the size of this subtree is less than n
        # #     // so the input parameter was out of range
        raise "TTree: index out of range"
        # #     throw std::range_error("TTree: index out of range");
        # # }
        pass

    def find_leaf(self, n: int, path: list[Nesbo]) -> InternalNode.Entry:
        pass
        #     /**
        #      * Given an integer n, returns the leaf containing the n-th bit of the bitvector,
        #      *      as well as the numbers of bits and ones preceding this leaf node.
        #      * @param n  an integer that satisfies 0 <= n < (number of bits in this subtree)
        #      * @return  a TTree::InternalNode::Entry with:
        #      *      b = number of bits preceding this node in the TTree
        #      *      o = number of ones preceding this node in the TTree
        #      *      P = a pointer to the leaf node containing the n-th bit
        #      * Note that P is strictly just a pointer to a `TTreeNode`, as defined by the `entry` struct
        #      * But the union is always of the `LeafNode` variant.
        #      */
        #     InternalNode::Entry findLeaf(unsigned long, vector<Nesbo> *path = nullptr);
        # # InternalNode::Entry TTree::findLeaf(unsigned long n, vector<Nesbo> *path) {
        if path is None:  # or len(path) == 0:
            pass
            # #     if (path == nullptr) {
            current = self
            # #         auto *current = this;
            bitsBefore: int = 0
            # #         unsigned long bitsBefore = 0;
            onesBefore: int = 0
            # #         unsigned long onesBefore = 0;
            while not current.isLeaf:
                pass
                # #         while (!current->isLeaf) {
                record = current.find_child(n - bitsBefore)
                # #             auto record = current->findChild(n - bitsBefore);
                bitsBefore += record.b
                # #             bitsBefore += record.b;
                onesBefore += record.o
                # #             onesBefore += record.o;
                current = current.node.internalNode.entries[record.i].p
                # #             current = current->node.internalNode->entries[record.i].P;
                # #         }
            return InternalNode.Entry().init_with_bop(bitsBefore, onesBefore, current)
            # #         return {bitsBefore, onesBefore, current};
        else:
            pass
            # #     } else {
            return self.find_leaf2(n, path)
            # #         return findLeaf2(n, *path);
        # #     }
        # # }
        pass

    def find_leaf2(self, n: int, path: list[Nesbo]) -> InternalNode.Entry:
        pass
        #     InternalNode::Entry findLeaf2(unsigned long, vector<Nesbo> &);
        # # InternalNode::Entry TTree::findLeaf2(unsigned long n, vector<Nesbo> &path) {
        current = None
        # #     TTree *current = nullptr;
        bitsBefore: int = 0
        onesBefore: int = 0
        # #     unsigned long bitsBefore = 0, onesBefore = 0;
        if path is None or len(path) == 0:
            pass
            # #     if (path.empty()) {
            # #         // If the path is empty, we have to do a regular findLeaf and store the
            # #         // results in the path vector
            current = self
            # #         current = this;
            bitsBefore = 0
            # #         bitsBefore = 0;
            onesBefore = 0
            # #         onesBefore = 0;
        else:
            pass
            # #     } else {
            # #         // Start at the end of the path (which is a leaf), and move up until we
            # #         // reach the subtree that also contains the desired bit
            k = len(path) - 1
            # #         unsigned long k = path.size() - 1;
            nesbo = path[k]
            # #         Nesbo nesbo = path[k];
            # #         // While the subtree nested at `nesbo.node` is entirely before, or
            # #         // entirely after the bit we are looking for...
            while n < nesbo.bitsBefore or nesbo.bitsBefore + nesbo.size <= n:
                pass
                # #         while (n < nesbo.bitsBefore || nesbo.bitsBefore + nesbo.size <= n) {
                # #             // Go `up` one level if possible. Else, we are at the root so do a
                # #             // regular search from there
                path.pop()
                # #             path.pop_back();
                if k == 0:
                    pass
                    # #             if (k == 0) {
                    current = nesbo.node
                    # #                 current = nesbo.node;
                    break
                    # #                 break;
                    # #             }
                k -= 1
                # #             k -= 1;
                nesbo = path[k]
                # #             nesbo = path[k];
            # #         };
            # #         // If we didn't exit early, then set the start of the search path to the
            # #         // last entry of path that still exists
            if current is None:
                pass
                # #         if (current == nullptr) {
                current = nesbo.node.node.internalNode.entries[nesbo.index].p
                # #             current = nesbo.node->node.internalNode->entries[nesbo.index].P;
                bitsBefore = nesbo.bitsBefore
                # #             bitsBefore = nesbo.bitsBefore;
                onesBefore = nesbo.onesBefore
                # #             onesBefore = nesbo.onesBefore;
                # #         }
            # #     }
        # #     // Now, we do a regular search from `current` and add the intermediate
        # #     // nodes we visit to `path`
        while not current.isLeaf:
            pass
            # #     while (!current->isLeaf) {
            record = current.find_child(n - bitsBefore)
            # #         auto record = current->findChild(n - bitsBefore);
            bitsBefore += record.b
            # #         bitsBefore += record.b;
            onesBefore += record.o
            # #         onesBefore += record.o;
            next_ = current.node.internalNode.entries[record.i]
            # #         auto next = current->node.internalNode->entries[record.i];
            xxx = Nesbo(current, record.i, next_.b, bitsBefore, onesBefore)
            path.append(xxx)
            # #         path.emplace_back(current, record.i, next.b, bitsBefore, onesBefore);
            current = next_.p
            # #         current = next.P;
            # #     }
        xxx: InternalNode.Entry = InternalNode.Entry().init_with_bop(bitsBefore, onesBefore, current)
        return xxx
        # #     return {bitsBefore, onesBefore, current};
        # # }
        pass

    def rank1(self, n: int, path):
        pass
        #     /**
        #      * Performs the `rank` operation on the bitvector represented by this tree
        #      * @param n  an integer with 0 <= n < (numbers of bits in the tree)
        #      * @return the number of 1-bits in the tree up to position n
        #      */
        #     unsigned long rank1(unsigned long, vector<Nesbo> *path = nullptr);
        # # unsigned long TTree::rank1(unsigned long n, vector<Nesbo> *path) {
        entry = self.find_leaf(n, path)
        # #     auto entry = findLeaf(n, path);
        bv = entry.p.node.leafNode.bv
        # #     auto &bv = entry.P->node.leafNode->bv;
        return entry.o + bv.rank1(n - entry.b)
        # #     return entry.o + bv.rank1(n - entry.b);
        # # }
        pass

    def access(self, n: int, path):
        pass
        #     /**
        #      * Performs the `access` operation on this subtree
        #      * @param n the index of a bit in the `TTree`
        #      * @return the value of the `n`-th bit in the tree
        #      */
        #     bool access(unsigned long, vector<Nesbo> *path = nullptr);
        # # bool TTree::access(unsigned long n, vector<Nesbo> *path) {
        entry = self.find_leaf(n, path)
        # #     auto entry = findLeaf(n, path);
        return entry.p.node.leafNode.bv[n - entry.b]
        # #     return entry.P->node.leafNode->bv[n - entry.b];
        # # }
        pass

    def set_bit(self, n: int, b: bool, path: list[Nesbo]):
        pass
        #     /**
        #      * Sets the bit at position n to the value of b
        #      * @param n the index of the bit to set
        #      * @param b the value to set the bit to
        #      * @return true if this bit changed, e.g.
        #      *      if the previous value of the bit was unequal to b
        #      */
        #     bool setBit(unsigned long, bool, vector<Nesbo> *path = nullptr);
        # # bool TTree::setBit(unsigned long n, bool b, vector<Nesbo> *path) {
        # #     // Find the leaf node that contains this bit
        entry: InternalNode.Entry = self.find_leaf(n, path)
        # #     auto entry = findLeaf(n, path);
        bv: BitVector = entry.p.node.leafNode.bv
        # #     BitVector<> &bv = entry.P->node.leafNode->bv;
        changed: bool = bv.set(n - entry.b, b)
        # #     bool changed = bv.set(n - entry.b, b);
        if changed:
            # #     if (changed) {
            # #         // Change the one-counters all the way up from this leaf
            d = 1  # count up if 0 -> 1
            if not b:
                d = -1  # count down if 1 -> 0
            # #         long d = b ? 1 : -1;
            entry.p.updateCounters(0, d)
            # #         entry.P->updateCounters(0, d);
            if path is not None and len(path) != 0:
                # #         if (path != nullptr && !path->empty()) {
                for i in path:
                    # #             for (auto &i : *path) {
                    # #                 // For each entry where `bitsBefore` includes the bit we changed,
                    # #                 // we must also increase/decrease the value of `onesBefore` by 1
                    if n < i.bitsBefore:
                        # #                 if (n < i.bitsBefore) {
                        i.onesBefore += d
                        # #                     i.onesBefore += d;
                        # #                 }
                    # #             }
                    # #         }
                # #     }
        return changed
        # #     return changed;
        # # }
        pass

    def bits(self):
        pass
        #     /**
        #      * Gets the total number of bits covered by this subtree
        #      * @return the total number of bits in this node's subtree
        #      */
        #     unsigned long bits();
        # # unsigned long TTree::bits() {
        if self.isLeaf:
            # #     if (isLeaf) {
            return self.node.leafNode.bits()
            # #         return node.leafNode->bits();
        else:
            # #     } else {
            return self.node.internalNode.bits()
            # #         return node.internalNode->bits();
            # #     }
        # # }
        pass

    def ones(self):
        pass
        #     /**
        #      * Gets the total number of ones covered by this subtree
        #      * @return the total number of 1-bits in this node's subtree
        #      */
        #     unsigned long ones();
        # # unsigned long TTree::ones() {
        if self.isLeaf:
            # #     if (isLeaf) {
            return self.node.leafNode.ones()
            # #         return node.leafNode->ones();
        else:
            # #     } else {
            return self.node.internalNode.ones()
            # #         return node.internalNode->ones();
            # #     }
        # # }
        pass

    def updateCounters(self, dBits: int, dOnes: int):
        pass
        #     /**
        #      * Changes the values of the counters `b` and `o` of all the nodes whose
        #      * subtree contains this node. Used to update these counters after modifying
        #      * the underlying bitvector, or the structure of the tree.
        #      *
        #      * @param dBits the change in the number of bits (e.g. 4 when 4 bits were inserted)
        #      * @param dOnes the change in the number of ones (e.g. 2 when 2 zeros were set to ones)
        #      */
        #     void updateCounters(long, long);
        # # void TTree::updateCounters(long dBits, long dOnes) {
        current: TTree = self
        # #     TTree *current = this;
        # #     // Go up in the tree until we reach the root
        while current.parent is not None:
            # #     while (current->parent != nullptr) {
            # #         // Take the entry in `current`s parent that points to `current`,
            # #         // and update its `b` and `o` counters.
            try:
                parent = current.parent
                # #         auto parent = current->parent;
                entry = parent.node.internalNode.entries[current.indexInParent]
                # #         auto &entry = parent->node.internalNode->entries[current->indexInParent];
                entry.b += dBits
                # #         entry.b += dBits;
                entry.o += dOnes
                # #         entry.o += dOnes;
                current = parent
                # #         current = parent;
                # #     }
            except Exception as e:
                pass
            # # }
        pass

    def insertBlock(self, index, path: list[Nesbo] = None) -> TTree:
        pass
        #     /**
        #      * Inserts one block of k^2 bits at the specified position
        #      *
        #      * @return the new root if the tree's root changed, nullptr otherwise
        #      */
        #     TTree *insertBlock(long unsigned, vector<Nesbo> *path = nullptr);
        # # TTree *TTree::insertBlock(long unsigned index, vector<Nesbo> *path) {
        return self.insertBits(index, Parameters.BLOCK_SIZE, path)
        # #     return this->insertBits(index, BLOCK_SIZE, path);
        # # }
        pass

    def deleteBlock(self, index: int, path: list[Nesbo] = None) -> TTree:
        pass
        #     /**
        #      * Deletes a block of k^2 bits at the specified position
        #      *
        #      * @return the new root if the tree's root changed, nullptr otherwise
        #      */
        #     TTree *deleteBlock(long unsigned, vector<Nesbo> *path = nullptr);
        # # TTree *TTree::deleteBlock(long unsigned index, vector<Nesbo> *path) {
        return self.deleteBits(index, Parameters.BLOCK_SIZE, path)
        # #     return this->deleteBits(index, BLOCK_SIZE, path);
        # # }
        pass

    def memoryUsage(self):
        pass
        #     unsigned long memoryUsage();
        # # unsigned long TTree::memoryUsage() {
        result = sys.getsizeof(TTree)
        # #     unsigned long result = sizeof(TTree);
        if self.isLeaf:
            # #     if (isLeaf) {
            result += self.node.leafNode.bv.memory_usage()
            # #         result += node.leafNode->bv.memoryUsage();
        else:
            # #     } else {
            result += sys.getsizeof(InternalNode)
            # #         result += sizeof(InternalNode);
            entries = self.node.internalNode.entries
            # #         auto &entries = node.internalNode->entries;
            for entry in entries:
                # #         for (auto &entry : entries) {
                if entry.p is not None:
                    # #             if (entry.P != nullptr) {
                    result += entry.p.memoryUsage()
                    # #                 result += entry.P->memoryUsage();
                    # #             }
                # #         }
            # #     }
        return result
        # #     return result;
        # # }
        pass

    def insertBits(self, index: int, count: int, path: list[Nesbo]) -> TTree:
        pass
        # private:
        #     /**
        #      * Inserts the given number of bits (set to zero) at the given position in the tree
        #      *
        #      * @param index the position at which to insert bits
        #      * @param count the number of bits to insert
        #      *
        #      * @return nullptr in most cases, but a pointer to the new root node
        #      * if this insert operation created a new root (e.g. when the height of the
        #      * tree increases)
        #      *
        #      * This method assumes that after the insertion, the size of the relevant
        #      * leaf vector is at most one block over the maximum, which can only be
        #      * guaranteed when only a single block is inserted. That is why this method
        #      * is private, and insertBlock is public.
        #          */
        #     TTree *
        #     insertBits(long unsigned, long unsigned, vector<Nesbo> *path = nullptr);
        # # TTree *TTree::insertBits(long unsigned index, long unsigned count,
        # #                          vector<Nesbo> *path) {
        entry = self.find_leaf(index, path)
        # #     auto entry = findLeaf(index, path);
        leaf = entry.p
        # #     auto leaf = entry.P;
        bv = leaf.node.leafNode.bv
        # #     auto &bv = leaf->node.leafNode->bv;
        bv.insert(index - entry.b, count)
        # #     bv.insert(index - entry.b, count);
        leaf.updateCounters(count, 0)
        # #     leaf->updateCounters(count, 0);
        # #     // Split this node up into two if it exceeds the size limit
        xxx = leaf.checkSizeUpper()
        return xxx
        # #     return leaf->checkSizeUpper();
        # # }
        pass

    def deleteBits(self, index: int, count: int, path: list[Nesbo] = None):
        pass
        #     /**
        #      * Deletes the given number of bits in the subtree, assuming they are all in the
        #      * same leaf node. This will be the case if a group of k^2 bits are deleted
        #      *
        #      * @param index the position of the first bit to delete
        #      * @param count the number of bits to delete
        #      *
        #      * @return nullptr in most cases, but a pointer to the new root node
        #      * if this insert operation changed the root (e.g. when the height of the
        #      * tree decreases)
        #      *
        #      * This method assumes that after the deletion, the size of the relevant
        #      * leaf vector is at most one block under the minimum, which can only be
        #      * guaranteed when only a single block is deleted. That is why this method
        #      * is private, and deleteBlock is public.
        #      */
        #     TTree *
        #     deleteBits(long unsigned, long unsigned, vector<Nesbo> *path = nullptr);
        # # TTree *TTree::deleteBits(long unsigned index, long unsigned count,
        # #                          vector<Nesbo> *path) {
        entry = self.find_leaf(index, path)
        # #     auto entry = findLeaf(index, path);
        leaf = entry.p
        # #     auto leaf = entry.P;
        bv = leaf.node.leafNode.bv
        # #     auto &bv = leaf->node.leafNode->bv;
        start = index - entry.b
        # #     long unsigned start = index - entry.b;
        end = start + count
        # #     long unsigned end = start + count;
        deletedOnes = bv.range_rank1(start, end)
        # #     long unsigned deletedOnes = bv.rangeRank1(start, end);
        bv.erase(start, end)
        # #     bv.erase(start, end);
        leaf.updateCounters(-count, -deletedOnes)
        # #     leaf->updateCounters(-count, -deletedOnes);
        return leaf.checkSizeLower()
        # #     return leaf->checkSizeLower();
        # # }
        pass

    def checkSizeUpper(self) -> TTree:
        pass
        #     /**
        #      * Checks if this node satisfies the maximum size for an internal node
        #      * or leaf node. If not, tries to spill a node to a sibling, or if that
        #      * fails will split this node into two and recursively check the parent
        #      *
        #      * @return nullptr in most cases, but returns the new root if it has
        #      *         changed, e.g. if the height of the tree has increased
        #      */
        #     TTree *checkSizeUpper();
        # # TTree *TTree::checkSizeUpper() {
        if self.isLeaf and self.size() > Parameters.leafSizeMax:
            # #     if (isLeaf && size() > leafSizeMax) {
            if not self.trySpillLeaf():
                # #         if (!trySpillLeaf()) {
                return self.splitLeaf()
                # #             return splitLeaf();
                # #         }
        elif not self.isLeaf and self.size() > Parameters.nodeSizeMax:
            # #     } else if (!isLeaf && size() > nodeSizeMax) {
            if not self.trySpillInternal():
                # #         if (!trySpillInternal()) {
                return self.splitInternal()
                # #             return splitInternal();
                # #         }
            # #     }
        return None
        # #     return nullptr;
        # # }
        pass

    def checkSizeLower(self) -> TTree | None:
        pass
        #     /**
        #      * Checks if this node satisfies the minimum size for an internal node
        #      * or leaf node. If not, tries to steal a node from a sibling, or if that
        #      * fails will merge this node with one of the siblings, and recursively
        #      * check the parent
        #      *
        #      * @return nullptr in most cases, but returns the new root if it has
        #      *         changed, e.g. if the height of the tree has decreased
        #      */
        #     TTree *checkSizeLower();
        # # TTree *TTree::checkSizeLower() {
        isRoot: bool = (self.parent is None)
        # #     bool isRoot = (parent == nullptr);
        if isRoot and (self.size() >= 2 or self.isLeaf):
            # #     if (isRoot && (size() >= 2 || isLeaf)) {
            return None
            # #         return nullptr;
        elif self.isLeaf and self.size() < Parameters.leafSizeMin:
            # #     } else if (isLeaf && size() < leafSizeMin) {
            if not self.tryStealLeaf():
                # #         if (!tryStealLeaf()) {
                return self.mergeLeaf()
                # #             return mergeLeaf();
                # #         }
        elif not self.isLeaf and self.size() < Parameters.nodeSizeMin:
            # #     } else if (!isLeaf && size() < nodeSizeMin) {
            if not self.tryStealInternal():
                # #         if (!tryStealInternal()) {
                return self.mergeInternal()
                # #             return mergeInternal();
                # #         }
                # #     }
        return None
        # #     return nullptr;
        # # }
        pass

    def trySpillInternal(self) -> bool:
        pass
        #     /**
        #      * Tries to move a child of an internal node to a sibling, and returns
        #      * whether it succeeded
        #      *
        #      * @return true if a spill could be done, false if not (e.g. if both of
        #      *         the nodes siblings are of maximum size)
        #      */
        #     bool trySpillInternal();
        # # bool TTree::trySpillInternal() {
        if self.parent is None:
            # #     if (parent == nullptr) {
            return False
            # #         return false;
            # #     }
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        n: int = self.parent.size()
        # #     unsigned long n = parent->size();
        entries = self.parent.node.internalNode.entries
        # #     auto &entries = parent->node.internalNode->entries;
        if idx > 0 and entries[idx - 1].p.size() < Parameters.nodeSizeMax:
            # #     if (idx > 0 && entries[idx - 1].P->size() < nodeSizeMax) {
            self.moveLeftInternal()
            # #         this->moveLeftInternal();
            return True
            # #         return true;
        elif idx + 1 < n and entries[idx + 1].p.size() < Parameters.nodeSizeMax:
            # #     } else if (idx + 1 < n && entries[idx + 1].P->size() < nodeSizeMax) {
            self.moveRightInternal()
            # #         this->moveRightInternal();
            return True
            # #         return true;
        else:
            # #     } else {
            return False
            # #         return false;
            # #     }
        # # }
        pass

    def trySpillLeaf(self) -> bool:
        pass
        #     /**
        #      * Tries to move a child of a leaf node to a sibling, and returns
        #      * whether it succeeded
        #      *
        #      * @return true if a spill could be done, false if not (e.g. if both of
        #      *         the nodes siblings are of maximum size)
        #      */
        #     bool trySpillLeaf();
        # # bool TTree::trySpillLeaf() {
        if self.parent is None:
            # #     if (parent == nullptr) {
            return False
            # #         return false;
            # #     }
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        n: int = self.parent.size()
        # #     unsigned long n = parent->size();
        entries = self.parent.node.internalNode.entries
        # #     auto &entries = parent->node.internalNode->entries;
        if idx > 0 and entries[idx - 1].p.size() < Parameters.leafSizeMax:
            # #     if (idx > 0 && entries[idx - 1].P->size() < leafSizeMax) {
            self.moveLeftLeaf()
            # #         this->moveLeftLeaf();
            return True
            # #         return true;
        elif idx + 1 < n and entries[idx + 1].p.size() < Parameters.leafSizeMax:
            # #     } else if (idx + 1 < n && entries[idx + 1].P->size() < leafSizeMax) {
            self.moveRightLeaf()
            # #         this->moveRightLeaf();
            return True
            # #         return true;
        else:
            # #     } else {
            return False
            # #         return false;
            # #     }
        # # }
        pass

    def splitInternal(self) -> TTree:
        pass
        #     /**
        #      * Splits this node into two nodes of minimum size, and recursively
        #      * checks the rest of the tree for meeting size requirements
        #      *
        #      * @return nullptr in most cases, but returns the new root if this operation
        #      *         causes the tree's height to increase, which changes the root
        #      */
        #     TTree *splitInternal();
        # # TTree *TTree::splitInternal() {
        entries = self.node.internalNode.entries
        # #     auto &entries = this->node.internalNode->entries;
        n: int = self.size()
        # #     unsigned long n = this->size();
        mid: int = int(n / 2)
        # #     unsigned long mid = n / 2;
        newNode = TTree()
        # #     auto newNode = new TTree();
        newNode.parent = self.parent
        # #     newNode->parent = parent;
        newNode.isLeaf = False
        # #     newNode->isLeaf = false;
        newNode.node.leafNode = None
        # #     newNode->node.leafNode = nullptr;
        newNode.node.internalNode = InternalNode()
        # #     newNode->node.internalNode = new InternalNode();
        d_b: int = 0
        d_o: int = 0
        # #     unsigned long d_b = 0, d_o = 0; // Count bits/ones in right half
        for i in range(mid, n):
            # #     for (unsigned long i = mid; i < n; i++) {
            entry = entries[i]
            # #         auto entry = entries[i];
            if entry.p is not None:
                # #         if (entry.P != nullptr) {
                entry.p.parent = newNode
                # #             entry.P->parent = newNode;
                entry.p.indexInParent = i - mid
                if i - mid >= 3:
                    pass  # debug  # 2024/8/23
                # #             entry.P->indexInParent = i - mid;
                d_b += entry.b
                # #             d_b += entry.b;
                d_o += entry.o
                # #             d_o += entry.o;
                # #         }
            newNode.node.internalNode.entries[i - mid] = entry
            # #         newNode->node.internalNode->entries[i - mid] = entry;
            entries[i] = InternalNode.Entry()
            # #         entries[i] = InternalNode::Entry();
            # #     }
        newNode.node.internalNode.size = n - mid
        # #     newNode->node.internalNode->size = n - mid;
        self.node.internalNode.size = mid
        # #     node.internalNode->size = mid;
        if self.parent is None:
            # #     if (parent == nullptr) {
            newRoot = TTree().init_with_arguments(self, newNode)
            # #         auto *newRoot = new TTree(this, newNode);
            return newRoot
            # #         return newRoot;
        else:
            # #     } else {
            self.parent.node.internalNode.entries[self.indexInParent].b -= d_b
            # #         parent->node.internalNode->entries[indexInParent].b -= d_b;
            self.parent.node.internalNode.entries[self.indexInParent].o -= d_o
            # #         parent->node.internalNode->entries[indexInParent].o -= d_o;
            xxx = InternalNode.Entry().init_with_bop(d_b, d_o, newNode)
            self.parent.node.internalNode.insert(self.indexInParent + 1, xxx)
            # #         parent->node.internalNode->insert(indexInParent + 1,
            # #                                           {d_b, d_o, newNode});
            return self.parent.checkSizeUpper()
            # #         return parent->checkSizeUpper();
            # #     }
        # # }
        pass

    def splitLeaf(self) -> TTree:
        pass
        #     /**
        #      * Splits this node into two nodes of minimum size, and recursively
        #      * checks the rest of the tree for meeting size requirements
        #      *
        #      * @return nullptr in most cases, but returns the new root if this operation
        #      *         causes the tree's height to increase, which changes the root
        #      */
        #     TTree *splitLeaf();
        # # TTree *TTree::splitLeaf() {
        n: int = self.node.leafNode.bits()
        # #     unsigned long n = this->node.leafNode->bits();
        mid: int = int(n / 2)
        # #     unsigned long mid = n / 2;
        mid -= mid % Parameters.BLOCK_SIZE
        # #     mid -= mid % BLOCK_SIZE;
        left = self.node.leafNode.bv
        # #     auto &left = this->node.leafNode->bv;
        right = BitVector().init_with_range(left, mid, n)
        # #     auto right = BitVector<>(left, mid, n);
        left.erase(mid, n)
        # #     left.erase(mid, n);
        newNode: TTree = TTree().init_with_bv(right)
        # #     auto *newNode = new TTree(right);
        if self.parent is None:
            # #     if (parent == nullptr) {
            newRoot: TTree = TTree().init_with_arguments(self, newNode)
            # #         auto *newRoot = new TTree(this, newNode);
            return newRoot
            # #         return newRoot;
        else:
            # #     } else {
            idx: int = self.indexInParent
            # #         unsigned long idx = indexInParent;
            newNode.parent = self.parent
            # #         newNode->parent = parent;
            entry = InternalNode.Entry().init_with_p(newNode)
            # #         InternalNode::Entry entry(newNode);
            self.parent.node.internalNode.insert(self.indexInParent + 1, entry)
            # #         parent->node.internalNode->insert(indexInParent + 1, entry);
            self.parent.node.internalNode.entries[idx].b -= entry.b
            # #         parent->node.internalNode->entries[idx].b -= entry.b;
            self.parent.node.internalNode.entries[idx].o -= entry.o
            # #         parent->node.internalNode->entries[idx].o -= entry.o;
            return self.parent.checkSizeUpper()
            # #         return parent->checkSizeUpper();
            # #     }
        # # }
        pass

    def tryStealInternal(self):
        pass
        #     /**
        #      * Tries to steal a node from one of this node's siblings, and returns true
        #      * if it succeeds
        #      *
        #      * @return true if this operation could be done, false otherwise (e.g. if both
        #      *         of this node's siblings are of minimum size)
        #      */
        #     bool tryStealInternal();
        # # bool TTree::tryStealInternal() {
        if self.parent is None:
            # #     if (parent == nullptr) {
            return False
            # #         return false;
            # #     }
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        n: int = self.parent.size()
        # #     unsigned long n = parent->size();
        entries: list[InternalNode.Entry] = self.parent.node.internalNode.entries
        # #     auto &entries = parent->node.internalNode->entries;
        if idx > 0 and entries[idx - 1].p.size() > Parameters.nodeSizeMin:
            # #     if (idx > 0 && entries[idx - 1].P->size() > nodeSizeMin) {
            entries[idx - 1].p.moveRightInternal()
            # #         entries[idx - 1].P->moveRightInternal();
            return True
            # #         return true;
        elif idx + 1 < n and entries[idx + 1].p.size() > Parameters.nodeSizeMin:
            # #     } else if (idx + 1 < n && entries[idx + 1].P->size() > nodeSizeMin) {
            entries[idx + 1].p.moveLeftInternal()
            # #         entries[idx + 1].P->moveLeftInternal();
            return True
            # #         return true;
        else:
            # #     } else {
            return False
            # #         return false;
            # #     }
        # # }
        pass

    def tryStealLeaf(self):
        pass
        #     /**
        #      * Tries to steal a node from one of this node's siblings, and returns true
        #      * if it succeeds
        #      *
        #      * @return true if this operation could be done, false otherwise (e.g. if both
        #      *         of this node's siblings are of minimum size)
        #      */
        #     bool tryStealLeaf();
        # # bool TTree::tryStealLeaf() {
        if self.parent is None:
            # #     if (parent == nullptr) {
            return False
            # #         return false;
            # #     }
        idx = self.indexInParent
        # #     unsigned long idx = indexInParent;
        n = self.parent.size()
        # #     unsigned long n = parent->size();
        entries = self.parent.node.internalNode.entries
        # #     auto &entries = parent->node.internalNode->entries;
        if idx > 0 and entries[idx - 1].p.size() > Parameters.leafSizeMin:
            # #     if (idx > 0 && entries[idx - 1].P->size() > leafSizeMin) {
            entries[idx - 1].p.moveRightLeaf()
            # #         entries[idx - 1].P->moveRightLeaf();
            return True
            # #         return true;
        elif idx + 1 < n and entries[idx + 1].p.size() > Parameters.leafSizeMin:
            # #     } else if (idx + 1 < n && entries[idx + 1].P->size() > leafSizeMin) {
            entries[idx + 1].p.moveLeftLeaf()
            # #         entries[idx + 1].P->moveLeftLeaf();
            return True
            # #         return true;
        else:
            # #     } else {
            return False
            # #         return false;
            # #     }
        # # }
        pass

    def mergeInternal(self):
        pass
        #     /**
        #      * Merges this node with a sibling, and recursively checks the rest of the
        #      * tree for meeting size constraints
        #      *
        #      * @return nullptr usually, but returns the new root if it changed due to this
        #      *         operation, e.g. when the height of the tree changed
        #      */
        #     TTree *mergeInternal();
        # # TTree *TTree::mergeInternal() {
        # #     // If we are the root, and we are too small, then we have only one child
        if self.parent is None:
            # #     if (parent == nullptr) {
            # #         // Delete this, our only child should become the root
            child = self.node.internalNode.entries[0].p
            # #         TTree *child = node.internalNode->entries[0].P;
            # #         // Overwrite the pointer in the entry, so that it is not deleted
            self.node.internalNode.entries[0].p = None
            # #         node.internalNode->entries[0].P = nullptr;
            del self
            # #         delete this;
            child.parent = None
            # #         child->parent = nullptr;
            child.indexInParent = 0
            # #         child->indexInParent = 0;
            return child
            # #         return child;
            # #     }
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        # #     TTree *left = nullptr, *right = nullptr;
        if idx > 0:
            # #     if (idx > 0) {
            left: TTree = self.parent.node.internalNode.entries[idx - 1].p
            # #         left = parent->node.internalNode->entries[idx - 1].P;
            right: TTree = self
            # #         right = this;
            idx -= 1
            # #         idx--;
        else:
            # #     } else {
            left: TTree = self
            # #         left = this;
            right: TTree = self.parent.node.internalNode.entries[idx + 1].p
            # #         right = parent->node.internalNode->entries[idx + 1].P;
            # #     }
        # #     // Merge `left` and `right` into one node
        internalNode: InternalNode = left.node.internalNode
        # #     auto &internalNode = left->node.internalNode;
        n: int = right.size()
        # #     unsigned long n = right->size();
        d_b: int = 0
        d_o: int = 0
        # #     unsigned long d_b = 0, d_o = 0;
        for i in range(n):
            # #     for (unsigned i = 0; i < n; i++) {
            entry: InternalNode.Entry = right.node.internalNode.entries[i]
            # #         auto entry = right->node.internalNode->entries[i];
            right.node.internalNode.entries[i].p = None
            # #         right->node.internalNode->entries[i].P = nullptr;
            d_b += entry.b
            # #         d_b += entry.b;
            d_o += entry.o
            # #         d_o += entry.o;
            entry.p.parent = left
            # #         entry.P->parent = left;
            self.node.internalNode.append(entry)  # node?
            # #         internalNode->append(entry);
            # #     }
        # #     // Delete the right child, and update the b and o counters for left
        self.parent.node.internalNode.remove(idx + 1)
        # #     parent->node.internalNode->remove(idx + 1);
        self.parent.node.internalNode.entries[idx].b += d_b
        # #     parent->node.internalNode->entries[idx].b += d_b;
        self.parent.node.internalNode.entries[idx].o += d_o
        # #     parent->node.internalNode->entries[idx].o += d_o;
        return self.parent.checkSizeLower()
        # #     return parent->checkSizeLower();
        # # }
        pass

    def mergeLeaf(self):
        pass
        #     /**
        #      * Merges this node with a sibling, and recursively checks the rest of the
        #      * tree for meeting size constraints
        #      *
        #      * @return nullptr usually, but returns the new root if it changed due to this
        #      *         operation, e.g. when the height of the tree changed
        #      */
        #     TTree *mergeLeaf();
        # # TTree *TTree::mergeLeaf() {
        if self.parent is None:
            # #     if (parent == nullptr) {
            return None
            # #         return nullptr;
            # #     }
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        # #     TTree *left = nullptr, *right = nullptr;
        if idx > 0:
            # #     if (idx > 0) {
            left: TTree = self.parent.node.internalNode.entries[idx - 1].p
            # #         left = parent->node.internalNode->entries[idx - 1].P;
            right: TTree = self
            # #         right = this;
            idx -= 1
            # #         idx--;
        else:
            # #     } else {
            left: TTree = self
            # #         left = this;
            right: TTree = self.parent.node.internalNode.entries[idx + 1].p
            # #         right = parent->node.internalNode->entries[idx + 1].P;
            # #     }
        leftBits: BitVector = left.node.leafNode.bv
        # #     auto &leftBits = left->node.leafNode->bv;
        rightBits: BitVector = right.node.leafNode.bv
        # #     auto &rightBits = right->node.leafNode->bv;
        # #     // Append `right`s bits to `left`
        leftBits.append(rightBits, 0, rightBits.size())
        # #     leftBits.append(rightBits, 0, rightBits.size());
        # #     // Update the b and o for `left`, and delete `right`
        d_b: int = self.parent.node.internalNode.entries[idx + 1].b
        # #     unsigned long d_b = parent->node.internalNode->entries[idx + 1].b;
        d_o: int = self.parent.node.internalNode.entries[idx + 1].o
        # #     unsigned long d_o = parent->node.internalNode->entries[idx + 1].o;
        self.parent.node.internalNode.entries[idx].b += d_b
        # #     parent->node.internalNode->entries[idx].b += d_b;
        self.parent.node.internalNode.entries[idx].o += d_o
        # #     parent->node.internalNode->entries[idx].o += d_o;
        self.parent.node.internalNode.remove(idx + 1)
        # #     parent->node.internalNode->remove(idx + 1);
        return self.parent.checkSizeLower()
        # #     return parent->checkSizeLower();
        # # }
        pass

    def moveLeftInternal(self):
        pass
        #     /**
        #      * Moves a single child (the leftmost child) of this node to the end of the
        #      * node's left sibling
        #      */
        #     void moveLeftInternal();
        # # void TTree::moveLeftInternal() {
        parent: TTree = self.parent
        # #     TTree *parent = this->parent;
        idx: int = self.indexInParent
        # #     unsigned long idx = this->indexInParent;
        sibling: TTree = self.parent.node.internalNode.entries[idx - 1].p
        # #     TTree *sibling = parent->node.internalNode->entries[idx - 1].P;
        # #     // Move the first child of `this` to the end of the left sibling
        toMove: InternalNode.Entry = self.node.internalNode.pop_first()
        # #     InternalNode::Entry toMove = this->node.internalNode->popFirst();
        d_b: int = toMove.b
        # #     unsigned long d_b = toMove.b;
        d_o: int = toMove.o
        # #     unsigned long d_o = toMove.o;
        toMove.p.parent = sibling
        # #     toMove.P->parent = sibling;
        sibling.node.internalNode.append(toMove)
        # #     sibling->node.internalNode->append(toMove);
        # #     // Finally, update the parent's b and o counters for `this` and `sibling`
        # #     // The number of bits/ones in `toMove` is subtracted from `this`, but added to `sibling`
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx].o -= d_o
        # #     parent->node.internalNode->entries[idx].o -= d_o;
        self.parent.node.internalNode.entries[idx - 1].b += d_b
        # #     parent->node.internalNode->entries[idx - 1].b += d_b;
        self.parent.node.internalNode.entries[idx - 1].o += d_o
        # #     parent->node.internalNode->entries[idx - 1].o += d_o;
        # # }
        pass

    def moveRightInternal(self):
        pass
        #     /**
        #      * Moves the rightmost child of this node to the start of the right sibling
        #      */
        #     void moveRightInternal();
        # # void TTree::moveRightInternal() {
        idx: int = self.indexInParent
        # #     unsigned long idx = this->indexInParent;
        sibling: TTree = self.parent.node.internalNode.entries[idx + 1].p
        # #     TTree *sibling = parent->node.internalNode->entries[idx + 1].P;
        # #     // Move the last child of `this` to the start of the left sibling
        toMove: InternalNode.Entry = self.node.internalNode.pop_last()
        # #     InternalNode::Entry toMove = this->node.internalNode->popLast();
        d_b: int = toMove.b
        # #     unsigned long d_b = toMove.b;
        d_o: int = toMove.o
        # #     unsigned long d_o = toMove.o;
        toMove.p.parent = sibling
        # #     toMove.P->parent = sibling;
        sibling.node.internalNode.insert(0, toMove)
        # #     sibling->node.internalNode->insert(0, toMove);
        # #     // Finally, update the parent's b and o counters for `this` and `sibling`
        # #     // The number of bits/ones in `toMove` is subtracted from `this`, but added to `sibling`
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx].o -= d_o
        # #     parent->node.internalNode->entries[idx].o -= d_o;
        self.parent.node.internalNode.entries[idx + 1].b += d_b
        # #     parent->node.internalNode->entries[idx + 1].b += d_b;
        self.parent.node.internalNode.entries[idx + 1].o += d_o
        # #     parent->node.internalNode->entries[idx + 1].o += d_o;
        # # }
        pass

    def moveLeftLeaf(self):
        pass
        #     /**
        #      * Moves the leftmost block of k^2 bits to the end of the left sibling
        #      */
        #     void moveLeftLeaf();
        # # void TTree::moveLeftLeaf() {
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        sibling: TTree = self.parent.node.internalNode.entries[idx - 1].p
        # #     TTree *sibling = parent->node.internalNode->entries[idx - 1].P;
        # #     // Take the first k*k block of `this`, and append it to `sibling`
        right: BitVector = self.node.leafNode.bv
        # #     BitVector<> &right = node.leafNode->bv;
        left: BitVector = sibling.node.leafNode.bv
        # #     BitVector<> &left = sibling->node.leafNode->bv;
        d_b: int = Parameters.BLOCK_SIZE
        # #     unsigned long d_b = BLOCK_SIZE;
        d_o: int = right.rank1(Parameters.BLOCK_SIZE)
        # #     unsigned long d_o = right.rank1(BLOCK_SIZE);
        left.append(right, 0, Parameters.BLOCK_SIZE)
        # #     left.append(right, 0, BLOCK_SIZE);
        right.erase(0, Parameters.BLOCK_SIZE)
        # #     right.erase(0, BLOCK_SIZE);
        # #     // Update the parent's b and o counters
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx].o -= d_o
        # #     parent->node.internalNode->entries[idx].o -= d_o;
        self.parent.node.internalNode.entries[idx - 1].b += d_b
        # #     parent->node.internalNode->entries[idx - 1].b += d_b;
        self.parent.node.internalNode.entries[idx - 1].o += d_o
        # #     parent->node.internalNode->entries[idx - 1].o += d_o;
        # # }
        pass

    def moveRightLeaf(self):
        pass
        #     /**
        #      * Moves the rightmost block of k^2 bits to the start of the right sibling
        #      */
        #     void moveRightLeaf();
        # # void TTree::moveRightLeaf() {
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        sibling: TTree = self.parent.node.internalNode.entries[idx + 1].p
        # #     TTree *sibling = parent->node.internalNode->entries[idx + 1].P;
        # #     // Take the first k*k block of `this`, and append it to `sibling`
        left: BitVector = self.node.leafNode.bv
        # #     BitVector<> &left = node.leafNode->bv;
        right: BitVector = sibling.node.leafNode.bv
        # #     BitVector<> &right = sibling->node.leafNode->bv;
        hi: int = left.size()
        # #     unsigned long hi = left.size();
        lo: int = hi - Parameters.BLOCK_SIZE
        # #     unsigned long lo = hi - BLOCK_SIZE;
        d_b: int = Parameters.BLOCK_SIZE
        # #     unsigned long d_b = BLOCK_SIZE;
        d_o: int = left.range_rank1(lo, hi)
        # #     unsigned long d_o = left.rangeRank1(lo, hi);
        right.insert_range(0, left, lo, hi)
        # #     right.insert(0, left, lo, hi);
        left.erase(lo, hi)
        # #     left.erase(lo, hi);

        # #     // Update the parent's b and o counters
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx].o -= d_o
        # #     parent->node.internalNode->entries[idx].o -= d_o;
        self.parent.node.internalNode.entries[idx + 1].b += d_b
        # #     parent->node.internalNode->entries[idx + 1].b += d_b;
        self.parent.node.internalNode.entries[idx + 1].o += d_o
        # #     parent->node.internalNode->entries[idx + 1].o += d_o;
        # # }
        pass

    # };
    pass  # end of class TTree

# #endif //DK2TREE_TTREE_H

# --------------------------------------------------------------------
# """
# class TTree:
# # // Constructors and destructors for data types that can't be in TTree.h
