"""
LTree.py << LTree.cpp
2024/8/1, T. Masuda
Amagasa Laboratory, University of Tsukuba
"""
from __future__ import annotations
import sys
from typing import Optional
import copy

# //
# // Created by anneke on 18/12/18.
# //
# 
# #ifndef DK2TREE_LTREE_H
# #define DK2TREE_LTREE_H
# 
# #include "BitVector.h"
# #include <utility>
# #include "parameters.cpp"
from src.Parameters import Parameters
from src.BitVector import BitVector


class LRecord:
    pass
    # /// LRecord type containing the number pf preceding bits, and the
    # /// index of a child node in the parent's `entries` list
    # struct LRecord {

    def __init__(self, b: int, i: int):
        pass
        self.b: int = b  # the number pf preceding bits
        #       unsigned long b;
        self.i: int = i  # index of a child node in the parent's `entries` list
        #       unsigned long i;
        pass
    # };
    pass  # end of class Record

# /// The three main structs forming the tree
# /// `LTree` represents one node in the tree, which contains a pointer to its
# /// parent (or nullptr for the root) and the index in the parent,
# /// as well as either an LInternalNode containing entries (b, P) or a LLeafNode
# /// containing a BitVector
# struct LInternalNode;
# struct LLeafNode;
# struct LTree;


class LNesbo:
    pass
    # struct LNesbo {
    # public:

    def __init__(self):
        pass
        self.node: Optional[LTree] = None
        #     LTree *node;
        self.index: int = 0
        #     unsigned long index;
        self.size: int = 0
        #     unsigned long size;
        self.bits_before: int = 0
        #     unsigned long bitsBefore;
        pass

    def init_with_arguments(self, node: LTree, index: int, size: int, bits_before: int):
        pass
        #     LNesbo(LTree *node, unsigned long index, unsigned long size,
        #           unsigned long bitsBefore) :
        self.node: Optional[LTree] = node
        #       node(node),
        self.index: int = index
        #       index(index),
        self.size: int = size
        #       size(size),
        self.bits_before: int = bits_before
        #       bitsBefore(bitsBefore) {}
        return self
        pass
    # };
    pass  # end of class LNesbo


class LInternalNode:
    pass
    # /**
    #  * A struct for the internal nodes of the tree
    #  * This contains a list of Entries of the form <b, P>
    #  */
    # struct LInternalNode {

    class Entry:
        pass
#       struct Entry {

        def __init__(self):
            pass
            #         Entry() :
            self.b: int = 0
            #         unsigned long b;
            self.p: Optional[LTree] = None
            #         LTree *P;
            #                 b(0), P(nullptr) {}
            pass

        def init_with_p(self, p: LTree):
            pass
            #         explicit Entry(LTree *);
            # LInternalNode::Entry::Entry(LTree *P) :
            self.b: int = p.bits()
            #         b(P->bits()),
            self.p: LTree = p
            #         P(P) {}
            return self
            pass

        def init_with_arguments(self, b: int, p: LTree):
            pass
            #       /// Construct Entry from the three fields
            #       /// This method is not unused, but is used as an initializer list
            # Entry(unsigned long b, LTree *P) :
            self.__init__()
            self.b: int = b
            self.p: Optional[LTree] = p
            #           b(b), P(P) {}
            return self
            pass

        def remove(self):
            pass
            #         /**
            #          * This function is called in the destructor of `LInternalNode`, to
            #          * delete the child nodes. It is not part of the destructor of `Entry`,
            #          * since that causes problems when the struct is used in methods such as
            #          * `findLeaf`
            #          */
            #           void remove();
            # # void LInternalNode::Entry::remove() {
            self.p = None
            # #     delete P;
            # # }
            pass
        #     };
        pass  # end of class Entry

#     /// An array of pointers to the child nodes
#     /// This is one more than the maximum, so that we can split nodes
#     /// after insertion instead of before
    def __init__(self):
        pass
#       /// The default constructor creates an empty internal node
#       LInternalNode() :
#             size(0),
#             entries{Entry()} {}
        self.size: int = 0
#       /// The number of children this node has
#       unsigned long size;
        self.entries = [self.Entry() for _ in range(Parameters.nodeSizeMax + 1)]
#       Entry entries[nodeSizeMax + 1];
        self.left: Optional[LTree] = None
        self.right: Optional[LTree] = None
        pass

    # def init_with_arguments(self, left: LTree, right: LTree, parent: Optional[LTree] = None):
    def init_with_arguments(self, left: LTree, right: LTree):
        pass
        #     /**
        #      * Creates a new internal node with the given two children
        #      *
        #      * @param left the first child of this node
        #      * @param right the second child of this node
        #      * @param parent the parent node, which has this as its internal node
        #      *        the left and right LTrees have their parent and indexInParent
        #      *        set correctly as well
        #      */
        #     LInternalNode(LTree *left, LTree *right, LTree *parent = nullptr);  // .h
        # # LInternalNode::LInternalNode(LTree *left, LTree *right, LTree *parent) :
        # #         size(2),
        # #         entries{Entry(left), Entry(right), Entry()} {
        # #     left->parent = parent;
        # #     left->indexInParent = 0;
        # #     right->parent = parent;
        # #     right->indexInParent = 1;
        # # }
        self.size = 2
        #         size(2),
        # left.parent = parent
        #     left->parent = parent;
        left.indexInParent = 0
        #     left->indexInParent = 0;
        # right.parent = parent
        #     right->parent = parent;
        right.indexInParent = 1
        #     right->indexInParent = 1;
        self.entries = [self.Entry().init_with_p(left), self.Entry().init_with_p(right), self.Entry(), self.Entry()]
        #         entries{Entry(left), Entry(right), Entry()} {
        return self
        # }
        pass  # end of def init_with_arguments

#     /**
#      * When an internal node is dropped, clear the entries it points to
#      */
#     ~LInternalNode() {
#         for (auto &entry : entries) {
#             entry.remove();
#         }
#     }

#     /**
#      * Returns the total number of bits in this tree, by summing up the
#      * b-parts of the entries
#      */
    def bits(self):
        pass
        #     unsigned long bits();
        # # unsigned long LInternalNode::bits() {
        total: int = 0
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
        pass  # end of def bits

#     /**
#      * Takes the leftmost entry out of this node and returns it
#      * @return
#      */
    def popFirst(self):
        pass
        #     Entry popFirst();
        # # LInternalNode::Entry LInternalNode::popFirst() {
        # #     // Take the first entry out, move everything else left
        result: LInternalNode.Entry = self.entries[0]
        # #     LInternalNode::Entry result = this->entries[0];
        self.remove(0)
        # #     this->remove(0);
        return result
        # #     return result;
        # # }
        pass

#     /**
#      * Takes the rightmost entry out of this node and returns it
#      * @return
#      */
    def popLast(self):
        pass
        #     Entry popLast();
        # # LInternalNode::Entry LInternalNode::popLast() {
        # #     // Take the last entry out
        self.size -= 1
        # #     size--;
        result: LInternalNode.Entry = self.entries[self.size]
        # #     LInternalNode::Entry result = entries[size];
        self.entries[self.size] = LInternalNode.Entry()
        # #     entries[size] = Entry();
        return result
        # #     return result;
        # # }
        pass

#     /**
#      * Adds the given entry to this node at the specified position
#      */
    def insert(self, idx: int, entry: LInternalNode.Entry):
        pass
        try:
            if idx < 0 or idx > 3:
                pass
            #     void insert(unsigned long, Entry);
            # # void LInternalNode::insert(unsigned long idx, LInternalNode::Entry entry) {
            # #     // Move everything from idx onwards right
            i: int = self.size
            while i > idx:
                # #     for (unsigned long i = size; i > idx; i--) {
                if i < 1 or i > 3:
                    pass
                self.entries[i] = self.entries[i - 1]
                # #         entries[i] = entries[i - 1];
                self.entries[i].p.indexInParent = i
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
        except Exception as e:
            pass
        pass

#     /**
#      * Adds the given entry to the end of this node
#      */
    def append(self, entry: LInternalNode.Entry):
        pass
        #     void append(Entry);
        # # void LInternalNode::append(LInternalNode::Entry entry) {
        entry.p.indexInParent = self.size
        # #     entry.P->indexInParent = size;
        self.entries[self.size] = entry
        # #     entries[size] = entry;
        self.size += 1
        # #     size++;
        # # }
        pass

#     /**
#      * Removes the entry at the specified position from this node
#      */
    def remove(self, idx: int):
        pass
        #     void remove(unsigned long);
        # # void LInternalNode::remove(unsigned long idx) {
        self.size -= 1
        # #     size--;
        for i in range(idx, self.size):
            # #     for (unsigned long i = idx; i < size; i++) {
            self.entries[i] = self.entries[i + 1]
            # #         entries[i] = entries[i + 1];
            self.entries[i].p.indexInParent = i
            # #         entries[i].P->indexInParent = i;
            # #     }
        self.entries[self.size] = LInternalNode.Entry()
        # #     entries[size] = Entry();
        # # }
        pass
    # };
    pass  # end of class LInternalNode


class LLeafNode:
    pass
    # /** A leaf node, which consists of a BitVector<> (represented by vector<bool>) */
    # struct LLeafNode {
    #   BitVector<> bv;

    #   /**
    #    * Constructs a leaf with the given number of bits
    #    */
    #   explicit LLeafNode(unsigned long size) :
    def __init__(self):
        pass
        self.bv = BitVector()
        # self.bits = 0
        pass

    def init_with_size(self, size: int):
        self.bv.init_with_size(size)
        #       bv(size) {}
        return self
        pass

    def init_with_bitvector(self, bv):
        pass
        #   /**
        #    * Constructs a leaf node from the given bit vector
        #    * @param bv the BitVector<> to be moved into this leaf node
        #    */
        #   explicit LLeafNode(BitVector<> bv) :
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
        # # unsigned long LLeafNode::bits() {
        return self.bv.size()
        # #     return bv.size();
        # # }
        pass
    # };
    pass  # end of class LLeafNode


class LTree:
    pass
    # /** LTree is the struct representing a single node (leaf or internal) of the LTree */
    # /** A single node is either an internal or leaf node, as indicated by the isLeaf value */
    # struct LTree {

    def __init__(self):
        pass
        self.isLeaf: bool = True
        #     bool isLeaf;
        self.parent: Optional[LTree] = None
        #     LTree *parent = nullptr;
        self.indexInParent: int = 0
        #     unsigned long indexInParent = 0;
        self.node = LTree.Node()
        pass

    class Node:
        # # // Constructors and destructors for data types that can't be in LTree.h
        def __init__(self):
            # # LTree::Node::Node() {
            self.internalNode = None
            # #     this->internalNode = nullptr;
            self.leafNode = LLeafNode().init_with_size(0)
            # #     this->leafNode = new LLeafNode(0);
            # # }
            pass

        def init_with_ps(self, p1: LTree, p2: LTree):
            # # LTree::Node::Node(LTree *P1, LTree *P2) {
            self.leafNode = None
            # #     this->leafNode = nullptr;
            self.internalNode = LInternalNode().init_with_arguments(p1, p2)
            # #     this->internalNode = new LInternalNode(P1, P2);
            # # }
            return self
            pass

        def init_with_bv(self, bv: BitVector):
            # # LTree::Node::Node(BitVector<> bv) {
            self.internalNode = None
            # #     this->internalNode = nullptr;
            self.leafNode = LLeafNode().init_with_bitvector(bv)
            # #     this->leafNode = new LLeafNode(bv);
            # # }
            return self
            pass

        def init_with_size(self, size: int):
            # # LTree::Node::Node(unsigned long size) {
            self.internalNode = None
            # #     this->internalNode = nullptr;
            self.leafNode = LLeafNode().init_with_size(size)
            # #     this->leafNode = new LLeafNode(size);
            # # }
            return self
            pass

            #     union Node {
            #         LInternalNode *internalNode;
            #         LLeafNode *leafNode;
            #
            #         Node();
            #
            #         Node(LTree *, LTree *);
            #
            #         explicit Node(BitVector<>);
            #
            #         explicit Node(unsigned long);
            #     } node;
            #
    #     /**
    #      * Constructs an empty leaf node
    #      */

    #     LTree() :
    #             isLeaf(true),
    #             node() {}

    def init_with_arguments(self, left: LTree, right: LTree) -> LTree:
        #     /**
        #      * Constructs a node with the two given `LTree`s as children
        #      * @param left the first child of this node
        #      * @param right the second child of this node
        #      */
        #     LTree(LTree *left, LTree *right) :
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
        self.node = self.Node().init_with_ps(left, right)
        #             node(left, right) {
        #     }
        return self
        pass

    def init_with_bv(self, bv: BitVector):
        #     /**
        #      * Constructs a leaf node with the given bit vector
        #      * @param bv the bit vector to be moved into this leaf node
        #      */
        #     explicit LTree(BitVector<> bv) :
        self.isLeaf = True
        #             isLeaf(true),
        self.node = LTree.Node().init_with_bv(bv)
        #             node(bv) {}
        return self
        pass

    def init_with_size(self, size: int) -> LTree:
        #     /**
        #      * Constructs an all-zeros leaf node with the specified size
        #      * @param size the size of this leaf node in bits
        #      */
        #     explicit LTree(unsigned long size) :
        self.isLeaf = True
        #             isLeaf(true),
        self.node = LTree.Node().init_with_size(size)
        #             node(size) {}
        return self
        pass

    def __del__(self):
        #     /// The LTree destructor decides which variant of the union to destroy
        #     ~LTree();
        # # LTree::~LTree() {
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

    def depth(self) -> int:
        pass
        #     /**
        #      * Returns the depth of this node, which is the length of the path from this
        #      * node to the root
        #      * @return the depth of `this`
        #      */
        #     unsigned long depth();
        # # unsigned long LTree::depth() {
        if self.parent is None:
            # #     if (parent == nullptr) {
            return 0
            # #         return 0;
        else:
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
        # # unsigned long LTree::height() {
        if self.isLeaf:
            # #     if (isLeaf) {
            return 0
            # #         return 0;
        else:
            # #     } else {
            entries = self.node.internalNode.entries
            # #         auto &entries = node.internalNode->entries;
            max_: int = 0
            # #         unsigned long max = 0;
            for entry in entries:
                # #         for (auto &entry : entries) {
                depth: int = entry.p.height()
                # #             unsigned long depth = entry.P->height();
                if depth > max_:
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

    def size(self):
        pass
        #     /**
        #      * Gets the size of this node, either as the number of children (internal
        #      * node) or in number of k^2 blocks (leaf node)
        #      * @return the size of this node in number of children or bit blocks,
        #      *         depending on the type of the node
        #      */
        #     unsigned long size();
        # # unsigned long LTree::size() {
        if self.isLeaf:
            # #     if (isLeaf) {
            return int(self.node.leafNode.bits() / Parameters.BLOCK_SIZE)
            # #         return node.leafNode->bits() / BLOCK_SIZE;
        else:
            # #     } else {
            return self.node.internalNode.size
            # #         return node.internalNode->size;
            # #     }
        # # }
        pass

    def findChild(self, n: int) -> LRecord:
        pass
        #     /**
        #      * Given an integer n, returns the child node containing the n-th bit in this subtree,
        #      *      as well as the numbers of bits preceding it.
        #      * @param n  an integer that satisfies 0 <= n < (number of bits in this subtree)
        #      * @return  a LTree::LInternalNode::Entry with:
        #      *      b = number of bits preceding this node in the LTree
        #      *      i = the index of the relevant subtree in the `entries` of `this`
        #      */
        #     LRecord findChild(unsigned long);
        # # LRecord LTree::findChild(unsigned long n) {
        bitsBefore: int = 0
        # #     unsigned long bitsBefore = 0;
        # #
        node: LInternalNode = self.node.internalNode
        # #     LInternalNode *node = this->node.internalNode;
        # #     unsigned long i;
        ii: int = 0
        for i in range(node.size):
            ii = i
            # #     for (i = 0; i < node->size; i++) {
            entry = node.entries[i]
            # #         auto &entry = node->entries[i];
            if bitsBefore + entry.b > n:
                # #         if (bitsBefore + entry.b > n) {
                return LRecord(bitsBefore, i)
                # #             return {bitsBefore, i};
                # #         }
            bitsBefore += entry.b
            # #         bitsBefore += entry.b;
            # #     }
        ii += 1
        # #     // If the required bit is one after the last bit in this tree,
        # #     // return the last child anyway
        # #     // This is necessary for appending bits
        if ii == node.size and bitsBefore == n:
            # #     if (i == node->size && bitsBefore == n) {
            entry = node.entries[ii - 1]
            # #         auto &entry = node->entries[i - 1];
            return LRecord(bitsBefore - entry.b, ii - 1)
            # #         return {bitsBefore - entry.b, i - 1};
            # #     }
        # #     // If we reach this point, that means that the size of this subtree is less than n
        # #     // so the input parameter was out of range
        raise "LTree: index out of range"
        # #     throw std::range_error("LTree: index out of range");
        # # }
        pass

    def findLeaf(self, n: int, path) -> LInternalNode.Entry:
        pass
        #     /**
        #      * Given an integer n, returns the leaf containing the n-th bit of the BitVector<>,
        #      *      as well as the numbers of bits preceding this leaf node.
        #      * @param n  an integer that satisfies 0 <= n < (number of bits in this subtree)
        #      * @param path  a pointer to a vector of `LNesbo` entries that represent the last path took
        #      *        if the path is nullptr, a regular `findLeaf` is done. If the path is empty, a regular
        #      *        `findLeaf` is done but the result is stored in the path vector
        #      * @return  a LTree::LInternalNode::Entry with:
        #      *      b = number of bits preceding this node in the LTree
        #      *      P = a pointer to the leaf node containing the n-th bit
        #      * Note that P is strictly just a pointer to a `LTreeNode`, as defined by the `entry` struct
        #      * But the union is always of the `LLeafNode` variant.
        #      */
        #     LInternalNode::Entry findLeaf(unsigned long, vector<LNesbo> *path = nullptr);
        # # LInternalNode::Entry LTree::findLeaf(unsigned long n, vector<LNesbo> *path) {
        if path is None:
            # #     if (path == nullptr) {
            current = self
            # #         auto *current = this;
            bitsBefore: int = 0
            # #         unsigned long bitsBefore = 0;
            while not current.isLeaf:
                # #         while (!current->isLeaf) {
                record: LRecord = current.findChild(n - bitsBefore)
                # #             auto record = current->findChild(n - bitsBefore);
                bitsBefore += record.b
                # #             bitsBefore += record.b;
                current = current.node.internalNode.entries[record.i].p
                # #             current = current->node.internalNode->entries[record.i].P;
                # #         }
            return LInternalNode.Entry().init_with_arguments(bitsBefore, current)
            # #         return {bitsBefore, current};
        else:
            # #     } else {
            return self.findLeaf2(n, path)
            # #         return findLeaf2(n, *path);
            # #     }
        # # }
        pass

    def findLeaf2(self, n: int, path) -> LInternalNode.Entry:
        pass
        #     LInternalNode::Entry findLeaf2(unsigned long, vector<LNesbo> &);
        # # LInternalNode::Entry LTree::findLeaf2(unsigned long n, vector<LNesbo> &path) {
        current = None
        # #     LTree *current = nullptr;
        bitsBefore: int = 0
        # #     unsigned long bitsBefore = 0;
        if len(path) == 0:
            # #     if (path.empty()) {
            # #         // If the path is empty, we have to do a regular findLeaf and store the
            # #         // results in the path vector
            current = self
            # #         current = this;
            bitsBefore = 0
            # #         bitsBefore = 0;
        else:
            # #     } else {
            # #         // Start at the end of the path (which is a leaf), and move up until we
            # #         // reach the subtree that also contains the desired bit
            k: int = len(path) - 1
            # #         unsigned long k = path.size() - 1;
            nesbo: LNesbo = path[k]
            # #         LNesbo nesbo = path[k];
            # #         // While the subtree nested at `nesbo.node` is entirely before, or
            # #         // entirely after the bit we are looking for...
            while (n < nesbo.bits_before) or nesbo.bits_before + nesbo.size <= n:
                # #         while (n < nesbo.bitsBefore || nesbo.bitsBefore + nesbo.size <= n) {
                # #             // Go `up` one level if possible. Else, we are at the root so do a
                # #             // regular search from there
                path.pop()
                # #             path.pop_back();
                if k == 0:
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
                # #         if (current == nullptr) {
                current = nesbo.node.node.internalNode.entries[nesbo.index].p
                # #             current = nesbo.node->node.internalNode->entries[nesbo.index].P;
                bitsBefore = nesbo.bits_before
                # #             bitsBefore = nesbo.bitsBefore;
                # #         }
            # #     }
        # #     // Now, we do a regular search from `current` and add the intermediate
        # #     // nodes we visit to `path`
        while not current.isLeaf:
            # #     while (!current->isLeaf) {
            record = current.findChild(n - bitsBefore)
            # #         auto record = current->findChild(n - bitsBefore);
            bitsBefore += record.b
            # #         bitsBefore += record.b;
            next_ = current.node.internalNode.entries[record.i]
            # #         auto next = current->node.internalNode->entries[record.i];
            xxx = LNesbo().init_with_arguments(current, record.i, next_.b, bitsBefore)
            path.append(xxx)
            # #         path.emplace_back(current, record.i, next.b, bitsBefore);
            current = next_.p
            # #         current = next.P;
            # #     }
        return LInternalNode.Entry().init_with_arguments(bitsBefore, current)
        # #     return {bitsBefore, current};
        # # }
        pass

    def access(self, n: int, path: list[LNesbo]):
        pass
        #     /**
        #      * Performs the `access` operation on this subtree
        #      * @param n the index of a bit in the `LTree`
        #      * @return the value of the `n`-th bit in the tree
        #      */
        #     bool access(unsigned long, vector<LNesbo> *path = nullptr);
        # # bool LTree::access(unsigned long n, vector<LNesbo> *path) {
        entry = self.findLeaf(n, path)
        # #     auto entry = findLeaf(n, path);
        return entry.p.node.leafNode.bv[n - entry.b]
        # #     return entry.P->node.leafNode->bv[n - entry.b];
        # # }
        pass

    def set_bit(self, n: int, b: bool, path: list[LNesbo]):
        pass
        #     /**
        #      * Sets the bit at position n to the value of b
        #      * @param n the index of the bit to set
        #      * @param b the value to set the bit to
        #      * @return true if this bit changed, e.g.
        #      *      if the previous value of the bit was unequal to b
        #      */
        #     bool setBit(unsigned long, bool, vector<LNesbo> *path = nullptr);
        # # bool LTree::setBit(unsigned long n, bool b, vector<LNesbo> *path) {
        # #     // Find the leaf node that contains this bit
        entry = self.findLeaf(n, path)
        # #     auto entry = findLeaf(n, path);
        bv = entry.p.node.leafNode.bv
        # #     BitVector<> &bv = entry.P->node.leafNode->bv;
        changed = bv.set(n - entry.b, b)
        # #     bool changed = bv.set(n - entry.b, b);
        # #
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
        # # unsigned long LTree::bits() {
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

    def updateCounters(self, dBits: int):
        pass
        #     /**
        #      * Changes the values of the counter `b` of all the nodes whose
        #      * subtree contains this node. Used to update these counters after modifying
        #      * the underlying BitVector<>, or the structure of the tree.
        #      *
        #      * @param dBits the change in the number of bits (e.g. 4 when 4 bits were inserted)
        #      */
        #     void updateCounters(long);
        # # void LTree::updateCounters(long dBits) {
        current: LTree = self
        # #     LTree *current = this;
        # #     // Go up in the tree until we reach the root
        while current.parent is not None:
            # #     while (current->parent != nullptr) {
            # #         // Take the entry in `current`s parent that points to `current`,
            # #         // and update its `b` counter.
            parent = current.parent
            # #         auto parent = current->parent;
            entry = parent.node.internalNode.entries[current.indexInParent]
            # #         auto &entry = parent->node.internalNode->entries[current->indexInParent];
            entry.b += dBits
            # #         entry.b += dBits;
            # #
            current = parent
            # #         current = parent;
        # #     }
        # # }
        pass

    def insertBlock(self, index: int, path: list[LNesbo]):
        pass
        #     /**
        #      * Inserts one block of k^2 bits at the specified position
        #      *
        #      * @return the new root if the tree's root changed, nullptr otherwise
        #      */
        #     LTree *insertBlock(long unsigned, vector<LNesbo> *path = nullptr);
        # # LTree *LTree::insertBlock(long unsigned index, vector<LNesbo> *path) {
        return self.insertBits(index, Parameters.BLOCK_SIZE, path)
        # #     return this->insertBits(index, BLOCK_SIZE, path);
        # # }
        pass

    def deleteBlock(self, index: int, path: list[LNesbo]):
        #     /**
        #      * Deletes a block of k^2 bits at the specified position
        #      *
        #      * @return the new root if the tree's root changed, nullptr otherwise
        #      */
        #     LTree *deleteBlock(long unsigned, vector<LNesbo> *path = nullptr);
        # # LTree *LTree::deleteBlock(long unsigned index, vector<LNesbo> *path) {
        return self.deleteBits(index, Parameters.BLOCK_SIZE, path)
        # #     return this->deleteBits(index, BLOCK_SIZE, path);
        # # }
        pass

    def memoryUsage(self) -> int:
        pass
        #   unsigned long memoryUsage();
        #   unsigned long LTree::memoryUsage() {
        result: int = sys.getsizeof(LTree)
        #       unsigned long result = sizeof(LTree);
        if self.isLeaf:
            #       if (isLeaf) {
            result += self.node.leafNode.bv.memory_usage()
            #           result += node.leafNode->bv.memoryUsage();
        else:
            #       } else {
            result += sys.getsizeof(LInternalNode)
            #           result += sizeof(LInternalNode);
            entries = self.node.internalNode.entries
            #           auto &entries = node.internalNode->entries;
            for entry in entries:
                #           for (auto &entry : entries) {
                if entry.p is not None:
                    #               if (entry.P != nullptr) {
                    result += entry.p.memoryUsage()
                    #                   result += entry.P->memoryUsage();
                    #               }
                #           }
            #       }
        return result
        #       return result;
        #   }
        pass

    def insertBits(self, index: int, count: int, path: list[LNesbo]):
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
        #      */

        #     LTree *
        #     insertBits(long unsigned, long unsigned, vector<LNesbo> *path = nullptr);
        # # LTree *LTree::insertBits(long unsigned index, long unsigned count,
        # #                          vector<LNesbo> *path) {
        entry = self.findLeaf(index, path)
        # #     auto entry = findLeaf(index, path);
        leaf = entry.p
        # #     auto leaf = entry.P;
        bv = leaf.node.leafNode.bv
        # #     auto &bv = leaf->node.leafNode->bv;
        bv.insert(index - entry.b, count)
        # #     bv.insert(index - entry.b, count);
        leaf.updateCounters(count)
        # #     leaf->updateCounters(count);
        # #
        # #     // Split this node up into two if it exceeds the size limit
        return leaf.checkSizeUpper()
        # #     return leaf->checkSizeUpper();
        # # }
        pass

    def deleteBits(self, index: int, count: int, path: list[LNesbo]) -> LTree:
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
        #     LTree *
        #     deleteBits(long unsigned, long unsigned, vector<LNesbo> *path = nullptr);
        # # LTree *LTree::deleteBits(long unsigned index, long unsigned count,
        # #                          vector<LNesbo> *path) {
        entry = self.findLeaf(index, path)
        # #     auto entry = findLeaf(index, path);
        leaf = entry.p
        # #     auto leaf = entry.P;
        bv = leaf.node.leafNode.bv
        # #     auto &bv = leaf->node.leafNode->bv;
        start = index - entry.b
        # #     long unsigned start = index - entry.b;
        end = start + count
        # #     long unsigned end = start + count;
        bv.erase(start, end)
        # #     bv.erase(start, end);
        leaf.updateCounters(-count)
        # #     leaf->updateCounters(-count);
        return leaf.checkSizeLower()
        # #     return leaf->checkSizeLower();
        # # }
        pass

    def checkSizeUpper(self):
        pass
        #     /**
        #      * Checks if this node satisfies the maximum size for an internal node
        #      * or leaf node. If not, tries to spill a node to a sibling, or if that
        #      * fails will split this node into two and recursively check the parent
        #      *
        #      * @return nullptr in most cases, but returns the new root if it has
        #      *         changed, e.g. if the height of the tree has increased
        #      */
        #     LTree *checkSizeUpper();
        # # LTree *LTree::checkSizeUpper() {
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

    def checkSizeLower(self):
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
        #     LTree *checkSizeLower();
        # # LTree *LTree::checkSizeLower() {
        isRoot = not self.parent
        # #     bool isRoot = (parent == nullptr);
        if isRoot and (self.size() >= 2 or self.isLeaf):
            # #     if (isRoot && (size() >= 2 || isLeaf)) {
            return None
            # #         return nullptr;
        elif self.isLeaf and self.size() < Parameters.leafSizeMin:
            # #     } else if (isLeaf && size() < leafSizeMin) {
            if not self.tryStealLeaf():
                # #         if (!tryStealLeaf()) {
                return self.merge_leaf()
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

    def trySpillInternal(self):
        pass
        #     /**
        #      * Tries to move a child of an internal node to a sibling, and returns
        #      * whether it succeeded
        #      *
        #      * @return true if a spill could be done, false if not (e.g. if both of
        #      *         the nodes siblings are of maximum size)
        #      */
        #     bool trySpillInternal();
        # # bool LTree::trySpillInternal() {
        if not self.parent:
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

    def trySpillLeaf(self):
        pass
        #     /**
        #      * Tries to move a child of a leaf node to a sibling, and returns
        #      * whether it succeeded
        #      *
        #      * @return true if a spill could be done, false if not (e.g. if both of
        #      *         the nodes siblings are of maximum size)
        #      */
        #     bool trySpillLeaf();
        # # bool LTree::trySpillLeaf() {
        if not self.parent:
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

    def splitInternal(self):
        pass
        #     /**
        #      * Splits this node into two nodes of minimum size, and recursively
        #      * checks the rest of the tree for meeting size requirements
        #      *
        #      * @return nullptr in most cases, but returns the new root if this operation
        #      *         causes the tree's height to increase, which changes the root
        #      */
        #     LTree *splitInternal();
        # # LTree *LTree::splitInternal() {
        entries = self.node.internalNode.entries
        # #     auto &entries = this->node.internalNode->entries;
        n = self.size()
        # #     unsigned long n = this->size();
        mid = int(n / 2)
        # #     unsigned long mid = n / 2;
        newNode = LTree()
        # #     auto newNode = new LTree();
        newNode.parent = self.parent
        # #     newNode->parent = parent;
        newNode.isLeaf = False
        # #     newNode->isLeaf = false;
        newNode.node.leafNode = None
        # #     newNode->node.leafNode = nullptr;
        newNode.node.internalNode = LInternalNode()
        # #     newNode->node.internalNode = new LInternalNode();
        d_b = 0
        # #     unsigned long d_b = 0; // Count bits in right half
        for i in range(mid, n):
            # #     for (unsigned long i = mid; i < n; i++) {
            entry = entries[i]
            # #         auto entry = entries[i];
            if entry.p:
                # #         if (entry.P != nullptr) {
                entry.p.parent = newNode
                # #             entry.P->parent = newNode;
                entry.p.indexInParent = i - mid
                # #             entry.P->indexInParent = i - mid;
                d_b += entry.b
                # #             d_b += entry.b;
                # #         }
            newNode.node.internalNode.entries[i - mid] = entry
            # #         newNode->node.internalNode->entries[i - mid] = entry;
            entries[i] = LInternalNode.Entry()
            # #         entries[i] = LInternalNode::Entry();
            # #     }
        newNode.node.internalNode.size = n - mid
        # #     newNode->node.internalNode->size = n - mid;
        self.node.internalNode.size = mid
        # #     node.internalNode->size = mid;
        if self.parent is None:
            # #     if (parent == nullptr) {
            newRoot = LTree().init_with_arguments(self, newNode)
            # #         auto *newRoot = new LTree(this, newNode);
            return newRoot
            # #         return newRoot;
        else:
            # #     } else {
            self.parent.node.internalNode.entries[self.indexInParent].b -= d_b
            # #         parent->node.internalNode->entries[indexInParent].b -= d_b;
            entry = LInternalNode.Entry().init_with_arguments(d_b, newNode)
            self.parent.node.internalNode.insert(self.indexInParent + 1, entry)
            # #         parent->node.internalNode->insert(indexInParent + 1,
            # #                                           {d_b, newNode});
            return self.parent.checkSizeUpper()
            # #         return parent->checkSizeUpper();
            # #     }
        # # }
        pass

    def splitLeaf(self):
        pass
        #     /**
        #      * Splits this node into two nodes of minimum size, and recursively
        #      * checks the rest of the tree for meeting size requirements
        #      *
        #      * @return nullptr in most cases, but returns the new root if this operation
        #      *         causes the tree's height to increase, which changes the root
        #      */
        #     LTree *splitLeaf();
        # # LTree *LTree::splitLeaf() {
        n = self.node.leafNode.bits()
        # #     unsigned long n = this->node.leafNode->bits();
        mid = int(n / 2)
        # #     unsigned long mid = n / 2;
        mid -= mid % Parameters.BLOCK_SIZE
        # #     mid -= mid % BLOCK_SIZE;
        left = self.node.leafNode.bv
        # #     auto &left = this->node.leafNode->bv;
        right = BitVector().init_with_range(left, mid, n)
        # #     auto right = BitVector<>(left, mid, n);
        left.erase(mid, n)
        # #     left.erase(mid, n);
        newNode = LTree().init_with_bv(right)
        # #     auto *newNode = new LTree(right);
        if self.parent is None:
            # #     if (parent == nullptr) {
            newRoot = LTree().init_with_arguments(self, newNode)
            # #         auto *newRoot = new LTree(this, newNode);
            return newRoot
            # #         return newRoot;
        else:
            # #     } else {
            idx = self.indexInParent
            # #         unsigned long idx = indexInParent;
            newNode.parent = self.parent
            # #         newNode->parent = parent;
            entry = LInternalNode.Entry().init_with_p(newNode)
            # #         LInternalNode::Entry entry(newNode);
            self.parent.node.internalNode.insert(self.indexInParent + 1, entry)
            # #         parent->node.internalNode->insert(indexInParent + 1, entry);
            self.parent.node.internalNode.entries[idx].b -= entry.b
            # #         parent->node.internalNode->entries[idx].b -= entry.b;
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
        # # bool LTree::tryStealInternal() {
        if not self.parent:
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
        # # bool LTree::tryStealLeaf() {
        if not self.parent:
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
        try:
            #     /**
            #      * Merges this node with a sibling, and recursively checks the rest of the
            #      * tree for meeting size constraints
            #      *
            #      * @return nullptr usually, but returns the new root if it changed due to this
            #      *         operation, e.g. when the height of the tree changed
            #      */
            #     LTree *mergeInternal();
            # # LTree *LTree::mergeInternal() {
            # #     // If we are the root, and we are too small, then we have only one child
            if not self.parent:
                # #     if (parent == nullptr) {
                # #         // Delete this, our only child should become the root
                child = self.node.internalNode.entries[0].p
                # #         LTree *child = node.internalNode->entries[0].P;
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
            idx = self.indexInParent
            # #     unsigned long idx = indexInParent;
            # left = None
            # right = None
            # #     LTree *left = nullptr, *right = nullptr;
            if idx > 0:
                # #     if (idx > 0) {
                left = self.parent.node.internalNode.entries[idx - 1].p
                # #         left = parent->node.internalNode->entries[idx - 1].P;
                right = self
                # #         right = this;
                idx -= 1
                # #         idx--;
            else:
                # #     } else {
                left = self
                # #         left = this;
                right = self.parent.node.internalNode.entries[idx + 1].p
                # #         right = parent->node.internalNode->entries[idx + 1].P;
                # #     }
            # #
            # #     // Merge `left` and `right` into one node
            internalNode = left.node.internalNode
            # #     auto &internalNode = left->node.internalNode;
            n = right.size()
            # #     unsigned long n = right->size();
            d_b = 0
            # #     unsigned long d_b = 0;
            for i in range(n):
                # #     for (unsigned i = 0; i < n; i++) {
                entry: LInternalNode.Entry = copy.deepcopy(right.node.internalNode.entries[i])  # 2024/8/23
                # entry = right.node.internalNode.entries[i]
                # #         auto entry = right->node.internalNode->entries[i];
                right.node.internalNode.entries[i].p = None
                # #         right->node.internalNode->entries[i].P = nullptr;
                d_b += entry.b
                # #         d_b += entry.b;
                entry.p.parent = left
                # #         entry.P->parent = left;
                internalNode.append(entry)
                # #         internalNode->append(entry);
                # #     }
            # #     // Delete the right child, and update the b counter for left
            self.parent.node.internalNode.remove(idx + 1)
            # #     parent->node.internalNode->remove(idx + 1);
            self.parent.node.internalNode.entries[idx].b += d_b
            # #     parent->node.internalNode->entries[idx].b += d_b;
            return self.parent.checkSizeLower()
            # #     return parent->checkSizeLower();
            # # }
        except Exception as e:
            pass
        pass

    def merge_leaf(self):
        pass
        #     /**
        #      * Merges this node with a sibling, and recursively checks the rest of the
        #      * tree for meeting size constraints
        #      *
        #      * @return nullptr usually, but returns the new root if it changed due to this
        #      *         operation, e.g. when the height of the tree changed
        #      */
        #     LTree *mergeLeaf();
        # # LTree *LTree::mergeLeaf() {
        if not self.parent:
            # #     if (parent == nullptr) {
            return False
            # #         return nullptr;
            # #     }
        idx: int = self.indexInParent
        # #     unsigned long idx = indexInParent;
        # left = None
        # right = None
        # #     LTree *left = nullptr, *right = nullptr;
        if idx > 0:
            # #     if (idx > 0) {
            left = self.parent.node.internalNode.entries[idx - 1].p
            # #         left = parent->node.internalNode->entries[idx - 1].P;
            right = self
            # #         right = this;
            idx -= 1
            # #         idx--;
        else:
            # #     } else {
            left = self
            # #         left = this;
            right = self.parent.node.internalNode.entries[idx + 1].p
            # #         right = parent->node.internalNode->entries[idx + 1].P;
            # #     }
        leftBits = left.node.leafNode.bv
        # #     auto &leftBits = left->node.leafNode->bv;
        rightBits = right.node.leafNode.bv
        # #     auto &rightBits = right->node.leafNode->bv;
        # #     // Append `right`s bits to `left`
        leftBits.append(rightBits, 0, rightBits.size())
        # #     leftBits.append(rightBits, 0, rightBits.size());
        # #     // Update the b for `left`, and delete `right`
        d_b = self.parent.node.internalNode.entries[idx + 1].b
        # #     unsigned long d_b = parent->node.internalNode->entries[idx + 1].b;
        self.parent.node.internalNode.entries[idx].b += d_b
        # #     parent->node.internalNode->entries[idx].b += d_b;
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
        # # void LTree::moveLeftInternal() {
        # parent = self.parent
        # #     LTree *parent = this->parent;
        idx = self.indexInParent
        # #     unsigned long idx = this->indexInParent;
        sibling = self.parent.node.internalNode.entries[idx - 1].p
        # #     LTree *sibling = parent->node.internalNode->entries[idx - 1].P;

        # #     // Move the first child of `this` to the end of the left sibling
        toMove = self.node.internalNode.popFirst()
        # #     LInternalNode::Entry toMove = this->node.internalNode->popFirst();
        d_b = toMove.b
        # #     unsigned long d_b = toMove.b;
        toMove.p.parent = sibling
        # #     toMove.P->parent = sibling;
        sibling.node.internalNode.append(toMove)
        # #     sibling->node.internalNode->append(toMove);

        # #     // Finally, update the parent's b counter for `this` and `sibling`
        # #     // The number of bits in `toMove` is subtracted from `this`, but added to `sibling`
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx - 1].b += d_b
        # #     parent->node.internalNode->entries[idx - 1].b += d_b;
        # # }
        pass

    def moveRightInternal(self):
        pass
        #     /**
        #      * Moves the rightmost child of this node to the start of the right sibling
        #      */
        #     void moveRightInternal();
        # # void LTree::moveRightInternal() {
        idx = self.indexInParent
        # #     unsigned long idx = this->indexInParent;
        sibling = self.parent.node.internalNode.entries[idx + 1].p
        # #     LTree *sibling = parent->node.internalNode->entries[idx + 1].P;

        # #     // Move the last child of `this` to the start of the left sibling
        toMove = self.node.internalNode.popLast()
        # #     LInternalNode::Entry toMove = this->node.internalNode->popLast();
        d_b = toMove.b
        # #     unsigned long d_b = toMove.b;
        toMove.p.parent = sibling
        # #     toMove.P->parent = sibling;
        sibling.node.internalNode.insert(0, toMove)
        # #     sibling->node.internalNode->insert(0, toMove);

        # #     // Finally, update the parent's b counter for `this` and `sibling`
        # #     // The number of bits in `toMove` is subtracted from `this`, but added to `sibling`
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx + 1].b += d_b
        # #     parent->node.internalNode->entries[idx + 1].b += d_b;
        # # }
        pass

    def moveLeftLeaf(self):
        pass
        #     /**
        #      * Moves the leftmost block of k^2 bits to the end of the left sibling
        #      */
        #     void moveLeftLeaf();
        # # void LTree::moveLeftLeaf() {
        idx = self.indexInParent
        # #     unsigned long idx = indexInParent;
        sibling = self.parent.node.internalNode.entries[idx - 1].p
        # #     LTree *sibling = parent->node.internalNode->entries[idx - 1].P;
        # #     // Take the first k*k block of `this`, and append it to `sibling`
        right = self.node.leafNode.bv
        # #     BitVector<> &right = node.leafNode->bv;
        left = sibling.node.leafNode.bv
        # #     BitVector<> &left = sibling->node.leafNode->bv;
        d_b = Parameters.BLOCK_SIZE
        # #     unsigned long d_b = BLOCK_SIZE;
        left.append(right, 0, Parameters.BLOCK_SIZE)
        # #     left.append(right, 0, BLOCK_SIZE);
        right.erase(0, Parameters.BLOCK_SIZE)
        # #     right.erase(0, BLOCK_SIZE);

        # #     // Update the parent's b counter
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx - 1].b += d_b
        # #     parent->node.internalNode->entries[idx - 1].b += d_b;
        # # }
        pass

    def moveRightLeaf(self):
        pass
        #     /**
        #      * Moves the rightmost block of k^2 bits to the start of the right sibling
        #      */
        #     void moveRightLeaf();
        # # void LTree::moveRightLeaf() {
        idx = self.indexInParent
        # #     unsigned long idx = indexInParent;
        sibling = self.parent.node.internalNode.entries[idx + 1].p
        # #     LTree *sibling = parent->node.internalNode->entries[idx + 1].P;
        # #     // Take the first k*k block of `this`, and append it to `sibling`
        left = self.node.leafNode.bv
        # #     BitVector<> &left = node.leafNode->bv;
        right = sibling.node.leafNode.bv
        # #     BitVector<> &right = sibling->node.leafNode->bv;
        hi = left.size()
        # #     unsigned long hi = left.size();
        lo = hi - Parameters.BLOCK_SIZE
        # #     unsigned long lo = hi - BLOCK_SIZE;
        d_b = Parameters.BLOCK_SIZE
        # #     unsigned long d_b = BLOCK_SIZE;
        right.insert_range(0, left, lo, hi)
        # #     right.insert(0, left, lo, hi);
        left.erase(lo, hi)
        # #     left.erase(lo, hi);

        # #     // Update the parent's b counter
        self.parent.node.internalNode.entries[idx].b -= d_b
        # #     parent->node.internalNode->entries[idx].b -= d_b;
        self.parent.node.internalNode.entries[idx + 1].b += d_b
        # #     parent->node.internalNode->entries[idx + 1].b += d_b;
        # # }
        pass

# };
    pass  # end of class LTree

# #endif //DK2TREE_LTREE_H

# """
# LTree.cpp
# """
# //
# // Created by anneke on 18/12/18.
# //
#
# #include <iostream>
# #include <utility>
# #include "LTree.h"
