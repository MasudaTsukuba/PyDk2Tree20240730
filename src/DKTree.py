"""
DKTree.py << DKTree.cpp
2024/7/30, T. Masuda
Amagasa Laboratory, University of Tsukuba
"""
from __future__ import annotations
import sys

# //
# // Created by anneke on 05/02/19.
# //

# #include <iostream>
# #include <sstream>
# #include "DKTree.h"
# using namespace std;
from src.Parameters import Parameters
from src.TTree import TTree, Nesbo
from src.LTree import LTree, LNesbo


def long_pow(a: int, b: int):
    pass
    # unsigned long long_pow(unsigned long a, unsigned long b) {
    result: int = 1
    #     unsigned long result = 1;
    while b != 0:
        pass
        #     while (b != 0) {
        if (b % 2) != 0:
            pass
            #         if (b % 2 != 0) {
            result *= a
            #             result *= a;
            #         }
        b >>= 1
        #         b >>= 1;
        a *= a
        #         a *= a;
        #     }
    return result
#     return result;
# }
    pass


class VectorData:
    pass
    # // a class that contains a vector of entries in the matrix
    # // with a start and end such that the entries at index start <= i < end are all in the same block at when
    # // dividing the matrix size by k^(iteration-1), and the bit with offset 0 for this iteration is located at firstAt
    # class VectorData {
    # public:

    def __init__(self):
        pass
        self.entry: list[int] = []
        #     const vector<unsigned long> &entry;
        self.start: int = 0
        #     unsigned long start;
        self.end: int = 0
        #     unsigned long end;
        self.iteration: int = 1
        #     unsigned long iteration;
        self.firstAt: int = 0
        #     unsigned long firstAt;
        pass

    def init_with_aentry(self, a_entry: list[int]):
        pass
        #     // constructor for a first iteration, it includes the whole vector and always starts at the first bit of
        #     the ttree
        #     explicit VectorData(vector<unsigned long> &aEntry)
        #             : entry(aEntry), start(0), end(aEntry.size()), iteration(1), firstAt(0) {}
        self.entry: list[int] = a_entry
        self.end: int = len(a_entry)
        return self
        pass

    def init_with_arguments(self, vector_data: VectorData, a_start: int, a_end: int, a_iteration: int, a_first_at: int):
        pass
        #     VectorData(VectorData &vectorData, unsigned long aStart, unsigned long aEnd, unsigned long aIteration,
        #     unsigned long aFirstAt)
        #             : entry(vectorData.entry), start(aStart), end(aEnd), iteration(aIteration), firstAt(aFirstAt) {}
        # self = vector_data
        self.entry: list[int] = vector_data.entry
        self.start: int = a_start
        self.end: int = a_end
        self.iteration: int = a_iteration
        self.firstAt: int = a_first_at
        return self
        pass

        # };
    pass  # end of class VectorData


class DKTree:
    # DKTree::DKTree() : ttree(new TTree()), ltree(new LTree()), freeColumns(), firstFreeColumn(0),
    # matrixSize(long_pow(k, 4ul)) {
    def __init__(self):
        self.ttree: TTree = TTree()
        self.ltree: LTree = LTree()
        self.tPath: list[Nesbo] = []  # | None = None
        self.lPath: list[LNesbo] = []  # | None = None
        self.freeColumns = []
        self.firstFreeColumn: int = 0
        self.matrixSize: int = long_pow(Parameters.k, 4)
        self.ttree.insertBlock(0)
        #     ttree->insertBlock(0);
        # }

    def init_with_power(self, power: int):
        pass
        # DKTree::DKTree(unsigned long power) : ttree(new TTree()), ltree(new LTree()), freeColumns(),
        # firstFreeColumn(0), matrixSize(long_pow(k, power)) {
        self.matrixSize: int = long_pow(Parameters.k, power)
        self.ttree.insertBlock(0)
        #     ttree->insertBlock(0);
        # }
        pass

    def __del__(self):
        pass
        # DKTree::~DKTree() {
        del self.ttree
        #     delete ttree;
        del self.ltree
        #     delete ltree;
        # }
        pass

    def add_edge(self, row: int, column: int):
        pass
        # void DKTree::addEdge(unsigned long row, unsigned long column) {
        #     // tests if both positions exist
        self.checkArgument(row, "addEdge")
        #     checkArgument(row, "addEdge");
        self.checkArgument(column, "addEdge")
        #     checkArgument(column, "addEdge");

        iteration: int = 1
        #     unsigned long iteration = 1;
        #     // the first position is 0+offset of the first iteration
        position: int = self.calculateOffset(row, column, iteration)
        #     unsigned long position = calculateOffset(row, column, iteration);;
        c_entry: bool = True
        #     bool cEntry = true; // to get the loop started
        #     // will change the values of iteration, position, and cEntry
        iteration, position, c_entry = self.traverseToFirst0OrEndOfTTree(row, column, iteration, position, c_entry)
        #     traverseToFirst0OrEndOfTTree(row, column, iteration, position, cEntry);
        #     // if the last bit found in the ttree is a 1 then find the final result in the ltree and set it to 1
        if c_entry:
            #     if (cEntry) {
            ltree_position: int = position - self.ttree.bits()
            #         unsigned long ltreePosition = position - ttree->bits();
            self.ltree.set_bit(ltree_position, True, self.lPath)
            #         ltree->setBit(ltreePosition, true, &lPath);
        else:
            pass
            #     } else { // if not then change it to a 1 and insert new blocks where necessary
            self.ttree.set_bit(position, True, self.tPath)
            #         ttree->setBit(position, true, &tPath);
            iteration += 1
            #         iteration++;
            block_size: int = int(self.matrixSize / long_pow(Parameters.k, iteration))
            #         unsigned long blockSize = matrixSize / long_pow(k, iteration);
            while block_size > 1:
                #         while (blockSize > 1) {
                #             // position +1 since paper has rank including the position,
                #             but function is exclusive position
                insert_at: int = self.ttree.rank1(position + 1, self.tPath) * Parameters.BLOCK_SIZE
                #             unsigned long insertAt = ttree->rank1(position + 1, &tPath) * BLOCK_SIZE;
                self.insertBlockTtree(insert_at)
                #             insertBlockTtree(insertAt);
                offset: int = self.calculateOffset(row, column, iteration)
                #             unsigned long offset = calculateOffset(row, column, iteration);
                position: int = insert_at + offset
                #             position = insertAt + offset;
                self.ttree.set_bit(position, True, self.tPath)
                #             ttree->setBit(position, true, &tPath);
                iteration += 1
                #             iteration++;
                block_size: int = int(self.matrixSize / long_pow(Parameters.k, iteration))
                #             blockSize = matrixSize / long_pow(k, iteration);
                #         }
            #         // position +1 since paper has rank including the position, but function is exclusive position
            l_tree_insert_at: int = ((self.ttree.rank1(position + 1, self.tPath) * Parameters.BLOCK_SIZE)
                                     - self.ttree.bits())
            #         unsigned long lTreeInsertAt = (ttree->rank1(position + 1, &tPath) * BLOCK_SIZE) - ttree->bits();
            self.insertBlockLtree(l_tree_insert_at)
            #         insertBlockLtree(lTreeInsertAt);
            offset: int = self.calculateOffset(row, column, iteration)
            #         unsigned long offset = calculateOffset(row, column, iteration);
            position: int = l_tree_insert_at + offset
            #         position = lTreeInsertAt + offset;
            self.ltree.set_bit(position, True, self.lPath)
            #         ltree->setBit(position, true, &lPath);
            #     }
            # }
        pass

    def remove_edge(self, row: int, column: int):
        pass
        # void DKTree::removeEdge(unsigned long row, unsigned long column) {
        position_of_first: int = 0
        #     const unsigned long POSITION_OF_FIRST = 0;
        first_iteration: int = 1
        #     const unsigned long FIRST_ITERATION = 1;
        self.delete_this_edge(row, column, first_iteration, position_of_first)
        #     deleteThisEdge(row, column, FIRST_ITERATION, POSITION_OF_FIRST);
        # }

    def delete_this_edge(self, row: int, column: int, iteration: int, position_of_first: int) -> bool:
        pass
        # bool DKTree::deleteThisEdge(const unsigned long row, const unsigned long column,
        # const unsigned long iteration,
        #                             const unsigned long positionOfFirst) {
        offset: int = self.calculateOffset(row, column, iteration)
        #     unsigned long offset = calculateOffset(row, column, iteration);
        if position_of_first >= self.ttree.bits():
            #     if (positionOfFirst >= ttree->bits()) {
            return self.delete_l_tree_edge(position_of_first, offset)
            #         return deleteLTreeEdge(positionOfFirst, offset);
        elif self.ttree.access(position_of_first + offset, self.tPath):
            #     } else if (ttree->access(positionOfFirst + offset, &tPath)) {
            return self.delete_t_tree_edge(row, column, iteration, position_of_first, offset)
            #         return deleteTTreeEdge(row, column, iteration, positionOfFirst, offset);
        else:
            #     } else {
            # // the current ttree bit is already false, so no changes should be made, as it came here it parent should
            # // be true and therefore stay true
            return True
            #         return true;
            #     }
            # }
        pass

    def delete_t_tree_edge(self, row: int, column: int, iteration: int, position_of_first: int, offset: int) -> bool:
        pass
        # bool DKTree::deleteTTreeEdge(const unsigned long row, const unsigned long column,
        # const unsigned long iteration,
        #                              const unsigned long positionOfFirst, unsigned long offset) {
        #     // if the current position is true then check if after deleting the next edge any of
        #     its children are still true
        next_position_of_first: int = ((self.ttree.rank1(position_of_first + offset + 1, self.tPath))
                                       * Parameters.BLOCK_SIZE)
        #     unsigned long nextPositionOfFirst = (ttree->rank1(positionOfFirst + offset + 1, &tPath)) * BLOCK_SIZE;
        new_current_bit: bool = self.delete_this_edge(row, column, iteration + 1, next_position_of_first)
        #     bool newCurrentBit = deleteThisEdge(row, column, iteration + 1, nextPositionOfFirst);
        #     // if any of its children are still true this one will stay true and therefore so should its parent.
        if new_current_bit:
            #     if (newCurrentBit) {
            return True
            #         return true;
            #     }
        self.ttree.set_bit(position_of_first + offset, False, self.tPath)
        #     ttree->setBit(positionOfFirst + offset, false, &tPath);
        if iteration > 1:
            #     if (iteration > 1) {
            #         // if we aren't in the first iteration, see if any of the nodes in this block is still true
            only0s: bool = True
            #         bool only0s = true;
            i = 0
            while i < Parameters.BLOCK_SIZE and only0s:
                #         for (unsigned long i = 0; i < BLOCK_SIZE && only0s; i++) {
                if self.ttree.access(position_of_first + i, self.tPath):
                    #             if (ttree->access(positionOfFirst + i, &tPath)) {
                    only0s: bool = False
                    #                 only0s = false;
                    #             }
                    #         }
                i += 1
                #         // if all nodes are false this block can be deleted
            if only0s:
                #         if (only0s) {
                self.deleteBlockTtree(position_of_first)
                #             deleteBlockTtree(positionOfFirst);
                #         }
            return not only0s
            #         return !only0s; // if all bits in this block are false the parent should be false,
            #         else it should be true.
            #     }
        #     // we are in the first iteration, so no bits should be deleted
        return True
        #     return true;
        # }
        pass

    def delete_l_tree_edge(self, position_of_first: int, offset: int):
        pass
        # bool DKTree::deleteLTreeEdge(const unsigned long positionOfFirst, unsigned long offset) {
        #     // if the position is in the ltree, set the bit to false in the ltree
        l_tree_position_of_first: int = position_of_first - self.ttree.bits()
        #     unsigned long lTreePositionOfFirst = positionOfFirst - ttree->bits();
        l_tree_position: int = l_tree_position_of_first + offset
        #     unsigned long lTreePosition = lTreePositionOfFirst + offset;
        self.ltree.set_bit(l_tree_position, False, self.lPath)
        #     ltree->setBit(lTreePosition, false, &lPath);
        #     // check if there are any positive bits in this block
        only0s: bool = True
        #     bool only0s = true;
        i: int = 0
        while i < Parameters.BLOCK_SIZE and only0s:
            #     for (unsigned long i = 0; i < BLOCK_SIZE && only0s; i++) {
            if self.ltree.access(l_tree_position_of_first + i, self.lPath):
                #         if (ltree->access(lTreePositionOfFirst + i, &lPath)) {
                only0s = False
                #             only0s = false;
                #         }
            i += 1
            #     }
        #     // iff all bits in the block are 0 delete this block
        if only0s:
            #     if (only0s) {
            self.deleteBlockLtree(l_tree_position_of_first)
            #         deleteBlockLtree(lTreePositionOfFirst);
            #     }
        return not only0s
        #     return !only0s; // if this block is all false then its parent should be false,
        #     else the parent should be true
        # }
        pass

    def insert_entry(self) -> int:
        pass
        # unsigned long DKTree::insertEntry() {
        #     unsigned long insertedColumn;
        if len(self.freeColumns) != 0:
            #     if (!freeColumns.empty()) {
            #         // if a column in the middle was freed earlier then first use this column
            inserted_column: int = self.freeColumns[0]
            #         insertedColumn = freeColumns.front();
            del self.freeColumns[0]
            #         freeColumns.erase(freeColumns.begin());
        else:
            #     } else {
            if self.firstFreeColumn > self.matrixSize:
                #         if (firstFreeColumn > matrixSize) {
                self.increaseMatrixSize()
                #             increaseMatrixSize();
                #         }
            #         // if not use the last column
            inserted_column: int = self.firstFreeColumn
            #         insertedColumn = firstFreeColumn;
            self.firstFreeColumn += 1
            #         firstFreeColumn++;
            #     }
        #     // as unused places in the matrix have 0 everywhere no action on the bitvector is needed
        return inserted_column
        #     return insertedColumn;
        # }
        pass

    def delete_entry(self, a: int):
        pass
        # void DKTree::deleteEntry(unsigned long a) {
        self.checkArgument(a, "deleteEntry")
        #     checkArgument(a, "deleteEntry");
        allOthers: list[int] = [i for i in range(self.firstFreeColumn)]
        #     vector<unsigned long> allOthers;
        thisOne: list[int] = [a]
        #     vector<unsigned long> thisOne{a};
        #     for (unsigned long i = 0; i < firstFreeColumn; i++) {
        #         allOthers.push_back(i);
        #     }
        for i in range(len(self.freeColumns) - 1, -1, -1):
            #     for (long i = freeColumns.size() - 1; i >= 0; i--) {
            del allOthers[self.freeColumns[i]]
            #         allOthers.erase(allOthers.begin() + freeColumns[i]);
            #     }

        thisA: VectorData = VectorData().init_with_aentry(thisOne)
        #     VectorData thisA(thisOne);
        others: VectorData = VectorData().init_with_aentry(allOthers)
        #     VectorData others(allOthers);
        self.delete_edges(thisA, others)
        #     deleteEdges(thisA, others);
        self.delete_edges(others, thisA)
        #     deleteEdges(others, thisA);
        if a == self.firstFreeColumn - 1:
            #     if (a == firstFreeColumn - 1) {
            self.firstFreeColumn -= 1
            #         firstFreeColumn--;
        else:
            #     } else {
            self.freeColumns.append(a)
            #         freeColumns.push_back(a);
            self.freeColumns.sort()
            #         sort(freeColumns.begin(), freeColumns.end());
            #     }
        # }
        pass

    def delete_edges(self, rows: VectorData, columns: VectorData) -> bool:
        pass
        # bool DKTree::deleteEdges(VectorData &rows, VectorData &columns) {
        if rows.firstAt != columns.firstAt or rows.iteration != columns.iteration:
            #     if (rows.firstAt != columns.firstAt || rows.iteration != columns.iteration) {
            #         std::stringstream error;
            error = "deleteEdges: rows and columns asynch\n"
            #         error << "deleteEdges: rows and columns asynch\n";
            raise Exception(error)
            #         throw std::invalid_argument(error.str());
            #     }
        k: int = Parameters.k
        partitionSize: int = int(self.matrixSize / long_pow(k, rows.iteration))
        #     const unsigned long partitionSize = matrixSize / long_pow(k, rows.iteration);
        only0s: bool = True
        #     bool only0s = true;
        if partitionSize > 1:
            #     if (partitionSize > 1) { // we are looking at ttree stuff
            #         // sort the rows and columns according to which offsets they belong
            rowStart: list[int] = [-1 for _ in range(k)]
            #         int rowStart[k];
            rowEnd: list[int] = [-1 for _ in range(k)]
            #         int rowEnd[k];
            columnStart: list[int] = [-1 for _ in range(k)]
            #         int columnStart[k];
            columnEnd: list[int] = [-1 for _ in range(k)]
            #         int columnEnd[k];

            # for i in range(k):
            #     #         for (int i = 0; i < k; i++) {
            #     rowStart[i] = -1
            #     #             rowStart[i] = -1;
            #     rowEnd[i] = -1
            #     #             rowEnd[i] = -1;
            #     columnStart[i] = -1
            #     #             columnStart[i] = -1;
            #     columnEnd[i] = -1
            #     #             columnEnd[i] = -1;
            #     #         }
            self.splitEntriesOnOffset(rows, partitionSize, rowStart, rowEnd)
            #         splitEntriesOnOffset(rows, partitionSize, rowStart, rowEnd);
            self.splitEntriesOnOffset(columns, partitionSize, columnStart, columnEnd)
            #         splitEntriesOnOffset(columns, partitionSize, columnStart, columnEnd);
            #         // for each offset, there can be a relation if there is at least one row and one column and
            #         if its value is not 0.
            for offset in range(Parameters.BLOCK_SIZE - 1, -1, -1):
                #         for (int offset = BLOCK_SIZE - 1; offset >= 0; offset--) {
                rowOffset: int = int(offset / k)
                #             unsigned long rowOffset = offset / k;
                #             // std::cout << "rowOffset " << rowOffset << "\n";
                columnOffset: int = offset % k
                #             unsigned long columnOffset = offset % k;
                #             // std::cout << "columnOffset " << columnOffset << "\n";
                currentNode: int = rows.firstAt + offset
                #             unsigned long currentNode = rows.firstAt + offset;
                nodeSubtreeHasEdges: bool = self.ttree.access(currentNode, self.tPath)
                #             bool nodeSubtreeHasEdges = ttree->access(currentNode, &tPath);
                if nodeSubtreeHasEdges:
                    #             if (nodeSubtreeHasEdges) {
                    if not (rowStart[rowOffset] == -1 or columnStart[columnOffset] == -1):
                        #                 if (!(rowStart[rowOffset] == -1 || columnStart[columnOffset] == -1)) {
                        #                     // there can only be a relation if there is at least 1 element in both
                        #                     of them
                        #                     // rank function is exclusive so +1
                        nextNode: int = self.ttree.rank1(currentNode + 1, self.tPath) * Parameters.BLOCK_SIZE
                        #                     unsigned long nextNode = ttree->rank1(currentNode + 1, &tPath)
                        #                     * BLOCK_SIZE;
                        #                     // if there are edges in this subtree find the edges stored in the child
                        #                     nodes
                        nextIteration: int = rows.iteration + 1
                        #                     unsigned long nextIteration = rows.iteration + 1;
                        rowData: VectorData = VectorData().init_with_arguments(rows, rowStart[rowOffset], rowEnd[rowOffset],
                                                                   nextIteration, nextNode)
                        #                     VectorData rowData(rows, rowStart[rowOffset], rowEnd[rowOffset],
                        #                     nextIteration, nextNode);

                        columnData: VectorData = VectorData().init_with_arguments(columns, columnStart[columnOffset],
                                                                      columnEnd[columnOffset], nextIteration, nextNode)
                        #                     VectorData columnData(columns, columnStart[columnOffset],
                        #                     columnEnd[columnOffset], nextIteration, nextNode);
                        #                     // check if there are still edges left in its child nodes after deleting
                        #                     the edges from the rowData columnData
                        stillHasEdges: bool = self.delete_edges(rowData, columnData)
                        #                     bool stillHasEdges = deleteEdges(rowData, columnData);
                        #                     // if there still are edges in its child nodes then this node stays 1
                        #                     and so should its parent
                        if stillHasEdges:
                            #                     if (stillHasEdges) {
                            only0s = False
                            #                         only0s = false;
                        else:
                            #                     } else {
                            #                         // if there are no edges in its child nodes this edge can be set
                            #                         to 0
                            self.ttree.set_bit(currentNode, False, self.tPath)
                            #                         ttree->setBit(currentNode, false, &tPath);
                            #                     }
                    else:
                        #                 } else {
                        #                     // if this offset is 1 and there is no edge to delete, its parent also
                        #                     should know there are still edges
                        only0s = False
                        #                     only0s = false;
                        #                 }
                    #             }
                #         }
            if only0s and rows.iteration > 1:
                #         if (only0s && rows.iteration > 1) {
                #             // if there are no more edges in this block it can be deleted
                self.deleteBlockTtree(rows.firstAt)
                #             deleteBlockTtree(rows.firstAt);
                #         }
            else:
                #     } else { // we look at ltree stuff
                only0s = self.deleteEdgesFromLTree(rows, columns)
                #         only0s = deleteEdgesFromLTree(rows, columns);
                #     }
        return not only0s
        #     return !only0s;
        # }
        pass

    def deleteEdgesFromLTree(self, rows: VectorData, columns: VectorData):
        pass
        # bool DKTree::deleteEdgesFromLTree(VectorData &rows, VectorData &columns) {
        only0s: bool = True
        #     bool only0s = true;
        partitionSize: int = int(self.matrixSize / long_pow(Parameters.k, rows.iteration))
        #     const unsigned long partitionSize = matrixSize / long_pow(k, rows.iteration);
        if partitionSize > 1:
            #     if (partitionSize > 1) {
            #         std::stringstream error;
            #         error << "findEdgesInLTree: not lTree iteration\n";
            raise Exception("findEdgesInLTree: not lTree iteration\n")
            #         throw std::invalid_argument(error.str());
            #     }
        ltreeposition: int = rows.firstAt - self.ttree.bits()
        #     unsigned long ltreeposition = rows.firstAt - ttree->bits();
        for i in range(rows.start, rows.end):
            #     for (unsigned long i = rows.start; i < rows.end; i++) {
            for j in range(columns.start, columns.end):
                #         for (unsigned long j = columns.start; j < columns.end; j++) {
                offset: int = self.calculateOffset(rows.entry[i], columns.entry[j], rows.iteration)
                #             unsigned long offset = calculateOffset(rows.entry[i], columns.entry[j], rows.iteration);
                nodePosition: int = ltreeposition + offset
                #             unsigned long nodePosition = ltreeposition + offset;
                self.ltree.set_bit(nodePosition, False, self.lPath)
                #             ltree->setBit(nodePosition, false, &lPath);
                #         }
            pass
            #     }
        for offset in range(Parameters.BLOCK_SIZE):
            #     for (unsigned long offset = 0; offset < BLOCK_SIZE && only0s; offset++) {
            if self.ltree.access(ltreeposition + offset, self.lPath):
                #         if (ltree->access(ltreeposition + offset, &lPath)) {
                only0s = False
                break
                #             only0s = false;
                #         }
        pass
        #     }
        if only0s:
            #     if (only0s) {
            self.deleteBlockLtree(ltreeposition)
            #         deleteBlockLtree(ltreeposition);
            #     }
        return only0s
        #     return only0s;
        # }
        pass

    def report_edge(self, a: int, b: int):
        pass
        # bool DKTree::reportEdge(unsigned long a, unsigned long b) {
        #     // tests if both positions exist
        self.checkArgument(a, "reportEdge")
        #     checkArgument(a, "reportEdge");
        self.checkArgument(b, "reportEdge")
        #     checkArgument(b, "reportEdge");
        iteration: int = 1
        #     unsigned long iteration = 1;
        tmax: int = self.ttree.bits()
        #     unsigned long tmax = ttree->bits();
        #     // the first position is 0+offset of the first iteration
        position: int = self.calculateOffset(a, b, iteration)
        #     unsigned long position = calculateOffset(a, b, iteration);;
        centry: bool = True
        #     bool centry = true; // to get the loop started
        iteration, position, centry = self.traverseToFirst0OrEndOfTTree(a, b, iteration, position, centry)
        #     traverseToFirst0OrEndOfTTree(a, b, iteration, position, centry);
        #     // if the last bit found in the ttree is a 1 then find the final result in the ltree
        if centry:
            #     if (centry) {
            ltreePosition: int = position - tmax
            #         unsigned long ltreePosition = position - tmax;
            centry = self.ltree.access(ltreePosition, self.lPath)
            #         centry = ltree->access(ltreePosition, &lPath);
            #     }
        return centry
        #     return centry;
        # }
        pass

    def report_all_edges(self, a: list[int], b: list[int]):
        pass
        # vector<std::pair<unsigned long, unsigned long>> DKTree::reportAllEdges(const vector<unsigned long> &A,
        #                                                                        const vector<unsigned long> &B) {
        rows_a: list[int] = a
        #     vector<unsigned long> rowsA(A);
        columns_b: list[int] = b
        #     vector<unsigned long> columnsB(B);
        self.sortAndCheckVector(rows_a)
        #     sortAndCheckVector(rowsA);
        self.sortAndCheckVector(columns_b)
        #     sortAndCheckVector(columnsB);
        rows: VectorData = VectorData().init_with_aentry(rows_a)
        #     VectorData rows(rowsA);
        columns: VectorData = VectorData().init_with_aentry(columns_b)
        #     VectorData columns(columnsB);
        findings: list[tuple[int, int]] = []
        #     vector<std::pair<unsigned long, unsigned long>> findings;
        self.findAllEdges(rows, columns, findings)
        #     findAllEdges(rows, columns, findings);
        #     return findings;
        # }
        return findings
        pass

    def findAllEdges(self, rows: VectorData, columns: VectorData, findings: list[tuple[int, int]]):
        pass
        # void
        # DKTree::findAllEdges(VectorData &rows, VectorData &columns,
        #                      vector<std::pair<unsigned long, unsigned long>> &findings) {
        if rows.firstAt != columns.firstAt or rows.iteration != columns.iteration:
            #     if (rows.firstAt != columns.firstAt || rows.iteration != columns.iteration) {
            #         std::stringstream error;
            #         error << "findAllEdges: rows and columns asynch\n";
            raise Exception("findAllEdges: rows and columns asynch\n")
            #         throw std::invalid_argument(error.str());
            #     }
        k: int = Parameters.k
        partitionSize: int = int(self.matrixSize / long_pow(k, rows.iteration))
        #     const unsigned long partitionSize = matrixSize / long_pow(k, rows.iteration);
        if partitionSize > 1:
            #     if (partitionSize > 1) { // we are looking at ttree stuff
            #         // sort the rows and columns according to which offsets they belong
            rowStart: list[int] = [-1 for _ in range(k)]
            #         int rowStart[k];
            rowEnd: list[int] = [-1 for _ in range(k)]
            #         int rowEnd[k];
            columnStart: list[int] = [-1 for _ in range(k)]
            #         int columnStart[k];
            columnEnd: list[int] = [-1 for _ in range(k)]
            #         int columnEnd[k];
            #         for (int i = 0; i < k; i++) {
            #             rowStart[i] = -1;
            #             rowEnd[i] = -1;
            #             columnStart[i] = -1;
            #             columnEnd[i] = -1;
            #         }
            self.splitEntriesOnOffset(rows, partitionSize, rowStart, rowEnd)
            #         splitEntriesOnOffset(rows, partitionSize, rowStart, rowEnd);
            self.splitEntriesOnOffset(columns, partitionSize, columnStart, columnEnd)
            #         splitEntriesOnOffset(columns, partitionSize, columnStart, columnEnd);
            #         // for each offset, there can be a relation if there is at least one row and one column and
            #         if its value is not 0.
            for offset in range(Parameters.BLOCK_SIZE):
                #         for (unsigned long offset = 0; offset < BLOCK_SIZE; offset++) {
                rowOffset: int = int(offset / k)
                #             unsigned long rowOffset = offset / k;
                columnOffset: int = offset % k
                #             unsigned long columnOffset = offset % k;
                if not (rowStart[rowOffset] == -1 or columnStart[columnOffset] == -1):
                    #             if (!(rowStart[rowOffset] == -1 || columnStart[columnOffset] == -1)) {
                    #                 // there can only be a relation if there is at least 1 element in both of them
                    currentNode: int = rows.firstAt + offset
                    #                 unsigned long currentNode = rows.firstAt + offset;
                    nodeSubtreeHasEdges = self.ttree.access(currentNode, self.tPath)
                    #                 bool nodeSubtreeHasEdges = ttree->access(currentNode, &tPath);
                    if nodeSubtreeHasEdges:
                        #                 if (nodeSubtreeHasEdges) {
                        #                     // rank function is exclusive so +1
                        nextNode: int = self.ttree.rank1(currentNode + 1, self.tPath) * Parameters.BLOCK_SIZE
                        #                     unsigned long nextNode = ttree->rank1(currentNode + 1, &tPath)
                        #                     * BLOCK_SIZE;
                        #                     // if there are edges in this subtree find the edges stored in the child
                        #                     nodes
                        nextIteration:int = rows.iteration + 1
                        #                     unsigned long nextIteration = rows.iteration + 1;
                        rowData: VectorData = VectorData().init_with_arguments(rows, rowStart[rowOffset], rowEnd[rowOffset],
                                                                   nextIteration, nextNode)
                        #                     VectorData rowData(rows, rowStart[rowOffset], rowEnd[rowOffset],
                        #                     nextIteration, nextNode);
                        columnData: VectorData = VectorData().init_with_arguments(columns, columnStart[columnOffset],
                                                                      columnEnd[columnOffset], nextIteration, nextNode)
                        #                     VectorData columnData(columns, columnStart[columnOffset],
                        #                     columnEnd[columnOffset],
                        #                     nextIteration,
                        #                                           nextNode);
                        self.findAllEdges(rowData, columnData, findings)
                        #                     findAllEdges(rowData, columnData, findings);
                    #                 }
                #             }
            #         }
        else:
            #     } else { // we look at ltree stuff
            self.findEdgesInLTree(rows, columns, findings)
            #         findEdgesInLTree(rows, columns, findings);
            #     }
        # }
        pass

    def findEdgesInLTree(self, rows: VectorData, columns: VectorData, findings: list[tuple[int, int]]):
        pass
        # void DKTree::findEdgesInLTree(const VectorData &rows, const VectorData &columns,
        #                               vector<pair<unsigned long, unsigned long>> &findings) {
        k: int = Parameters.k
        partitionSize: int = int(self.matrixSize / long_pow(k, rows.iteration))
        #     const unsigned long partitionSize = matrixSize / long_pow(k, rows.iteration);
        if partitionSize > 1:
            #     if (partitionSize > 1) {
            #         std::stringstream error;
            #         error << "findEdgesInLTree: not lTree iteration\n";
            raise Exception("findEdgesInLTree: not lTree iteration\n")
            #         throw std::invalid_argument(error.str());
            #     }
        ltreeposition: int = rows.firstAt - self.ttree.bits()
        #     unsigned long ltreeposition = rows.firstAt - ttree->bits();
        for i in range(rows.start, rows.end):
            #     for (unsigned long i = rows.start; i < rows.end; i++) {
            for j in range(columns.start, columns.end):
                #         for (unsigned long j = columns.start; j < columns.end; j++) {
                offset: int = self.calculateOffset(rows.entry[i], columns.entry[j], rows.iteration)
                #             unsigned long offset = calculateOffset(rows.entry[i], columns.entry[j], rows.iteration);
                nodePosition: int = ltreeposition + offset
                #             unsigned long nodePosition = ltreeposition + offset;
                hasEdge: bool = self.ltree.access(nodePosition, self.lPath)
                #             bool hasEdge = ltree->access(nodePosition, &lPath);
                if hasEdge:
                    #             if (hasEdge) {
                    edge: tuple[int, int] = (rows.entry[i], columns.entry[j])
                    #                 pair<unsigned long, unsigned long> edge(rows.entry[i], columns.entry[j]);
                    findings.append(edge)
                    #                 findings.push_back(edge);
                    #             }
                #         }
            #     }
        # }
        pass

    def splitEntriesOnOffset(self, entries: VectorData, partitionSize: int, entryStart: list[int], entryEnd: list[int]):
        pass
        # void DKTree::splitEntriesOnOffset(const VectorData &entries, const unsigned long partitionSize,
        # int *entryStart, int *entryEnd) const {
        formerPartitionSize: int = partitionSize * Parameters.k
        #     const unsigned long formerPartitionSize = partitionSize * k;
        offsetStarted: int = 0
        #     unsigned long offsetStarted;
        for i in range(entries.start, entries.end):
            #     for (unsigned long i = entries.start; i < entries.end; i++) {
            entryInBlock: int = entries.entry[i] % formerPartitionSize
            #         unsigned long entryInBlock = entries.entry[i] % formerPartitionSize;
            entryOffset: int = int(entryInBlock / partitionSize)
            #         unsigned long entryOffset = (entryInBlock / partitionSize);
            if entryStart[entryOffset] == -1:
                #         if (entryStart[entryOffset] == -1) {
                entryStart[entryOffset] = i
                #             entryStart[entryOffset] = i;
                if i > entries.start:
                    #             if (i > entries.start) {
                    entryEnd[offsetStarted] = i
                    #                 entryEnd[offsetStarted] = i;
                    #             }
                offsetStarted = entryOffset
                #             offsetStarted = entryOffset;
                #         }
            if i == entries.end - 1:
                #         if (i == entries.end - 1) {
                entryEnd[offsetStarted] = i + 1
                #             entryEnd[offsetStarted] = i + 1;
                #         }
            #     }
        # }
        pass

    def printtt(self):
        pass
        # void DKTree::printtt() {
        print("ttree:")
        #     cout << "ttree:" << endl;
        self.printttree(self.ttree, 1)
        #     printttree(ttree);
        print()
        #     printf("\n");
        print("ltree:")
        #     cout << "ltree:" << endl;
        self.printltree(self.ltree, 1)
        #     printltree(ltree);
        print()
        #     printf("\n");
        # }
        pass

    def printttree(self, tree: TTree, depth: int):
        pass
        # void DKTree::printttree(TTree *tree, unsigned long depth) {
        prefix: str = ''
        #     std::string prefix;
        for i in range(depth):
            #     for (unsigned long i = 0; i < depth; i++) {
            prefix += "| "
            #         prefix += "| ";
            #     }
        if tree.isLeaf:
            #     if (tree->isLeaf) {
            bv = tree.node.leafNode.bv
            #         auto &bv = tree->node.leafNode->bv;
            print(prefix)
            #         printf("%s", prefix.c_str());
            for b in bv.data:
                #         for (auto b : bv.data) {
                print(b)
                #             printf("%i", (bool) b);
                #         }
            print()
            #         printf("\n");
        else:
            #     } else {
            print(prefix, tree.node.internalNode.bits(), tree.node.internalNode.ones())
            #         printf("%s (%lu bits, %lu ones)\n", prefix.c_str(), tree->node.internalNode->bits(),
            #         tree->node.internalNode->ones());
            for entry in tree.node.internalNode.entries:
                #         for (auto &entry : tree->node.internalNode->entries) {
                if entry.p is None:
                    #             if (entry.P == nullptr) {
                    break
                    #                 break;
                    #             }
                self.printttree(entry.p, depth + 1)
                #             printttree(entry.P, depth + 1);
                #         }
            #     }
        # }
        pass

    def printltree(self, tree: LTree, depth: int):
        pass
        # void DKTree::printltree(LTree *tree, unsigned long depth) {
        prefix: str = ''
        #     std::string prefix;
        for i in range(depth):
            #     for (unsigned long i = 0; i < depth; i++) {
            prefix += "| "
            #         prefix += "| ";
            #     }
        if tree.isLeaf:
            #     if (tree->isLeaf) {
            bv = tree.node.leafNode.bv
            #         auto &bv = tree->node.leafNode->bv;
            print(prefix)
            #         printf("%s", prefix.c_str());
            for b in bv.data:
                #         for (auto b : bv.data) {
                print(b)
                #             printf("%i", (bool) b);
                #         }
            print()
            #         printf("\n");
        else:
            #     } else {
            print(prefix, tree.node.internalNode.bits())
            #         printf("%s (%lu bits)\n", prefix.c_str(), tree->node.internalNode->bits());
            for entry in tree.node.internalNode.entries:
                #         for (auto &entry : tree->node.internalNode->entries) {
                if entry.p is None:
                    #             if (entry.P == nullptr) {
                    break
                    #                 break;
                    #             }
                self.printltree(entry.p, depth + 1)
                #             printltree(entry.P, depth + 1);
                #         }
            #     }
        # }
        pass

    def increaseMatrixSize(self):
        pass
        # void DKTree::increaseMatrixSize() {
        FIRST_BIT: int = 0
        #     const unsigned long FIRST_BIT = 0;
        #     // if the matrix is full, increase the size by multiplying with k
        self.matrixSize *= Parameters.k
        #     matrixSize *= k;
        #     // position +1 since paper has rank including the position, but function is exclusive position
        if self.ttree.rank1(Parameters.BLOCK_SIZE, self.tPath) > 0:
            #     if (ttree->rank1(BLOCK_SIZE, &tPath) > 0) {
            #         // if there already is a 1 somewhere in the matrix, add a new block
            #         // in front of the bitvector and set the first bit to 1
            self.insertBlockTtree(FIRST_BIT)
            #         insertBlockTtree(FIRST_BIT);
            self.ttree.set_bit(FIRST_BIT, True, self.tPath)
            #         ttree->setBit(FIRST_BIT, true, &tPath);
            #     }
        # }
        pass

    def calculateOffset(self, row: int, column: int, iteration: int):
        pass
        # unsigned long
        # DKTree::calculateOffset(const unsigned long row, const unsigned long column, const unsigned long iteration) {
        #     // first remove the rows and columns not belonging to the current block
        formerPartitionSize: int = int(self.matrixSize / long_pow(Parameters.k, iteration - 1))
        #     unsigned long formerPartitionSize = matrixSize / long_pow(k, iteration - 1);
        rowInBlock: int = row % formerPartitionSize
        #     unsigned long rowInBlock = row % formerPartitionSize;
        columnInBlock: int = column % formerPartitionSize
        #     unsigned long columnInBlock = column % formerPartitionSize;
        #     //calculate the offset, each row partition adds k to the offset, each column partition 1
        partitionSize: int = int(self.matrixSize / long_pow(Parameters.k, iteration))
        #     unsigned long partitionSize = matrixSize / long_pow(k, iteration);
        if partitionSize == 0:
            #     if (partitionSize == 0) {
            raise Exception("partition size is 0\n")
            #         throw std::invalid_argument("partition size is 0\n");
            #     }
        xxx = int(rowInBlock / partitionSize)
        rowOffset: int = int(Parameters.k * xxx)
        #     unsigned long rowOffset = k * (rowInBlock / partitionSize);
        columnOffset: int = int(columnInBlock / partitionSize)
        #     unsigned long columnOffset = columnInBlock / partitionSize;
        return rowOffset + columnOffset
        #     return rowOffset + columnOffset;
        # }
        pass

    def checkArgument(self, a: int, functionName: str):
        pass
        # void DKTree::checkArgument(unsigned long a, std::string functionName) {
        if a >= self.firstFreeColumn:
            #     if (a >= firstFreeColumn) {
            #         std::stringstream error;
            error = (f"{functionName}: invalid argument {a}, position not occupied in matrix, "
                     + f"firstfreecolumn = {self.firstFreeColumn}\n")
            #         error << functionName << ": invalid argument " << a << ", position not occupied in matrix,
            #         firstfreecolumn = "<< firstFreeColumn <<"\n";
            raise Exception(error)
            #         throw std::invalid_argument(error.str());
        elif a < 0:
            pass
            #     } else if (a < 0) {
            error = f"{functionName}: invalid argument {a}, position does not exist\n"
            #         std::stringstream error;
            #         error << functionName << ": invalid argument " << a << ", position does not exist\n";
            raise Exception(error)
            #         throw std::invalid_argument(error.str());
            #     } else {
        else:
            pass
            for fc in self.freeColumns:
                pass
                #         for (auto &fc: freeColumns) {
                if fc == a:
                    pass
                    #             if (fc == a) {
                    error = f"{functionName}: invalid argument {a}, position was deleted from matrix\n"
                    #                 std::stringstream error;
                    #                 error << functionName << ": invalid argument " << a <<
                    #                 ", position was deleted from matrix\n";
                    raise Exception(error)
                    #                 throw std::invalid_argument(error.str());
                    #             }
                    #         }
                    #     }
        # }
        pass

    def traverseToFirst0OrEndOfTTree(self, row, column, iteration, position, cEntry) -> tuple[int, int, bool]:
        pass
        # void DKTree::traverseToFirst0OrEndOfTTree(unsigned long row, unsigned long column, unsigned long &iteration,
        #                                           unsigned long &position, bool &cEntry) {
        tmax: int = self.ttree.bits()
        #     unsigned long tmax = ttree->bits();
        while cEntry and position < tmax:
            #     // while the current position is a 1 and the end of the ttree is not reached, access the next bit
            #     while (cEntry && position < tmax) {
            cEntry = self.ttree.access(position, self.tPath)
            #         cEntry = ttree->access(position, &tPath);
            if cEntry:
                #         if (cEntry) {
                iteration += 1
                #             iteration++;
                offset: int = self.calculateOffset(row, column, iteration)
                #             unsigned long offset = calculateOffset(row, column, iteration);
                # // position +1 since paper has rank including the position, but function is exclusive position
                positionOfFirst: int = self.ttree.rank1(position + 1, self.tPath) * Parameters.BLOCK_SIZE
                #             unsigned long positionOfFirst = ttree->rank1(position + 1, &tPath) * BLOCK_SIZE;
                position = positionOfFirst + offset
                #             position = positionOfFirst + offset;
                #         }
            #     }
        # }
        return iteration, position, cEntry
        pass

    def sortAndCheckVector(self, elements: list[int]):
        pass
        # void DKTree::sortAndCheckVector(vector<unsigned long> &elements) {
        if len(elements) == 0:
            #     if (elements.empty()) {
            #         std::stringstream error;
            error: str = "sortAndCheckVector: invalid argument, empty input vector \n"
            #         error << "sortAndCheckVector: invalid argument, empty input vector \n";
            raise Exception(error)
            #         throw std::invalid_argument(error.str());
            #     }
        #     // sort and delete doubles
        elements.sort()
        #     sort(elements.begin(), elements.end());
        elements = list(set(elements))
        #     elements.erase(unique(elements.begin(), elements.end()), elements.end());
        functionName: str = "reportAllEdges"
        #     std::string functionName = "reportAllEdges";
        for element in elements:
            #     for (auto element:elements) {
            self.checkArgument(element, functionName)
            #         checkArgument(element, functionName);
            #     }
        # }
        pass

    def insertBlockTtree(self, position: int):
        pass
        # void DKTree::insertBlockTtree(unsigned long position) {
        newRoot: TTree = self.ttree.insertBlock(position, self.tPath)
        #     TTree *newRoot = ttree->insertBlock(position, &tPath);
        if newRoot is not None:
            #     if (newRoot != nullptr) {
            self.ttree = newRoot
            #         ttree = newRoot;
            #     }
        self.tPath = []
        #     tPath.clear();
        # }
        pass

    def insertBlockLtree(self, position: int):
        pass
        # void DKTree::insertBlockLtree(unsigned long position) {
        newRoot: LTree = self.ltree.insertBlock(position, self.lPath)
        #     LTree *newRoot = ltree->insertBlock(position, &lPath);
        if newRoot is not None:
            #     if (newRoot != nullptr) {
            self.ltree = newRoot
            #         ltree = newRoot;
            #     }
        self.lPath = []
        #     lPath.clear();
        # }
        pass

    def deleteBlockTtree(self, position: int):
        pass
        # void DKTree::deleteBlockTtree(unsigned long position) {
        newRoot: TTree = self.ttree.deleteBlock(position, self.tPath)
        #     TTree *newRoot = ttree->deleteBlock(position, &tPath);
        if newRoot is not None:
            #     if (newRoot != nullptr) {
            self.ttree = newRoot
            #         ttree = newRoot;
            #     }
            self.tPath = []
        #     tPath.clear();
        # }
        pass

    def deleteBlockLtree(self, position: int):
        pass
        # void DKTree::deleteBlockLtree(unsigned long position) {
        newRoot: LTree = self.ltree.deleteBlock(position, self.lPath)
        #     LTree *newRoot = ltree->deleteBlock(position, &lPath);
        if newRoot is not None:
            #     if (newRoot != nullptr) {
            self.ltree = newRoot
            #         ltree = newRoot;
            #     }
        #     lPath.clear();
        # }
        pass

    def memory_usage(self) -> int:
        pass
        # unsigned long DKTree::memoryUsage() {
        base: int = sys.getsizeof(DKTree)
        #     unsigned long base = sizeof(DKTree),
        tSize = self.ttree.memoryUsage()
        #         tSize = ttree->memoryUsage(),
        tBits = self.ttree.bits()
        #         tBits = ttree->bits(),
        lSize = self.ltree.memoryUsage()
        #         lSize = ltree->memoryUsage(),
        lBits = self.ltree.bits()
        #         lBits = ltree->bits(),
        tPathSize = len(self.tPath) * sys.getsizeof(Nesbo)
        #         tPathSize = tPath.size() * sizeof(Nesbo),
        lPathSize = len(self.lPath) * sys.getsizeof(LNesbo)
        #         lPathSize = lPath.size() * sizeof(LNesbo),
        freeSize = len(self.freeColumns) * sys.getsizeof(int)
        #         freeSize = freeColumns.size() * sizeof(unsigned long);
        print(f"  TTree: {tSize}\n")
        #     printf("  TTree: %lu\n", tSize);
        print(f"    Of which bitvector: {(tBits + 7) / 8}\n")
        #     printf("    Of which bitvector: %lu\n", (tBits + 7) / 8);
        print(f"  LTree: {lSize}\n")
        #     printf("  LTree: %lu\n", lSize);
        print(f"    Of which bitvector: {(lBits + 7) / 8}\n")
        #     printf("    Of which bitvector: %lu\n", (lBits + 7) / 8);
        print(f"  TPath: {tPathSize}\n")
        #     printf("  TPath: %lu\n", tPathSize);
        print(f"  LPath: {lPathSize}\n")
        #     printf("  LPath: %lu\n", lPathSize);
        print(f"  Free Columns: {freeSize}\n")
        #     printf("  Free Columns: %lu\n", freeSize);
        return base + tSize + lSize + tPathSize + lPathSize + freeSize
        #     return base + tSize + lSize + tPathSize + lPathSize + freeSize;
        # }
        pass

    def insert_entries(self, count: int):  # 2024/8/6
        for i in range(count):
            self.insert_entry()
