#!/usr/bin/python3


from PyQt5.QtCore import QLineF, QPointF

import math
import time
import random
import numpy as np

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        # call either the banded or non banded alignment function
        if not self.banded:
            score, alignment1, alignment2 = self.compute_alignment(seq1, seq2)
        else:
            score, alignment1, alignment2 = self.compute_alignment_banded(seq1, seq2)

        alignment1 = alignment1.format(len(seq1), align_length, ',BANDED' if banded else '')
        alignment2 = alignment2.format(len(seq2), align_length, ',BANDED' if banded else '')

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}


    # O(kn) where k is the bandwidth (7) and n is the number of rows.
    def compute_alignment_banded(self, upper_seq, left_seq):
        # only care about the first characters up to the max of the sequence
        # O(n)
        upper_seq = upper_seq[:self.MaxCharactersToAlign]
        left_seq = left_seq[:self.MaxCharactersToAlign]

        # If the lengths of the two sequences are more than 3 bp longer in length, we fail with banded algorithm.
        if abs(len(upper_seq) - len(left_seq)) > (MAXINDELS * 2) + 1:  # O(1) constant look up)
            return math.inf, "No Alignment Possible", "No Alignment Possible"

        # Set ROW and COLUMN length
        # column length to the bandwidth
        columns = (MAXINDELS * 2) + 1  # O(1)

        # Set Row Length:
        # for the way the 7 column table is set up, the longer sequence needs to be on the left.
        # O(1) constant time operations, lookups
        if len(left_seq) < len(upper_seq):
            rows = len(upper_seq) + 1  # update the row length with the longer sequence
            upper_seq, left_seq = left_seq, upper_seq  # swap the two sequences
        else:
            rows = len(left_seq) + 1  # row length setting

        # initialize table to infinity
        # O(k*n) time, O(k*n) space
        dp_table = [[(math.inf, "null") for x in range(columns)] for y in range(rows)]

        # Start point is actually [0][3] in a banded 2D array.
        dp_table[0][3] = (0, "start_pt")

        # initialize the edges in a banded array, the top edge will be different.
        # O(1) constant time base case initialization
        for i in range(1, MAXINDELS + 1):
            val = (INDEL * i)
            dp_table[0][i + MAXINDELS] = (val, "left")
        for i in range(1, MAXINDELS + 1):
            val = (INDEL * i)
            dp_table[i][MAXINDELS - i] = (val, "top")

        # Total time and space complexity is O(kn)
        # O(k*n) where n is the number of rows in the table, k is a constant (7)
        for i in range(0, rows):
            for j in range(0, columns):

                # check if we are allowed to view this cell
                if dp_table[i][j][1] == "null":
                    nw_value = math.inf

                    # These cases handle the banded properties of our 2D array and correctly look up where the top,
                    # left, and diagonal are.
                    if j == 6:
                        top = math.inf  # if we have reached the end of the band, there is nothing on top.
                        left = dp_table[i][j - 1][0]  # the left value will be directly left (normal)
                        diagonal = dp_table[i - 1][j][0]  # the diagonal value will be one row up directly above it.
                    elif j == 0:
                        top = dp_table[i - 1][j + 1][0]  # we go one column left and and one row up to find the the top
                        left = math.inf  # we can't go left if its the start
                        diagonal = dp_table[i - 1][j][0]  # the diagonal value will be one row up directly above it.
                    else:
                        top = dp_table[i - 1][j + 1][0]  # we go one column left and and one row up to find the the top
                        left = dp_table[i][j - 1][0]  # the left value will be directly left (normal)
                        diagonal = dp_table[i - 1][j][0]  # the diagonal value will be one row up directly above it.

                    # if diagonal is an option
                    if diagonal != math.inf:
                        # if we aren't out of bounds of our banded 2D array
                        if len(upper_seq) < i - (MAXINDELS - j):
                            dp_table[i][j] = (math.inf, dp_table[i][j][1])
                            continue  # if the value is "out of bounds," restart the for loop

                        # Check if we have a match or a substitution
                        else:
                            # get the correct indices for our letters
                            if left_seq[i - 1] == upper_seq[i - (4 - j)]:
                                nw_value = diagonal + MATCH
                                dp_table[i][j] = (dp_table[i][j][0], "diagonal")
                            else:
                                nw_value = diagonal + SUB
                                dp_table[i][j] = (nw_value, "diagonal")

                    # if top
                    if top != math.inf:
                        if (top + INDEL) <= nw_value:
                            nw_value = top + INDEL
                            dp_table[i][j] = (nw_value, "top")

                    # if left
                    if left != math.inf:
                        if (left + INDEL) <= nw_value:
                            nw_value = left + INDEL
                            dp_table[i][j] = (nw_value, "left")

                    # Update with new needleman wunsch score
                    dp_table[i][j] = (nw_value, dp_table[i][j][1])

        # Find which column we are in (final column location)
        # O(1) constant time (max MAXINDELS iterations)
        final_col_location = MAXINDELS
        while dp_table[rows - 1][final_col_location][0] == math.inf:
            final_col_location -= 1

        # Get the final alignments from the sequences O(kn)
        align_row, align_col = self.alignment_banded(upper_seq, left_seq, dp_table, rows, final_col_location)

        # Get final needleman wunsch score from the table.
        nw_value = dp_table[rows - 1][final_col_location][0]

        # Return the final variables O(1)
        return nw_value, align_row, align_col

    # Time complexity is (k*n) and space complexity is O(k*n) (same as above)
    # Where k is the constant bandwidth 7 and n is the number of rows.
    def alignment_banded(self, upper_seq, left_seq, dp_table, num_rows, col_val):
        align_row = ""
        align_col = ""

        # Select the starting pointer
        curr_val_ptr = dp_table[num_rows - 1][col_val][1]

        # iterates until it finds the start point, at worst is O(n+ m)
        while curr_val_ptr != "start_pt":

            if curr_val_ptr == "diagonal":
                align_row = upper_seq[num_rows - 1 - (4 - col_val)] + align_row
                align_col = left_seq[num_rows - 2] + align_col
                num_rows -= 1  # we moved a single row

            elif curr_val_ptr == "left":
                align_row = upper_seq[num_rows - 1 - (4 - col_val)] + align_row
                align_col = "-" + align_col
                col_val -= 1  # we moved a column

            elif curr_val_ptr == "top":
                align_row = "-" + align_row
                align_col = left_seq[num_rows - 2] + align_col
                # we need to move both a column and a row
                col_val += 1
                num_rows -= 1

            curr_val_ptr = dp_table[num_rows - 1][col_val][1]

        if align_row == "":
            align_row = "No Alignment Possible"
        else:
            align_row = align_row[:100]
        if align_col == "":
            align_col = "No Alignment Possible"
        else:
            align_col = align_col[:100]

        return align_row, align_col

    # Unrestricted algorithm for calculating alignment and score
    # Time and space complexity are both O(n*m)
    # Where m is num rows and n is num_cols
    def compute_alignment(self, left_seq, top_seq):
        # Get the dimensions of the table
        num_cols = min(len(top_seq), self.MaxCharactersToAlign) + 1
        num_rows = min(len(left_seq), self.MaxCharactersToAlign) + 1

        # Only consider up to the num of characters
        top_seq = top_seq[:self.MaxCharactersToAlign]
        left_seq = left_seq[:self.MaxCharactersToAlign]

        # Initialize the 2D array of tuples.
        # O(n + m)
        dp_table = [[(0, "null") for i in range(num_cols)] for j in range(num_rows)]

        # Set the base cases for the top row and left row, and the starting point
        for i in range(1, num_rows):
            val = i * INDEL
            dp_table[i][0] = (val, "top")
        for j in range(1, num_cols):
            val = j * INDEL
            dp_table[0][j] = (val, "left")
        dp_table[0][0] = (0, "start_pt")

        # Iterate through every value in our 2D array row by row.
        # O(n * m) time complexity
        for i in range(1, num_rows):
            for j in range(1, num_cols):
                # grab the current letters from the sequence
                top_let = top_seq[j - 1]
                left_let = left_seq[i - 1]

                # For the current value in the array, grab the correct surrounding vals in the dp_table
                top = dp_table[i - 1][j]
                diagonal = dp_table[i - 1][j - 1]
                left = dp_table[i][j - 1]

                # Diagonal case is either a substitution or a match
                if left_let == top_let:
                    nw_val = diagonal[0] + MATCH
                    dp_table[i][j] = (nw_val, "diagonal")
                else:
                    nw_val = diagonal[0] + SUB
                    dp_table[i][j] = (nw_val, "diagonal")

                # Top case criterion
                if top[0] + INDEL <= nw_val:
                    nw_val = top[0] + INDEL
                    dp_table[i][j] = (nw_val, "top")

                # Left case criterion
                if left[0] + INDEL <= nw_val:
                    nw_val = left[0] + INDEL
                    dp_table[i][j] = (nw_val, "left")

        # Run traceback with the alignment function to figure out the optimal path
        align1, align2 = self.find_alignment(top_seq, left_seq, dp_table)

        # Grab the needleman wunsch optimal value from the lower right corner
        nw_val = dp_table[num_rows - 1][num_cols - 1][0]

        return nw_val, align1, align2

    # O(nm) time
    # O(n + m) space to hold the two strings
    def find_alignment(self, top_sequence, left_sequence, dp_table):
        num_rows = len(left_sequence)
        num_cols = len(top_sequence)
        align_col = ""
        align_row = ""

        # While we have not reached the starting point:
        current_point = dp_table[num_rows][num_cols][1]
        while current_point != "start_pt":

            # Criteria for adding chars to our alignment sequences if we are diagonal
            if current_point == "diagonal":
                align_col = top_sequence[num_cols - 1] + align_col
                align_row = left_sequence[num_rows - 1] + align_row
                num_cols -= 1
                num_rows -= 1

            # Adding alignment chars if our previous value was located up.
            elif current_point == "top":
                align_col = "-" + align_col
                align_row = left_sequence[num_rows - 1] + align_row
                num_rows -= 1

            # Criteria for adding alignment chars if the previous value was left of the current value.
            elif current_point == "left":
                align_col = top_sequence[num_cols - 1] + align_col
                align_row = "-" + align_row
                num_cols -= 1

            # Grab the new current point and restart the loop
            current_point = dp_table[num_rows][num_cols][1]

        # If we couldn't align anything, then we return no alignment possible.  Otherwise we truncate and return only
        # 100 chars.
        if align_col == "":
            align_col = "No Alignment Possible"
        else:
            align_col = align_col[:100]
        if align_row == "":
            align_row = "No Alignment Possible"
        else:
            align_row = align_row[:100]

        return align_row, align_col
