from PyQt5.QtCore import QLineF, QPointF, QObject

import time

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = .1


# Global useful functions
# Time Complexity: O(1) because the function is just mathematical operations and has no loops
# Space Complexity: O(1) constant because no further information is required to be stored
def get_slope(point_one, point_two):
    return (point_two.y() - point_one.y()) / (point_two.x() - point_one.x())

# Time Complexity: O(1)
# Space Complexity: O(1) this is a simple lookup of a value within an object
def sort_key(point):
    return point.x()


class ConvexHullSolver(QObject):
    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()

        # SORT THE POINTS HERE USING PYTHON BUILT IN SORT FUNCTION
        # Defined a sort_key to sort the point objects by globally.  We want to sort by x value.
        # Time Complexity: 0(n * log(n)) (This sort function implements Python Timsort which is reported as O(n*log(n))
        # Space Complexity: O(n) From several sources on the net, it only needs space that scales linearly with n
        points.sort(key=sort_key)

        t2 = time.time()
        self.showText('Time Elapsed (Sorting): {:3.3f} sec'.format(t2 - t1))
        t3 = time.time()

        # CALL TO THE DIVIDE AND CONQUER FUNCTION
        final_hull = self.divide_and_conquer(points, pause, view)
        t4 = time.time()

        # CONVERT LIST OF POINTS TO A LIST OF LINES
        final_line_list = []
        for i in range(len(final_hull)):
            final_line_list.append(QLineF(final_hull[i], final_hull[(i + 1) % len(final_hull)]))
        self.showHull(final_line_list, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.12f} sec'.format(t4 - t3))

    # DIVIDE AND CONQUER FUNCTION
    # Overall Time Complexity of O(n*log(n))
    # Overall Space Complexity of O(n) // CHECK ON THIS
    # Master Theorem: T(n) = 2T(n/2) + O(n)
    # Branching factor = 2, half the work done per level, O(n) merge
    def divide_and_conquer(self, points, pause, view):
        # Handle the base case of 1 point
        length = len(points)  # O(1)
        if length == 1:  # 0(1)
            return points

        # Recursively split into left and right branches
        midpoint = length // 2  # O(1)
        left_branch = self.divide_and_conquer(points[:midpoint], pause, view)  # O(n) time complexity
        right_branch = self.divide_and_conquer(points[midpoint:], pause, view)  # O(n) time complexity

        # Simplest merge case: merge two branches of size 1
        if len(left_branch) == 1 and len(right_branch) == 1:
            left_branch.append(right_branch[0])
            return left_branch

        # find the rightmost point of the left branch and the leftmost point of the right branch
        left_pt = left_branch.index(max(left_branch, key=sort_key))  # O(n) has to iterate through at most n/2 points
        right_pt = right_branch.index(min(right_branch, key=sort_key)) # O(n) "   "

        # METHOD TO FIND UPPER TANGENT POINTS:
        # O(n) time
        # O(n) space (we have to go through and check at most n points as possible tangent points
        # Start the rightmost point of the left branch and leftmost point of the right branch.
        current_left_upper = left_pt
        current_right_upper = right_pt

        # Get the current slope
        current_slope = get_slope(left_branch[left_pt], right_branch[right_pt])  # O(1)
        left_changed = True
        right_changed = True
        # While we've changed either the left or right points
        while left_changed or right_changed:
            left_changed = False
            right_changed = False
            # Move the left point counterclockwise while the slope becomes more negative
            while True:
                new_slope = get_slope(left_branch[(current_left_upper - 1 % len(left_branch))],
                                      right_branch[current_right_upper])
                if new_slope < current_slope:
                    current_slope = new_slope
                    left_changed = True
                    current_left_upper = (current_left_upper - 1) % len(left_branch)
                else:
                    break
            # Move the right point clockwise while the slope becomes more positive
            while True:
                new_slope = get_slope(left_branch[current_left_upper],
                                      right_branch[(current_right_upper + 1) % len(right_branch)])
                if new_slope > current_slope:
                    current_slope = new_slope
                    right_changed = True
                    current_right_upper = (current_right_upper + 1) % len(right_branch)
                else:
                    break

        # METHOD TO FIND LOWER TANGENT POINTS
        # O(n) time
        # O(n) space (we have to go through and check at most n points as possible tangent points)
        current_left_lower = left_pt
        current_right_lower = right_pt
        # Get the new slope
        current_slope_lower = get_slope(left_branch[left_pt], right_branch[right_pt])
        left_changed = True
        right_changed = True
        while left_changed or right_changed:
            left_changed = False
            right_changed = False
            # Move the left point clockwise while the the slope becomes more positive
            while True:
                new_slope = get_slope(left_branch[(current_left_lower + 1) % len(left_branch)],
                                      right_branch[current_right_lower])
                if new_slope > current_slope_lower:
                    current_slope_lower = new_slope
                    left_changed = True
                    current_left_lower = (current_left_lower + 1) % len(left_branch)
                else:
                    break
            # Move the right point counterclockwise while the slope becomes more negative.
            while True:
                new_slope = get_slope(left_branch[current_left_lower],
                                      right_branch[(current_right_lower - 1) % len(right_branch)])
                if new_slope < current_slope_lower:
                    current_slope_lower = new_slope
                    right_changed = True
                    current_right_lower = (current_right_lower - 1) % len(right_branch)
                else:
                    break

        # Step to show recursive hulls and tangent lines
        # O(n)
        if pause:
            self.iterate_recursive_steps(left_branch, right_branch, current_left_upper, current_left_lower,
                                         current_right_upper, current_right_lower)

        # MERGE STEP:
        # Create a new hull by moving through the appropriate indices of the array.
        # O(n) time because we are iterating through at most every point given
        # O(n) space because we have to store at most n total points
        return_hull = []
        point_to_add = current_left_lower  # start with the lower left tangent point
        return_hull.append(left_branch[point_to_add])

        # move clockwise around the left branch (more positive) until we run into the upper left tangent point
        while point_to_add != current_left_upper:
            point_to_add = (point_to_add + 1) % len(left_branch)
            return_hull.append(left_branch[point_to_add])

        # manually add the index of the upper right tangent point
        point_to_add = current_right_upper
        return_hull.append(right_branch[point_to_add])

        # move clockwise around the right branch until we run into the lower right tangent point
        while point_to_add != current_right_lower:
            point_to_add = (point_to_add + 1) % len(right_branch)
            return_hull.append(right_branch[point_to_add])

        return return_hull

    def iterate_recursive_steps(self, left_branch, right_branch, left_upper, left_lower, right_upper, right_lower):
        left_branch_hull = []
        right_branch_hull = []
        for i in range(len(left_branch)):
            left_branch_hull.append(QLineF(left_branch[i], left_branch[(i + 1) % len(left_branch)]))
        for i in range(len(right_branch)):
            right_branch_hull.append(QLineF(right_branch[i], right_branch[(i + 1) % len(right_branch)]))
        upper_tangent_line = QLineF(left_branch[left_upper], right_branch[right_upper])
        lower_tangent_line = QLineF(left_branch[left_lower], right_branch[right_lower])
        self.showHull(left_branch_hull, RED)
        self.showHull(right_branch_hull, RED)
        self.showTangent([upper_tangent_line, lower_tangent_line], BLUE)
        self.eraseHull(left_branch_hull)
        self.eraseHull(right_branch_hull)
        self.eraseTangent([upper_tangent_line, lower_tangent_line])
