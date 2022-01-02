#!/usr/bin/python3
import traceback
from copy import deepcopy

from PyQt5.QtCore import QLineF, QPointF

import numpy as np
from TSPClasses import *
import heapq


class HeapQueue:

    def __init__(self):
        self.list = []
        self.heapify_list()

    def heapify_list(self):
        heapq.heapify(self.list)

    def heap_push(self, val):
        heapq.heappush(self.list, val)

    def heap_pop(self):
        return heapq.heappop(self.list)

    def get_length(self):
        return len(self.list)


class Node:

    def __init__(self, lower_bound, current_city, cities_to_visit, reduced_matrix_current, node_order, total_cost):
        self.current_city = current_city
        self.cities_to_visit = cities_to_visit
        self.reduced_matrix_current = reduced_matrix_current
        self.node_order = node_order
        self.total_cost = total_cost
        self.lower_bound = lower_bound

    def get_lower_bound(self):
        return self.lower_bound

    def set_lower_bound(self, val):
        self.lower_bound = val

    def get_current_city(self):
        return self.current_city

    def set_current_city(self, val):
        self.current_city = val

    def get_cities_to_visit(self):
        return self.cities_to_visit

    def set_cities_to_visit(self, val):
        self.cities_to_visit = val

    def get_reduced_cost_matrix(self):
        return self.reduced_matrix_current

    def set_reduced_cost_matrix(self, val):
        self.reduced_matrix_current = val

    def get_node_order(self):
        return self.node_order

    def set_node_order(self, val):
        self.node_order = val

    def get_total_cost(self):
        return self.total_cost

    def set_total_cost(self, val):
        self.total_cost = val

    # Criteria to sort in priority queue (sort off the lower bound)
    def __lt__(self, other):
        return self.lower_bound < other.lower_bound


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None
        self.bottom_counter = 0

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    # greedy tour
    def greedy(self, time_allowance=60.0):
        # Initialize the variables we need to report
        results = {}
        cities = self._scenario.getCities()
        num_cities = len(cities)
        found_tour = False
        count = 0
        bssf = None
        start_time = time.time()

        # Start the greedy tour
        while not found_tour and time.time() - start_time < time_allowance:
            best_route = None
            min_total_dist = np.inf

            # Iterate through each city sequentially to find a greedy tour
            for current_city_index in range(num_cities):
                path = []
                remaining_cities = list(range(num_cities))
                path.append(cities[current_city_index])
                remaining_cities.remove(current_city_index)
                while len(path) != num_cities:
                    min_cost = np.inf
                    dest_city = None

                    # always pick the next closest city
                    for city in remaining_cities:
                        new_cost = cities[current_city_index].costTo(cities[city])
                        if new_cost < min_cost:
                            min_cost = new_cost
                            dest_city = city
                    if dest_city is None:  # tie breaking condition, just continue and pick the first element.
                        dest_city = remaining_cities[0]
                    path.append(cities[dest_city])
                    remaining_cities.remove(dest_city)
                    current_city_index = dest_city
                count += 1
                new_route = TSPSolution(path)
                new_total_dist = new_route.cost
                if new_total_dist < min_total_dist:
                    min_total_dist = new_total_dist
                    best_route = new_route
            bssf = best_route

            # We must know if we found a valid tour before we return
            if bssf.cost < np.inf:
                found_tour = True
            else:
                bssf.cost = np.inf

        end_time = time.time()

        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        results['bssf'] = bssf
        results['time'] = end_time - start_time
        results['cost'] = bssf.cost
        results['count'] = count
        results['soln'] = bssf

        return results

    # Branch and Bound Tour
    # Worst case O(n!) time, but will likely be much better due to pruning
    def branchAndBound(self, time_allowance=60.0):

        # Set some initial variables we need to keep track of
        num_solutions = 0
        bssf_updates = 0
        pruned_subproblems = 0

        # we will by default have one item on the queue at least and 1 state
        max_priority_queue_length = 1
        total_states_created = 1

        # Get the cities from the GUI
        cities = self._scenario.getCities()
        self.cities = cities

        # Get a starting value from running a greedy algorithm
        # My greedy implementation is at worst O(n^2) time and O(n) space to hold the cities
        bssf = self.greedy(time_allowance=time_allowance)['soln']
        self.current_lowest_cost = bssf.cost

        # get the initial reduced cost matrix and a starting lower bound
        lower_bound, initial_reduced_matrix = self.make_starting_matrix(cities)

        # make a node representing the first subproblem
        # Finding a reduced cost matrix is O(n^2) time and space
        first_node = Node(lower_bound, cities[0], cities[1:], initial_reduced_matrix, [cities[0].get_index()],
                          lower_bound)

        # Create the priority queue and push the first subproblem onto the queue
        priority_queue = HeapQueue()
        priority_queue.heap_push(first_node)

        time_start = time.time()

        # While there are still items on the queue and we are under the time limit
        while priority_queue.get_length() and time.time() - time_start < time_allowance:

            # update the max length that we achieved on the priority queue
            max_priority_queue_length = max(priority_queue.get_length(), max_priority_queue_length)

            # pop is O(logn) time because we have to bubble up/down (reheapify)
            # O(1) space
            popped_node = priority_queue.heap_pop()
            # O(1) lookup
            if popped_node.get_total_cost() < self.current_lowest_cost:
                # Worst case we have to visit every city -- O(n) time
                for city in popped_node.get_cities_to_visit():
                    if popped_node.get_current_city().costTo(city) != np.inf:

                        # reduced cost matrix is O(n^2) time and space
                        new_node = self.make_matrix(city, popped_node.get_reduced_cost_matrix(),
                                                    popped_node)

                        # if we have run out of cities to visit
                        if len(new_node.get_cities_to_visit()) == 0:
                            # This function is at worst O(n) for the number of cities
                            path = self.get_cities_from_indexes(new_node.get_node_order())
                            bssf = TSPSolution(path)

                            if bssf.cost < self.current_lowest_cost:
                                self.current_lowest_cost = min(bssf.cost, self.current_lowest_cost)
                                bssf_updates += 1
                            num_solutions += 1
                            print("*** Solution found! Score: {} ***".format(self.current_lowest_cost))

                        # else we still have cities to visit
                        else:
                            # If the node has a lower bound that is still interesting
                            if new_node.get_total_cost() < self.current_lowest_cost:
                                # O(logn) time for pushing a new node to the priority queue
                                priority_queue.heap_push(new_node)
                                total_states_created += 1

                            # prune (do not add to our queue)
                            else:
                                pruned_subproblems += 1
                                total_states_created += 1

            # Prune (this is a node we used to think is interesting, but is no longer)
            else:
                pruned_subproblems += 1
                total_states_created += 1

        end_time = time.time()

        results = {'soln': bssf, 'max': max_priority_queue_length, 'total': total_states_created,
                   'pruned': pruned_subproblems, 'cost': self.current_lowest_cost, 'time': end_time - time_start,
                   'count': num_solutions}

        return results

    # For some reason, keeping node indices rather than the city objects as the array is much faster (likely due to
    # less work when deep copying)
    # O(n) time and space
    def get_cities_from_indexes(self, indices):
        cities = []
        for i in indices:
            cities.append(self.cities[i])
        return cities

    # O(n^2) time and space to visit and update the cells of an n by n matrix where n is the num of cities
    def make_starting_matrix(self, list_of_cities):
        # O(n^2) time and space to make and initialize an n by n matrix
        matrix = np.full((len(list_of_cities), len(list_of_cities)), fill_value=np.inf)
        for i, city in enumerate(list_of_cities):
            for j, dest_city in enumerate(list_of_cities):
                if i == j:
                    continue  # skip
                distance = city.costTo(dest_city)
                matrix[i][j] = distance

        reduction_total = 0

        # O(n^2) time to perform matrix reduction
        for row in range(matrix.shape[0]):
            row_min = np.min(matrix[row])
            matrix[row] = matrix[row] - row_min
            reduction_total += row_min

        for col in range(matrix.shape[1]):
            col_min = np.min(matrix[:, col])
            matrix[:, col] = matrix[:, col] - col_min
            reduction_total += col_min

        return reduction_total, matrix

    # O(1) time and space
    def get_value(self, score, cities_visited):
        if len(cities_visited):
            return score / len(cities_visited)
        else:
            return score  # we don't care about the first level, we still have to look at it.

    # O(n^2) time and space to visit and update the cells of an n by n matrix where n is the num of cities
    def make_matrix(self, next_city_to_visit, current_matrix, current_node):
        node = deepcopy(current_node)
        new_matrix = current_matrix.copy()

        sum_to_reduce = 0
        initial_cost_to_city = new_matrix[node.get_current_city().get_index()][next_city_to_visit.get_index()]

        # O(1) operations to visit the column and row that we are working with
        new_matrix[node.get_current_city().get_index()] = np.inf
        new_matrix[:, next_city_to_visit.get_index()] = np.inf
        new_matrix[next_city_to_visit.get_index()][node.get_current_city().get_index()] = np.inf

        # the following is the actual matrix reduction, which is O(n^2) time and O(1) space since it uses the same
        # matrix already in use
        for row in range(new_matrix.shape[0]):
            min_for_row = np.min(new_matrix[row])
            if np.isinf(min_for_row):
                continue
            new_matrix[row] = new_matrix[row] - min_for_row
            sum_to_reduce += min_for_row

        for col in range(new_matrix.shape[1]):
            min_for_col = np.min(new_matrix[:, col])
            if np.isinf(min_for_col):
                continue
            new_matrix[:, col] = new_matrix[:, col] - min_for_col
            sum_to_reduce += min_for_col

        # O(n) time to visit to all of the cities to delete the city we have now visited
        unvisited_cities = node.get_cities_to_visit()
        index_to_delete = None
        for i, city in enumerate(unvisited_cities):
            if city.get_index() == next_city_to_visit.get_index():
                index_to_delete = i
                break
        del unvisited_cities[index_to_delete]

        # O(1) mathematical operations and lookups to calculate new lower bound
        new_cost = node.get_total_cost() + initial_cost_to_city + sum_to_reduce

        # return the new subproblem
        new_node = Node(self.get_value(node.get_total_cost(), node.get_node_order()), next_city_to_visit,
                        unvisited_cities, new_matrix, node.get_node_order() + [next_city_to_visit.get_index()],
                        new_cost)

        return new_node

    def fancy(self, time_allowance=60.0):
        pass
