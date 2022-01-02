#!/usr/bin/python3
import math

from CS312Graph import *
import time


class NetworkRoutingSolver:
    def __init__(self):
        pass

    def initializeNetwork(self, network):
        assert (type(network) == CS312Graph)
        self.network = network


    # O(V + E)
    # we have to check at most all of the vertices and all of the edges
    # We can view the number of edges as a constant multiple of the number of vertices
    # Therefore we can O(V) time and space.
    def getShortestPath(self, destIndex):
        self.dest = destIndex

        # Initialization here. O(1) time and space
        path_edges = []
        total_length = 0
        nodes = self.network.nodes
        previous_nodes = self.dijkstras_results["prev"]
        node_distances = self.dijkstras_results["cost"]
        current_index = destIndex

        # This section could be up to O(V + E), or O(V)
        while current_index != self.source:
            node_previous = previous_nodes[current_index]
            if node_previous is None:
                return {'cost': math.inf, 'path': []}
            length_current_edge = self.findCurrentEdgeLength(node_previous, current_index)
            if length_current_edge is None:
                return {'cost': math.inf, 'path': []}
            path_edges.append(
                (nodes[current_index].loc, nodes[node_previous].loc, '{:.0f}'.format(length_current_edge)))
            total_length += length_current_edge
            current_index = node_previous
        return {'cost': node_distances[destIndex], 'path': path_edges}

    # O(1) time and spcae helper function to get the proper corresponding edge that matches (since all vertices have a
    # fixed number of edges)
    def findCurrentEdgeLength(self, current_node, previous_node):
        for edge in self.network.nodes[current_node].neighbors:
            if edge.dest.node_id == previous_node:
                return edge.length

    # This time complexity is the same as dijkstra's
    def computeShortestPaths(self, srcIndex, use_heap=False):
        self.source = srcIndex
        t1 = time.time()

        #runs dijkstra's algorithm and stores the results in a variable in self named "dijkstras_results"
        self.dijkstras_results = self.run_dijkstras(self.network.nodes, self.source, use_heap)

        t2 = time.time()
        return (t2 - t1)

    # If dijkstras is implemented with an array, we see O(V^2) complexity for time and O(V) for space
    # If implemented with a minHeap, we see O(V*log*V) for time and O(V) for space
    def run_dijkstras(self, graph, source_node, use_heap):
        dist = {}
        prev = {}

        # O(V) looping through all vertices
        for index, node in enumerate(graph):
            dist[node.node_id] = math.inf
            prev[node.node_id] = None

        dist[source_node] = 0

        if use_heap:
            p_queue = MinHeapQueue(dist)
        else:
            p_queue = ArrayQueue(source_node, dist)

        iteration_num = 0

        # The while loop condition is O(1) since we are just querying a static variable
        while p_queue.size() != 0:
            # Delete min is the main place where the time complexity changes
            # For a heap, the complexity of delete_min() is O(log(V)). It is called V times, so O(V*log(V)
            # For an array implementation, the complexity is O(V). It is called V times, so O(V^2)
            u = p_queue.delete_min(dist)

            # for every edge E, (We can say V, since E is a constant multiple of V)
            for edge in graph[u].neighbors:
                if (edge.length + dist[u]) < dist[edge.dest.node_id]:
                    dist[edge.dest.node_id] = dist[u] + edge.length
                    prev[edge.dest.node_id] = u

                    # decrease_key() is O(1) for an array because we don't need to do anything
                    # for a minHeap it is log(V)
                    p_queue.decrease_key(edge.dest.node_id, dist[edge.dest.node_id])
            iteration_num += 1
        return {"cost": dist, "prev": prev}


class MinHeapQueue:

    # Time complexity is O(V*log(V)) because we must insert V times at O(log(V)) cost
    # (every vertex that is in our graph)
    # Space complexity is O(V) because we store each vertex
    def __init__(self, nodes):
        self.heap_array = []
        self.distances = []
        self.index_keys = list(nodes)

        for node_id, value in nodes.items():
            self.insert(node_id, value)

    # We can check the size in O(1) time. O(1) space
    def size(self):
        return len(self.heap_array)

    # Insert is O(1) at the end of the array, but we have to call bubble_up which is log(V) therefore, the time
    # complexity is O(log(V))
    # Space Complexity is O(1) since it doesn't grow with the input. Each vertex is the same size.
    def insert(self, node, value):
        self.heap_array.append(node)
        self.distances.append(value)
        self.index_keys[node] = self.size() - 1
        self.bubble_up(self.size() - 1)

    # Time complexity is log(V) because we are recursively traveling through the levels of the tree.
    # Space Complexity is O(1) since nothing is stored
    def bubble_up(self, last_node_num):
        still_bubbling = True
        while still_bubbling and last_node_num != 0:
            if self.distances[last_node_num] < self.get_parent_values(last_node_num):
                parent_index = self.get_parent_index(last_node_num)
                self.swap(last_node_num, parent_index)
                still_bubbling = True
                last_node_num = parent_index
            else:
                still_bubbling = False
        return

    # Time complexity is log(V) because we are recursively traveling down through the levels of the tree.
    # Space Complexity is O(1) since nothing is stored
    def bubble_down(self, node_num):
        while node_num < (len(self.heap_array) // 2):
            min_child, child_index = self.find_min_child_dist_index(node_num)
            if min_child < self.distances[node_num]:
                self.swap(child_index, node_num)
                node_num = child_index
            else:
                return

    # This function is O(1) time complexity since it just looks up values of nodes
    # O(1) space as nothing is stored.
    def find_min_child_dist_index(self, node_num):
        left_child_dist = math.inf
        right_child_dist = math.inf
        if ((2 * node_num) + 1) < len(self.distances):
            left_child_dist = self.distances[(2 * node_num) + 1]
        if ((2 * node_num) + 2) < len(self.distances):
            right_child_dist = self.distances[(2 * node_num) + 2]
        if left_child_dist <= right_child_dist:
            return left_child_dist, (2 * node_num) + 1
        else:
            return right_child_dist, (2 * node_num) + 2

    # This function is simply just lookup, so O(1) time, but we call bubble_up which is log(V) so this time complexity
    # is O(log(V))
    # O(1) space -- nothing is stored and the space remains constant
    def decrease_key(self, node_id, value):
        index_node = self.index_keys[node_id]
        if index_node is None:
            return
        self.distances[index_node] = value
        self.bubble_up(index_node)

    # O(1) time except for bubble_down which is log(V) so O(log(V)) time
    # O(1) space nothing is stored
    def delete_min(self, dist):
        min_node = self.heap_array[0]
        if len(self.heap_array) == 1:
            self.heap_array.pop(0)
            self.distances.pop(0)
            return min_node

        # update the master index map
        self.index_keys[self.heap_array[-1]] = 0
        self.index_keys[self.heap_array[0]] = None

        # update the actual distance and node to do complete the swap
        final_node_dist = self.distances[-1]
        final_node = self.heap_array[-1]
        self.heap_array[0] = final_node
        self.distances[0] = final_node_dist

        # delete last element
        self.heap_array.pop(len(self.heap_array) - 1)
        self.distances.pop(len(self.heap_array) - 1)

        # bubble down to re-order the heap properly, O(logV) time
        self.bubble_down(0)
        return min_node

    # O(1) time
    # O(1) space, nothing stored
    def get_parent_values(self, node_number):
        return self.distances[(node_number - 1) // 2]

    # O(1) time, just lookups and assignments
    # O(1) space, no increase in data stored.  Tree is still the same length as before
    def swap(self, index_child, index_parent):
        # swap the key values
        parent_node = self.heap_array[index_parent]
        child_node = self.heap_array[index_child]
        self.index_keys[child_node] = index_parent
        self.index_keys[parent_node] = index_child

        # swap distance values
        child_node = self.distances[index_child]
        self.distances[index_child] = self.distances[index_parent]
        self.distances[index_parent] = child_node

        # swap the positions in the heap array
        child_node = self.heap_array[index_child]
        self.heap_array[index_child] = self.heap_array[index_parent]
        self.heap_array[index_parent] = child_node

    # O(1) time and space helper function to make the code a little neater.
    def get_parent_index(self, node_num):
        return (node_num - 1) // 2


class ArrayQueue:

    # O(V) time complexity to initialize, because we iterate through all vertices V and insert is O(1)
    # O(V) space complexity as we have to store V vertices with their edges (a constant multiple of the # of vertices)
    def __init__(self, first_node, distances):
        self.array = []
        keys = distances.keys()
        for node_id in keys:
            self.insert(node_id)

    # O(1) time, inserting at the end of an array is constant time
    # O(1) space, nothing more is stored besides the vertex.
    def insert(self, k):
        self.array.append(k)
        return

    # O(1) time and space.  Does nothing in array implementation
    def decrease_key(self, k, val):
        pass

    # Time complexity is O(V) because we most iterate at most through V vertices in our array.
    # Space complexity is O(1) since we aren't storing anything
    def delete_min(self, distances):
        min_value = math.inf
        min_index = 0
        for i in range(len(self.array)):
            if distances[self.array[i]] < min_value:
                min_value = distances[self.array[i]]
                min_index = i
        index_to_return = self.array[min_index]
        self.array.pop(min_index)
        return index_to_return

    # O(1) time and space, constant lookup and nothing is stored.
    def size(self):
        return len(self.array)
