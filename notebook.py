# %% [markdown]
# # 19. Graphs
#

# %% [markdown]
# ## Implementation
#
# Basic operations that we can perform on a graph:
#
# 1. Create a graph
# 2. Display graph vertices
# 3. Display graph edges
# 4. Add a vertex
# 5. Add an edge
#

# %%
# Dict representing a graph, where key is the name of a vertex and values are the names of the vertices it connects to.

graph_elements = {
    "a": ["b", "c"],
    "b": ["a", "d"],
    "c": ["a", "d"],
    "d": ["e"],
    "e": ["d"],
}

print(graph_elements)

# %%


class Graph:
    def __init__(self, gdict=None):
        if gdict is None:
            gdict = []
        self.gdict = gdict

    def get_vertices(self):
        """
        2. Display graph vertices

        - we need to get keys of the graph dict,
          because they represent vertices
        """
        return list(self.gdict.keys())

    def get_edges(self):
        """
        3. Display graph edges

        - we have to find each of the pairs of vertices which have an edge in between them.
        - create an empty list of edges
        - iterate through the edge values associated with each of the vertices
        - build a list containing the distinct group of edges found from the vertices
        """
        edge_names = []
        for vertex in self.gdict:
            for next_vertex in self.gdict[vertex]:
                if {next_vertex, vertex} not in edge_names:
                    edge_names.append({vertex, next_vertex})
        return edge_names

    def add_vertex(self, vertex):
        """
        4. Add a vertex
        - we need to add another additional key to the graph dictionary.
        """
        if vertex not in self.gdict:
            self.gdict[vertex] = []

    def add_edge(self, edge):
        """
        5. Adding an edge
        - treat the new vertex as a tuple
        - validate if the edge is already present
        - if not then add the edge
        """
        edge = set(edge)
        (vrtx1, vrtx2) = tuple(edge)
        if vrtx1 in self.gdict:
            self.gdict[vrtx1].append(vrtx2)
        else:
            self.gdict[vrtx1] = [vrtx2]


# %%

g = Graph(graph_elements)
print(g.get_vertices())
print(g.get_edges())

g.add_vertex("f")
print(g.get_vertices())

g.add_edge({"a", "e"})
g.add_edge({"f", "x"})
print(g.get_edges())

print(g.gdict)

# %%
## Exercises

# %%
"""
GENERATE SHORTEST PATH

    a    b
   / \  / \
  |   \/   \
  |    d    \
  |   /\     \
  |  |_|      |
c|____________|e

"""


graph_dict = {
    "a": ["c"],
    "b": ["d", "e"],
    "c": ["e"],
    "d": ["a", "d", "b"],
    "e": ["b", "c"],
}


def find_shortest_path(graph, start, end, path=None):
    if not path:
        path = []
    path = path + [start]
    if start == end:
        return path
    shortest = None
    for node in graph[start]:
        if node not in path:
            new_path = find_shortest_path(graph, node, end, path)
            if new_path:
                if not shortest or len(new_path) < len(shortest):
                    shortest = new_path
    return shortest


# Find the shortest path between d and c nodes
print(find_shortest_path(graph_dict, "d", "c"))


"""
Flow of execution: 

(d, c, None) 
  (a, c, [d])
    (c, c, [d, a]) -> [d, a, c]   # set as shortest
    
  (b, c, [d])
    (e, c, [d, b])
      (c, c, [d, b, e])  -> [d, a, b, c]  # check if shortest
        
"""

# %%

from collections import defaultdict

"""
A slightly different implementation of a graph. 
It also has a method that allows us to generate all possible paths 
from one vertex to another within the graph.

"""


class Graph:
    def __init__(self, vertices):
        # No. of vertices
        self.v = vertices

        # default dictionary to store graph
        self.graph = defaultdict(list)

    # function to add an edge to graph
    def add_edge(self, u, v):
        self.graph[u].append(v)

    def print_all_paths(self, start, destination):

        # Mark all the vertices as not visited
        visited = [False] * self.v

        # Create an array to store paths
        path = []

        # Call the recursive helper function to print all paths
        self._print_all_paths(start, destination, visited, path)

    def _print_all_paths(self, vertex, destination, visited, path):
        """
        A recursive function to print all paths from 'start' to 'destination'.
        visited[] keeps track of vertices in current path.
        path[] stores actual vertices and path_index is current
        index in path[]
        """
        # Mark the current node as visited and store in path
        visited[vertex] = True
        path.append(vertex)

        if vertex == destination:
            print(path)
        else:
            # Recur for all the vertices adjacent to this vertex
            for i in self.graph[vertex]:
                if not visited[i]:
                    self._print_all_paths(i, destination, visited, path)

        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[vertex] = False


g = Graph(4)
g.add_edge(0, 1)
g.add_edge(0, 2)
g.add_edge(0, 3)
g.add_edge(2, 0)
g.add_edge(2, 1)
g.add_edge(1, 3)


start_vertex = 2
dest_vertex = 3
print(
    "Following are all different paths from vertex {} to vertex {}:".format(
        start_vertex, dest_vertex
    )
)
g.print_all_paths(start_vertex, dest_vertex)
