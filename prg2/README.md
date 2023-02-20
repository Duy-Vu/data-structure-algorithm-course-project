# DOCUMENTATION 

## Datastructure used in the assignment:
  1. `WAYS`: an `unordered_map` with key is the `WayID` and value is a pointer to `struct Way`. In `struct Way`, there is a `vector way_coords` which stores a list of Coord in this way `WaydID`.     
  2. `CROSSROADS`: an `unordered_map` with key is the `Coord` and value is a pointer to `struct Crossroad`. In `struct Crossroad`, there is a `vector ways_from` which stores a list of pair of WaydID and Coord of the next crossroad this way `WaydID` leading to.   


## Reasons for choosing the data structures
1. The 2 global `unordered_map`(s): Fast access time  
  Amortized constant accessing time is needed when implementing other commands like `route_least_crossroads`, `route_shortest_distance`, and `trim_ways` so that it doesn't become the bottleneck of the algorithms. 

2. The 2 `struct`(s):
    - Struct is very convenient in storing different properties of an object by grouping all of them together.
    - It is easy to use and implement, no need to worry about private or public interface like in `Class`.

3. Pointer to struct is used when storing parent way leading to the `Crossroad`, so it doesn't requires look up when we need (reconstructing the path), though the `WayID` is used instead of pointer to `struct Way` since we don't need to look up for it when reconstructing, and we can also use the ID as output value immediately.

4. The `vector` in each `struct`: 
    - `way_coords`: To store the `Coord`(s) of a `Way`, for the command `get_way_coords`
    - `ways_from`: To store neighbors `Crossroad`(s) and `Way`(s) for each `Crossroad`.


## Possible disadvantages and discussion:
1. `vector ways_from` in `struct Crossroad` can store pointers pointing to other neighbor pairs of `Way` and `Crossroad` instead of `WayID` and `Coord` since the performance of the command `ways_from` may not be so important anyway while many algorithms now have to rely on that accessing time of `unordered_map CROSSROADS` and `WAYS` which requires looking up every time, and can be linear in worst case and become the bottleneck of the algorithms.
2. `route_shortest_distance`: A* may be better with euclidian distance is used as heuristic function. 
3. `trim_ways` (Prim's vs Kruskal's): Prim's algorithm is faster with dense graph, i.e many more edges than vertices, wheras Kruskal's performs better in sparse graphs. And it seems that we have more dense graphs than sparse ones.  
