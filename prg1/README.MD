# DOCUMENTATION 

## Datastructure used in the assignment:
  1. `PLACES`: an `unordered_map` with key is the `PlaceID` and value is a `struct Place`.   
  2. `AREAS`: an `unordered_map` with key is the `AreaID` and value is a `struct Area`. In `struct Area`, there are 2 `vector`(s). `vector region` stores a list of coordinates for that area, and `vector childrenID` stores a list of children areas' IDs.   
  3. An `unordered_set visited` stores IDs of visted areas when implementing operation `all_subareas_in_area`.

## Reasons for choosing the data structures
1. The 2 global `unordered_map`(s):
    - For accessing operations `place_name_type`, `place_coord`, or `area_name` which are called more often than other operations, using `unordered_map` 
allow these operations to run in constant time. There are some rare cases where harsh collision exists, the operation might run in linear time, but still,
these cases are almost negligible. In short, running these operations N times will only cost linear time - O(N), and this can be shown by the running the following performance test:
      ```
      > read "perftest-access.txt"
      ** Commands from 'perftest-access.txt'
      > # Test the performance of stop_name/stop_coord/region_name
      > perftest place_name_type;place_coord;area_name 20 5000 10;30;100;300;1000;3000;10000;30000;100000;300000;1000000
      Timeout for each N is 20 sec. 
      For each N perform 5000 random command(s) from:
      place_name_type place_coord area_name 
      <Result table>
      > 
      ** End of commands from 'perftest-access.txt'
      ```
      Result table

      |       N |    add (sec) |   cmds (sec) |  total (sec) |     N log N |
      | ------- |------------- | ------------ | ------------ | ----------- |
      |      10 |            0 |            0 |            0 | 23.02585093 |
      |      30 |            0 |    0.0009973 |    0.0009973 | 102.0359214 |
      |     100 |            0 |    0.0019946 |    0.0019946 | 460.5170186 |
      |     300 |     0.000997 |    0.0019944 |    0.0029914 | 1711.134742 |
      |    1000 |     0.000997 |            0 |     0.000997 | 6907.755279 |
      |    3000 |    0.0029917 |    0.0019951 |    0.0049868 |  24019.1027 | 
      |   10000 |    0.0099734 |    0.0009973 |    0.0109707 | 92103.40372 |
      |   30000 |    0.0289222 |    0.0009974 |    0.0299196 | 309268.5798 |
      |  100000 |    0.0987364 |    0.0029927 |     0.101729 | 1151292.546 |
      |  300000 |     0.298203 |     0.001995 |     0.300198 | 3783461.326 |
      | 1000000 |      1.02725 |    0.0039889 |      1.03124 | 13815510.56 |

      _I add 1 extra collumn at the end of the table to show asymptotic complexity of the accessing operations runs in linear time._

      With large N, N=10000, 30000, 100000, 300000, 1000000, it can be seen from the table above that the rate of increasing in total absolute time
      is more similar to O(N) than either O(N log N) or O(N^2).  

    - Besides accessing, inserting and removing operations like `add_place`, `add_area`, or `remove_place` also have the same asymptotic time complexity.
    - Other operations using same implementation of the above mentioned operations like `change_place_name`, or `change_place_coord` may also gain the same advantages, and thus have the same asymptotic time.

2. The 2 `struct`(s) as the value according to each key storing in the 2 `unordered_map`(s)
    - Struct is very convenient in storing different properties of an object by grouping all of them together.
    - It is easy to use and implement, no need to worry about private or public interface like in `Class`.

3. The 2 `vector`(s) in the value field of the `struct` Area 
    - The order of elements in those vectors is not important
    - In operation like `get_area_coords` in which the return type is `vector`, using `vector` help the algorithm avoid any extra work to convert from any other container to `vector`.
    - There's no operation requiring these vector to search for element inside it. This is important to consider since searching for element in a `vector` is slow - O(N), compared to `unordered_set` or `set` which run in O(log N) or O(1) on average respectively.

4. `all_subareas_in_area`: Depth-first search and `unordered_set` 
    1. Depth-first search (DFS) algorithm: 
        - One of the most simple, efficient algorithm to find all connected components in a graph. In this assignment, the problem is even simpler. 
Each `Area` have a list of children `Area` and so do for each children `Area`, and no child area can exist in more than 1 list of any children `Area`. DFS uses stack (or vector in 
this assignment) to store the children `Area`(s) and remove it when each of this `Area` is visited. All visited `Area` is add into an `unordered_set` which contains all children `Area`(s).
    2. `unordered_set`: 
        - In DFS, there's a step of checking whether an element has already been visited before by searching that element inside the container containing all children `Area`(s). 
With `unordered_set`, the operation runs in constant time most of the time.

## Possible disadvantages and discussion:
1. For sorting operation like `places_alphabetically` or `places_coord_order`, the time complexity is O(N log N) if `map` is not also used to store elements. However, doing so will make the insertion
operation runs in O(log N) instead of O(1), and running N operations will be O(N log N) instead of O(N), even though when it comes to sorting operation, the operation will run in O(N) instead of O(N log N).
Without any further specification on which operation is more important, not only using some extra `map`(s) doesn't bring any real value in terms of total time complexity, but it also adds some cost on memory consumption. 
2. Some may argue that the value in the 2 `unordered_map` above should be pointer pointing to `struct` objects. Pointer or smart pointer is indeed very helpful in memory management as it allows us to get value by 
refering to a memory location instead of hard copying the whole `struct`. However, in the project assignment, pointers are also used in operation's implementation, especially in the `for` loop. 
Hence, no unnecessary copying exists in the implementation.
3. Operation that can be improved but don't have enough time: 
      - `places_closest_to`: Using heap may help the operation runs in O(N) time instead of O(N log N)
      - `all_subareas_in_area`: Maybe `unordered_set` is unnecessary since the step of checking whether an `Area` has already been visited before may always result in False. This is because no child area can exist in more than 1 list of any children `Area`.
