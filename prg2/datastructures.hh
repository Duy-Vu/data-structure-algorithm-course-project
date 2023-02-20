// Datastructures.hh

#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH
#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <cmath>
#include <utility>
#include <limits>
#include <functional>

// Types for IDs
using PlaceID = long long int;
using AreaID = long long int;
using Name = std::string;
using WayID = std::string;
using Distance = int;       // Type for a distance (in metres)

// Return values for cases where required thing was not found
PlaceID const NO_PLACE = -1;
AreaID const NO_AREA = -1;
WayID const NO_WAY = "!!No way!!";

// Return value for cases where integer values were not found
int const NO_VALUE = std::numeric_limits<int>::min();
int const MAX_VALUE = std::numeric_limits<int>::max();

// Return value for cases where Duration is unknown
Distance const NO_DISTANCE = NO_VALUE;

// Return value for cases where name values were not found
Name const NO_NAME = "!!NO_NAME!!";

// Enumeration for different place types
enum class PlaceType { OTHER=0, FIREPIT, SHELTER, PARKING, PEAK, BAY, AREA, NO_TYPE };

// Enumeration to mark if a node has not been discovered, been discorved but not yet processed, or fully processed
enum class ColorType { WHITE=0, GRAY, BLACK };

// Type for a coordinate (x, y)
struct Coord
{
    int x;
    int y;
};

struct Place
{
    std::string name;
    Coord position;
    PlaceType placetype;
};

struct Area
{
    std::string name;
    std::vector<Coord> region;
    AreaID parentID;
    std::vector<AreaID> childrenID;
};

struct Way
{
    std::vector<Coord> way_coords;
    Distance length;
    bool mst;   // for minimum spanning tree
};

struct Crossroad
{
    Coord coord;
    std::vector<std::pair<WayID, Coord>> ways_from; // Get all neighbor paths
    std::pair<Crossroad*, WayID> way_to;            // For backtrack the path
    Distance distance_from_start;
    ColorType color;
};

// Define the coordinate of the origin
Coord const ORIGIN = Coord{0,0};

// Return value for cases where coordinates were not found
Coord const NO_COORD = {NO_VALUE, NO_VALUE};

// Define - for calculating Euclidean distance between 2 Coord
inline float operator-(Coord c1, Coord c2)
{
    return sqrt( pow(c1.x-c2.x, 2) + pow(c1.y-c2.y, 2) );
}

// Define < for Coord so that it can be used for comparison
inline bool operator<(Coord c1, Coord c2)
{
    double dist_from_org1 = c1 - ORIGIN;
    double dist_from_org2 = c2 - ORIGIN;
    if (dist_from_org1 == dist_from_org2) {
        return c1.y < c2.y;
    } else {
        return dist_from_org1 < dist_from_org2;
    }
}

// Defining == and hash function for Coord so that it can be used
// as key for std::unordered_map/set, if needed
inline bool operator==(Coord c1, Coord c2) { return c1.x == c2.x && c1.y == c2.y; }
inline bool operator!=(Coord c1, Coord c2) { return !(c1==c2); }

struct CoordHash
{
    std::size_t operator()(Coord xy) const
    {
        auto hasher = std::hash<int>();
        auto xhash = hasher(xy.x);
        auto yhash = hasher(xy.y);
        // Combine hash values (magic!)
        return xhash ^ (yhash + 0x9e3779b9 + (xhash << 6) + (xhash >> 2));
    }
};


class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    // Copy of phase 1 operations

    // Estimate of performance: O(1)
    // The size() method runs in constant time.
    int place_count();

    // Estimate of performance: O(n)
    // Going through each pair of the 2 unordered maps to delete.
    void clear_all();

    // Estimate of performance: O(n)
    // Going through each pair of the unordered map and get the key of each pair.
    std::vector<PlaceID> all_places();

    // Estimate of performance: average O(1), worst O(n)
    // Adding element to the unordered map takes constant time in most cases,
    // but in the worst case where collision exists, it can take linear time.
    bool add_place(PlaceID id, Name const& name, PlaceType type, Coord xy);

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take constant
    // time and retrieving it (if it can be found) also takes constant time.
    // However, if hash collision exists, it takes linear time to either find
    // or access the element.
    std::pair<Name, PlaceType> get_place_name_type(PlaceID id);  // Performance is critical

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take constant
    // time and retrieving it (if it can be found) also takes constant time.
    // However, if hash collision exists, it takes linear time to either find
    // or access the element.
    Coord get_place_coord(PlaceID id);  // Performance is critical

    // Implement the operations below only after implementing the ones above

    // Estimate of performance: O(n log n)
    // Creating a return vector containing all ID takes linear time, and sorting
    // it using std::sort algorithm which has time complexity of O(n log n)
    std::vector<PlaceID> places_alphabetically();

    // Estimate of performance: O(n log n)
    // Creating a return vector containing all ID takes linear time, and sorting
    // it using std::sort algorithm which has time complexity of O(n log n)
    std::vector<PlaceID> places_coord_order();

    // Estimate of performance: O(n)
    // Going through each element in the unordered map and
    // checking if the element has the same specified name property
    std::vector<PlaceID> find_places_name(Name const& name);    // Performance is not critical

    // Estimate of performance: O(n)
    // Going through each element in the unordered map and
    // checking if the element has the same specified name property
    std::vector<PlaceID> find_places_type(PlaceType type);      // Performance is not critical

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take
    // constant time and retrieving it to modify the field 'name' (if it
    // can be found) also takes constant time. However, if hash collision
    // exists, it takes linear time to either find or access the element.
    bool change_place_name(PlaceID id, Name const& newname);

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take
    // constant time and retrieving it to modify the field 'position' (if
    // it can be found) also takes constant time. However, if hash collision
    // exists, it takes linear time to either find or access the element.
    bool change_place_coord(PlaceID id, Coord newcoord);

    // Implement the operations below only after implementing the ones above

    // Estimate of performance: average O(1), worst O(n)
    // Adding element to the unordered map takes constant time in most cases,
    // but in the worst case where collision exists, it can take linear time.
    bool add_area(AreaID id, Name const& name, std::vector<Coord> coords);

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take constant
    // time and retrieving it (if it can be found) also takes constant time.
    // However, if hash collision exists, it takes linear time to either find
    // or access the element.
    Name get_area_name(AreaID id);      // Performance is critical

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take constant
    // time and retrieving it (if it can be found) also takes constant time.
    // However, if hash collision exists, it takes linear time to either find
    // or access the element.
    std::vector<Coord> get_area_coords(AreaID id);

    // Estimate of performance: O(n)
    // Going through each pair of the unordered map and get the key of each pair.
    std::vector<AreaID> all_areas();

    // Estimate of performance: average O(1), worst O(n).
    // Finding and accessing either struct 'Area' to modify the field 'childrenID'
    // or 'parentID' in unordered map takes the same amount of constant time. However,
    // if hash collision exists, it takes linear time to do so.
    bool add_subarea_to_area(AreaID id, AreaID parentid);

    // Estimate of performance: O(n) in worst and average cases, and O(1) in best case
    // From the current area with given ID, find the parent area and its parent, and so
    // on until we cannot find parent area anymore. The worst case happens when all areas
    // are ancestor areas of the area we are asking for. Best case happens when the given
    // ID is not belong to any area or the area doesn't have any parent area.
    std::vector<AreaID> subarea_in_areas(AreaID id);

    // Estimate of performance: O(1)
    // Do nothing
    void creation_finished();

    // Estimate of performance: average and worst O(n), best O(1)
    // Using Depth-First search, continuously adding all direct children areas to the stack
    // and poping them out to search for the next children areas and marking as visited. Note
    // that after adding all possible children area, in the rest of iterations of the while
    // loop, only steps of popping and marking as visited takes place. The time complexity of
    // DFS is O(V + E) where V are numbers of vertices (areas in our cases), and E are edges.
    // But since no area has more than 1 parent area, max(E) = V - 1, and thus the time complexity
    // DFS is O(n). Creating a return vector from the unordered set of visited areas also takes
    // linear time. The best case happens when the asking ID is not found or the area doesn't have
    // any child area.
    std::vector<AreaID> all_subareas_in_area(AreaID id);

    // Estimate of performance: O(n log n)
    // Creating a return vector of areas' IDs satisfying the condition for same
    // asking type takes linear time, and sorting it by std::sort takes O(n log n).
    std::vector<PlaceID> places_closest_to(Coord xy, PlaceType type);

    // Estimate of performance: average O(1), worst O(n)
    // Finding and removing element in the unordered map takes constant time on average,
    // but in the worst case where collision exists, it can take linear time.
    bool remove_place(PlaceID id);      // Performance is not critical

    // Estimate of performance: average and worst O(n), best O(1)
    // Getting ancestor vector by calling subarea_in_areas() twice takes linear time, and
    // then checking for the first diffrent ancestor from those 2 vector from the last position
    // also takes linear time. The best case happens when either areas doesn't have any parent
    // area or the highest ancestors of those areas are different.
    AreaID common_area_of_subareas(AreaID id1, AreaID id2);


    // Phase 2 operations

    // Estimate of performance: O(n)
    // Going through each pair of the unordered map and get the key of each pair.
    std::vector<WayID> all_ways();

    // Estimate of performance:
    // Short rationale for estimate:
    bool add_way(WayID id, std::vector<Coord> coords);

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take constant
    // time and accessing its attribute (if it can be found) also takes constant time.
    // However, if hash collision exists, it takes linear time to either find or access
    // the element.
    std::vector<std::pair<WayID, Coord>> ways_from(Coord xy);

    // Estimate of performance: average O(1), worst O(n)
    // With unordered map, finding if the element based on a key take constant
    // time and accessing its attribute (if it can be found) also takes constant time.
    // However, if hash collision exists, it takes linear time to either find access
    // the element.
    std::vector<Coord> get_way_coords(WayID id);

    // Estimate of performance: O(n)
    // Going through each pair of the 2 unordered maps to delete.
    void clear_ways();

    // Estimate of performance: O(V + E)
    // The algorithm calls route_least_crossroads function which its performance is explained below
    std::vector<std::tuple<Coord, WayID, Distance>> route_any(Coord fromxy, Coord toxy);

    // Non-compulsory operations

    // Estimate of performance: best O(1), worst O(n)
    // Finding a way in unordered_map WAYS to delete it takes constant time, but then finding a way
    // which is in a pair of the vector ways_to of 2 crossroads takes linear time in worst cases,
    // though in real cases, the vectors may not contain many ways, and thus this finding part is
    // nearly constant.
    bool remove_way(WayID id);

    // Estimate of performance: O(V + E)
    // The algorithm is BFS. The while loop go through each vertex in the queue, in total V nodes,
    // and in processing step, iteratively add neighbor nodes color as WHITE only (to avoid add vertex
    // that are already in the queue), in total E times. The constructing path part is linear, so it
    // has the same asymptotic performance with the algorithm's.
    std::vector<std::tuple<Coord, WayID, Distance>> route_least_crossroads(Coord fromxy, Coord toxy);

    // Estimate of performance: O(V + E)
    // The algorithm is DFS. The while loop goes through each vertex in the stack, O(V) times, and each
    // vertex are added to to the stack only when vertices are WHITE (same as BFS, to avoid add same node
    // to the stack), O(E) times. The constructing cycle part is linear, so it has the same asymptotic
    // performance with the algorithm's.
    std::vector<std::tuple<Coord, WayID>> route_with_cycle(Coord fromxy);

    // Estimate of performance: O((V + E)log(V))
    // It is Dijkstra's algorithm. The while loop makes at most O(V) times and the for loop at most O(E)
    // times in total, and on each round of the while and the for loop, nodes are added and sorted, so it
    // takes O(log(V)).
    std::vector<std::tuple<Coord, WayID, Distance>> route_shortest_distance(Coord fromxy, Coord toxy);

    // Estimate of performance: O((V + E)log(V))
    // It's a Prim’s algorithm to solve minimum spanning tree(s) problem. The outermost for loop is just to
    // make sure we can process all nodes even when the graph is disconnected. The while loop makes at most
    // O(V) times and the for loop at most O(E) times in total, and on each round of the while and the for
    // loop, nodes are added and sorted, so it takes O(log(V)). And the final trimming part takes O(E) time.
    Distance trim_ways();

private:
    std::unordered_map<PlaceID,Place> PLACES;
    std::unordered_map<AreaID, Area> AREAS;

    std::unordered_map<WayID, Way*> WAYS;
    std::unordered_map<Coord, Crossroad*, CoordHash> CROSSROADS;

    // Estimate of performance: O(1)
    // The size() method runs in constant time.
    int way_count();

    // Estimate of performance: O(n)
    // Calculate the Euclidean distance for each pair of coordinates in the list of coordinates of the way
    int calculate_distance(const std::vector<Coord>& coords);

    // Estimate of performance: average O(1), worst O(n)
    // Finding or accessing in unordermap takes constant time on average.
    // If hash collision exists, it takes linear time to do the same thing.
    void add_crossroads(const Coord& coord1, const std::pair<WayID, Coord>& way_to);

    // Estimate of performance: O(n)
    // Get the parent Crossroad from the current Crossroad iteratively until the parent Crossroad is the start Crossroad
    void construct_path(Crossroad* cur_crossroad, std::vector<std::tuple<Coord, WayID, Distance>>& path);

    // Estimate of performance: O(n)
    // Get the parent Crossroad from the current Crossroad iteratively until the parent Crossroad is the start Crossroad
    void construct_cycle(Crossroad* cur_crossroad, std::vector<std::tuple<Coord, WayID>>& path);
};

#endif // DATASTRUCTURES_HH
