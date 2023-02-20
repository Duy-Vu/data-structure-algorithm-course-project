// Datastructures.cc

#include "datastructures.hh"
#include <random>
#include <unordered_set>
#include <queue>
#include <stack>
#include <algorithm>
#include <cmath>
#include <QtDebug>

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

Datastructures::Datastructures()
{
    PLACES = {};
    AREAS = {};
    WAYS = {};
    CROSSROADS = {};
}
Datastructures::~Datastructures()
{
    clear_all();
    clear_ways();
}

// Copy of phase 1 operations

int Datastructures::place_count()
{
    return PLACES.size();
}

void Datastructures::clear_all()
{
    PLACES.clear();
    AREAS.clear();
}

std::vector<PlaceID> Datastructures::all_places()
{
    std::vector<PlaceID> places;
    places.reserve(place_count());  // Pre-allocate enough space for the return vector

    for(auto const& pair_place : PLACES) {
        places.push_back(pair_place.first);
    }
    return places;
}

bool Datastructures::add_place(PlaceID id, const Name& name, PlaceType type, Coord xy)
{
    auto [it, insert] = PLACES.emplace(id, Place{name, xy, type});
    return insert;
}

std::pair<Name, PlaceType> Datastructures::get_place_name_type(PlaceID id)
{
    if (PLACES.find(id) != PLACES.end()) {
        return {PLACES.at(id).name, PLACES.at(id).placetype};
    } else {
        return {NO_NAME, PlaceType::NO_TYPE};
    }
}

Coord Datastructures::get_place_coord(PlaceID id)
{
    if (PLACES.find(id) != PLACES.end()) {
        return PLACES.at(id).position;
    } else {
        return NO_COORD;
    }
}

bool Datastructures::add_area(AreaID id, const Name &name, std::vector<Coord> coords)
{
    auto [it, insert] = AREAS.emplace(id, Area{name, coords, NO_AREA, {}});
    return insert;
}

Name Datastructures::get_area_name(AreaID id)
{
    if (AREAS.find(id) != AREAS.end()) {
        return AREAS.at(id).name;
    } else {
        return NO_NAME;
    }
}

std::vector<Coord> Datastructures::get_area_coords(AreaID id)
{
    if (AREAS.find(id) != AREAS.end()) {
        return AREAS.at(id).region;
    }
    return {NO_COORD};
}

void Datastructures::creation_finished()
{
}

std::vector<PlaceID> Datastructures::places_alphabetically()
{
    std::vector<PlaceID> places;
    places.reserve(place_count());

    for (auto const& pair_place : PLACES) {
        places.push_back(pair_place.first);
    }
    std::sort(places.begin(), places.end(), [this](PlaceID id1, PlaceID id2){   // Capture this pointer to access variable PLACES
        return PLACES.at(id1).name < PLACES.at(id2).name;
    });
    return places;
}

std::vector<PlaceID> Datastructures::places_coord_order()
{
    std::vector<PlaceID> places;
    places.reserve(place_count());

    for (auto const& pair_place : PLACES) {
        places.push_back(pair_place.first);
    }
    std::sort(places.begin(), places.end(), [this](PlaceID id1, PlaceID id2){
        return PLACES.at(id1).position < PLACES.at(id2).position;
    });
    return places;
}

std::vector<PlaceID> Datastructures::find_places_name(Name const& name)
{
    std::vector<PlaceID> places = {};
    for(auto const& pair_place : PLACES) {
        if (pair_place.second.name == name) {
            places.push_back(pair_place.first);
        }
    }
    return places;
}

std::vector<PlaceID> Datastructures::find_places_type(PlaceType type)
{
    std::vector<PlaceID> places = {};
    for(auto const& pair_place : PLACES) {
        if (pair_place.second.placetype == type) {
            places.push_back(pair_place.first);
        }
    }
    return places;
}

bool Datastructures::change_place_name(PlaceID id, const Name& newname)
{
    if (PLACES.find(id) != PLACES.end()) {
        PLACES.at(id).name = newname;
        return true;
    }
    return false;
}

bool Datastructures::change_place_coord(PlaceID id, Coord newcoord)
{
    if (PLACES.find(id) != PLACES.end()) {
        PLACES.at(id).position = newcoord;
        return true;
    }
    return false;
}

std::vector<AreaID> Datastructures::all_areas()
{
    std::vector<AreaID> areas = {};
    areas.reserve(AREAS.size());

    for(auto const& pair_area : AREAS) {
        areas.push_back(pair_area.first);
    }
    return areas;
}

bool Datastructures::add_subarea_to_area(AreaID id, AreaID parentid)
{
    if (AREAS.find(id) != AREAS.end() && AREAS.find(parentid) != AREAS.end()) {
        Area child = AREAS.at(id);
        Area parent = AREAS.at(parentid);
        if (child.parentID == NO_AREA) {
            AREAS.at(id).parentID = parentid;
            AREAS.at(parentid).childrenID.push_back(id);
            return true;
        }
    }
    return false;
}

std::vector<AreaID> Datastructures::subarea_in_areas(AreaID id)
{
    // If area doesn't exist
    if (AREAS.find(id) == AREAS.end()) {
        return {NO_AREA};
    }

    // Otherwise
    std::vector<AreaID> all_parents_areas = {};
    AreaID parent_area = AREAS.at(id).parentID;
    while(parent_area != NO_AREA) {
        all_parents_areas.push_back(parent_area);
        parent_area = AREAS.at(parent_area).parentID;     // copy parentID to search for parent of this parent
    }
    return all_parents_areas;
}

std::vector<PlaceID> Datastructures::places_closest_to(Coord xy, PlaceType type)
{
    std::vector<PlaceID> places = {};

    if (type == PlaceType::NO_TYPE) {   // If type is not given
        for (auto const& pair_place : PLACES) {
            places.push_back(pair_place.first);
        }
    }
    else{   // If type is given
        for (auto const& pair_place : PLACES) {
            if (pair_place.second.placetype == type) {
                places.push_back(pair_place.first);
            }
        }
    }

    std::sort(places.begin(), places.end(), [this, &xy](PlaceID id1, PlaceID id2) {
        double d1 = PLACES.at(id1).position - xy;
        double d2 = PLACES.at(id2).position - xy;

        if (d1 == d2) {
           return PLACES.at(id1).position.y < PLACES.at(id2).position.y;
        } else {
            return d1 < d2;
        }
    });

    unsigned int max_size = 3;
    if (places.size() <= max_size) {
        return places;
    } else {
        return std::vector<PlaceID>(places.begin(), places.begin() + max_size);
    }
}

bool Datastructures::remove_place(PlaceID id)
{
    if (PLACES.find(id) != PLACES.end()) {
        PLACES.erase(id);
        return true;
    }
    return false;
}

std::vector<AreaID> Datastructures::all_subareas_in_area(AreaID id)
{
    if (AREAS.find(id) == AREAS.end()){
        return {NO_AREA};
    }

    // Using DFS to search all children areas
    std::unordered_set<AreaID> visited = {};
    std::vector<AreaID> stackID = {id};
    while (!stackID.empty()) {
        AreaID child_id = stackID.back();
        stackID.pop_back();     // O(1)
        if (visited.find(child_id) == visited.end()) {
            // Label as visited
            visited.emplace(child_id);
            // Add to stack
            auto first_child = AREAS.at(child_id).childrenID.begin();
            auto last_child = AREAS.at(child_id).childrenID.end();
            stackID.reserve(stackID.size() + std::distance(first_child, last_child));
            stackID.insert(stackID.end(), first_child, last_child);
        }
    }

    // Copy all IDs (except the one as the parameter of the function) from unordered_set to vector
    visited.erase(id);
    std::vector<AreaID> children_areas;
    children_areas.reserve(visited.size());
    children_areas.insert(children_areas.end(), visited.begin(), visited.end());
    return children_areas;
}

AreaID Datastructures::common_area_of_subareas(AreaID id1, AreaID id2)
{
    std::vector<AreaID> all_parents1 = subarea_in_areas(id1);
    std::vector<AreaID> all_parents2 = subarea_in_areas(id2);
    if (all_parents1.empty() || all_parents2.empty()) {
        return NO_AREA;
    }

    // Highest common ancestor is the last element and start to check from this common ancestor
    std::vector<AreaID>::const_reverse_iterator parent1 = all_parents1.crbegin();
    std::vector<AreaID>::const_reverse_iterator parent2 = all_parents2.crbegin();
    while(*parent1 == *parent2) {
        if (parent1 == all_parents1.crend()-1 || parent2 == all_parents2.crend()-1) {
            break;
        }
        parent1++;
        parent2++;
    }

    // Check if having highest common ancestor or not
    if (parent1 == all_parents1.crbegin()) {
        if(*parent1 != *parent2) {
            return NO_AREA;
        }
    }
    return *parent1;
}

// Phase 2 operations

int Datastructures::way_count()
{
    return WAYS.size();
}

int Datastructures::calculate_distance(const std::vector<Coord>& coords)
{
    int total_length = 0;
    for(auto cur_coord_it = coords.begin(); cur_coord_it != coords.end()-1; cur_coord_it++) {
        total_length += floor(*cur_coord_it - *(cur_coord_it+1));
    }
    return total_length;
}

void Datastructures::add_crossroads(const Coord& coord, const std::pair<WayID, Coord>& way_to)
{
    auto search_crossroad = CROSSROADS.find(coord);
    if(search_crossroad != CROSSROADS.end()) {
        search_crossroad->second->ways_from.push_back(way_to);
    }
    std::vector<std::pair<WayID, Coord>> list_way_to = {way_to};
    Crossroad* new_crossroad = new Crossroad({coord, list_way_to, std::make_pair(nullptr, NO_WAY), 0, ColorType::WHITE});
    CROSSROADS.emplace(coord, new_crossroad);
}

std::vector<WayID> Datastructures::all_ways()
{
    std::vector<WayID> ways;
    ways.reserve(way_count());

    for(auto const& pair_way : WAYS) {
        ways.push_back(pair_way.first);
    }
    return ways;
}

bool Datastructures::add_way(WayID id, std::vector<Coord> coords)
{
    auto [it, insert] = WAYS.emplace(id, new Way({coords, calculate_distance(coords), false}));
    if (insert){
        // Add path going in 2 directions
        add_crossroads(coords.front(), {id, coords.back()});
        add_crossroads(coords.back(), {id, coords.front()});
    }
    return insert;
}

std::vector<std::pair<WayID, Coord>> Datastructures::ways_from(Coord xy)
{
    auto search_crossroad = CROSSROADS.find(xy);
    if(search_crossroad == CROSSROADS.end()) {
        return {};
    }
    return search_crossroad->second->ways_from;
}

std::vector<Coord> Datastructures::get_way_coords(WayID id)
{
    auto search_way = WAYS.find(id);
    if(search_way == WAYS.end()) {
        return {NO_COORD};
    }
    return search_way->second->way_coords;
}

void Datastructures::clear_ways()
{
    for(auto const& way : WAYS){
        delete way.second;
    }
    WAYS.clear();

    for (const auto& crossroad : CROSSROADS) {
        delete crossroad.second;
    }
    CROSSROADS.clear();
}

std::vector<std::tuple<Coord, WayID, Distance> > Datastructures::route_any(Coord fromxy, Coord toxy)
{
    return route_least_crossroads(fromxy, toxy);
}

bool Datastructures::remove_way(WayID id)
{
    auto search_way = WAYS.find(id);
    if(search_way == WAYS.end()) {
        return false;
    }

    // Get coord before delete this pointer, then remove it from the unorder_map
    const Coord coord1 = search_way->second->way_coords.front();
    const Coord coord2 = search_way->second->way_coords.back();
    delete search_way->second;
    WAYS.erase(id);

    // Go to each coord, find it in the unordered_map, then access the way_from attribute, remove the pair <wayid, coord next>
    Crossroad* crossroad1 = CROSSROADS.at(coord1);
    Crossroad* crossroad2 = CROSSROADS.at(coord2);
    crossroad1->ways_from.erase(
        std::find_if(crossroad1->ways_from.begin(), crossroad1->ways_from.end(),
            [&id](std::pair<WayID, Coord> removing_way){
                return removing_way.first == id;
            }
        )
    );
    crossroad2->ways_from.erase(
        std::find_if(crossroad2->ways_from.begin(), crossroad2->ways_from.end(),
            [&id](std::pair<WayID, Coord> removing_way){
                return removing_way.first == id;
            }
        )
    );

    // If way list is empty, remove crossroad by delete crossroad pointer, then erase crossroad from unordered_map
    if (crossroad1->ways_from.empty()) {
        delete crossroad1;
        CROSSROADS.erase(coord1);
    }
    if (crossroad2->ways_from.empty()) {
        delete crossroad2;
        CROSSROADS.erase(coord2);
    }
    return true;
}

void Datastructures::construct_path(Crossroad* cur_crossroad, std::vector<std::tuple<Coord, WayID, Distance>>& path)
{
    Crossroad* prev_crossroad = cur_crossroad->way_to.first;
    Coord prev_coord = prev_crossroad->coord;
    Distance distance = prev_crossroad->distance_from_start;
    path.push_back(std::make_tuple(prev_coord, cur_crossroad->way_to.second, distance));
    if (distance != 0){
        construct_path(prev_crossroad, path);
    }
}

std::vector<std::tuple<Coord, WayID, Distance> > Datastructures::route_least_crossroads(Coord fromxy, Coord toxy)
{
    auto start_crossroad_pair = CROSSROADS.find(fromxy);
    auto end_crossroad_pair = CROSSROADS.find(toxy);
    if(start_crossroad_pair == CROSSROADS.end() || end_crossroad_pair == CROSSROADS.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }

    // BFS
    for(const auto& crossroad : CROSSROADS) {
        crossroad.second->color = ColorType::WHITE;
    }

    Crossroad* start_crossroad = start_crossroad_pair->second;
    start_crossroad->way_to = std::make_pair(nullptr, NO_WAY);
    start_crossroad->distance_from_start = 0;
    start_crossroad->color = ColorType::GRAY;

    std::queue<Crossroad*> frontier;
    frontier.push(start_crossroad);

    while (!frontier.empty()) {
        Crossroad* visiting_crossroad = frontier.front();
        frontier.pop();

        for(auto const& way_pair : visiting_crossroad->ways_from) {
            const WayID way_id = way_pair.first;
            Crossroad* next_crossroad = CROSSROADS.at(way_pair.second);

            // Avoid add nodes already in the queue
            if(next_crossroad->color == ColorType::WHITE) {
                next_crossroad->way_to = std::make_pair(visiting_crossroad, way_id);
                next_crossroad->distance_from_start = visiting_crossroad->distance_from_start + WAYS.at(way_id)->length;

                // Stop if find the destination
                if(next_crossroad->coord == toxy) {
                    std::vector<std::tuple<Coord, WayID, Distance>> path = {std::make_tuple(toxy, NO_WAY, next_crossroad->distance_from_start)};
                    construct_path(next_crossroad, path);
                    std::reverse(path.begin(), path.end());
                    return path;
                }

                next_crossroad->color = ColorType::GRAY;
                frontier.push(next_crossroad);
            }
        }
    }
    return {};
}

void Datastructures::construct_cycle(Crossroad* cur_crossroad, std::vector<std::tuple<Coord, WayID> > &path)
{
    Crossroad* prev_crossroad = cur_crossroad->way_to.first;
    if (prev_crossroad != nullptr){
        Coord prev_coord = prev_crossroad->coord;
        path.push_back(std::make_tuple(prev_coord, cur_crossroad->way_to.second));
        construct_cycle(prev_crossroad, path);
    }
}

std::vector<std::tuple<Coord, WayID> > Datastructures::route_with_cycle(Coord fromxy)
{
    auto start_crossroad_pair = CROSSROADS.find(fromxy);
    if(start_crossroad_pair == CROSSROADS.end()) {
        return {{NO_COORD, NO_WAY}};
    }

    // DFS
    for(const auto& crossroad : CROSSROADS) {
        crossroad.second->color = ColorType::WHITE;
    }

    Crossroad* start_crossroad = start_crossroad_pair->second;
    start_crossroad->way_to = std::make_pair(nullptr, NO_WAY);

    std::stack<Crossroad*> frontier;
    frontier.push(start_crossroad);

    while (!frontier.empty()) {
        Crossroad* visiting_crossroad = frontier.top();
        frontier.pop();

        if(visiting_crossroad->color == ColorType::WHITE) {
            for(auto const& way_pair : visiting_crossroad->ways_from) {
                const WayID way_id = way_pair.first;

                // Not proceed if the neighbor is immediately the one it came from earlier
                if(way_id == visiting_crossroad->way_to.second) {
                    continue;
                }

                Crossroad* next_crossroad = CROSSROADS.at(way_pair.second);
                if(next_crossroad->color == ColorType::WHITE) {
                    next_crossroad->way_to = std::make_pair(visiting_crossroad, way_id);
                    frontier.push(next_crossroad);
                }

                // Detect cycle
                else if(next_crossroad->color == ColorType::GRAY) {
                    std::vector<std::tuple<Coord, WayID>> cycle = {std::make_tuple(next_crossroad->coord, NO_WAY),
                                                                   std::make_tuple(visiting_crossroad->coord, way_id)};
                    construct_cycle(visiting_crossroad, cycle);
                    std::reverse(cycle.begin(), cycle.end());
                    return cycle;
                }
            }

            visiting_crossroad->color = ColorType::GRAY;
        }
    }

    return {};
}

std::vector<std::tuple<Coord, WayID, Distance> > Datastructures::route_shortest_distance(Coord fromxy, Coord toxy)
{
    auto start_crossroad_pair = CROSSROADS.find(fromxy);
    auto end_crossroad_pair = CROSSROADS.find(toxy);
    if(start_crossroad_pair == CROSSROADS.end() || end_crossroad_pair == CROSSROADS.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }

    // Dijkstra's algorithm
    for(const auto& crossroad : CROSSROADS) {
        crossroad.second->color = ColorType::WHITE;
        crossroad.second->way_to = std::make_pair(nullptr, NO_WAY);
        crossroad.second->distance_from_start = MAX_VALUE;
    }

    Crossroad* start_crossroad = start_crossroad_pair->second;
    start_crossroad->distance_from_start = 0;
    start_crossroad->color = ColorType::GRAY;

    // Defined comparator used in priority queue (smallest distance crossroad comes top)
    struct compare_crossroad_pointer{
        bool operator()(const Crossroad* cross_ptr1, const Crossroad* cross_ptr2){
            return cross_ptr1->distance_from_start > cross_ptr2->distance_from_start;
        }
    };
    std::priority_queue<Crossroad*, std::vector<Crossroad*>, compare_crossroad_pointer> frontier;
    frontier.push(start_crossroad);

    while (!frontier.empty()) {
        Crossroad* visiting_crossroad = frontier.top();
        frontier.pop();

        // Skip processing nodes that have been fully processed since there's no better way to this node
        if (visiting_crossroad->color == ColorType::BLACK) {
            continue;
        }

        // Stop when finally we have to process the destination node, i.e there's no better route to this node
        if(visiting_crossroad->coord == toxy) {
            std::vector<std::tuple<Coord, WayID, Distance>> path = {std::make_tuple(toxy, NO_WAY, visiting_crossroad->distance_from_start)};
            construct_path(visiting_crossroad, path);
            std::reverse(path.begin(), path.end());
            return path;
        }

        for(auto const& way_pair : visiting_crossroad->ways_from) {
            const WayID way_id = way_pair.first;
            Crossroad* next_crossroad = CROSSROADS.at(way_pair.second);

            // Skip node (explained above)
            if (visiting_crossroad->color == ColorType::BLACK ) {
                continue;
            }

            if(next_crossroad->color == ColorType::WHITE) {
                next_crossroad->color = ColorType::GRAY;
                frontier.push(next_crossroad);
            }

            // Relaxation Dijkstra
            int new_distance = visiting_crossroad->distance_from_start + WAYS.at(way_id)->length;
            if (next_crossroad->distance_from_start > new_distance) {
                next_crossroad->distance_from_start = new_distance;
                next_crossroad->way_to = std::make_pair(visiting_crossroad, way_id);

                // Since priority queue not support decrease-key, we have to manually push it to the priority queue
                frontier.push(next_crossroad);
            }
        }

        visiting_crossroad->color = ColorType::BLACK;
    }

    return {};
}

Distance Datastructures::trim_ways()
{
    Distance total_distance = 0;

    struct compare_way_pointer{
        bool operator()(const Way* way_ptr1, const Way* way_ptr2){
            return way_ptr1->length > way_ptr2->length;
        }
    };
    std::priority_queue<Way*, std::vector<Way*>, compare_way_pointer> all_ways;

    // Prim's algorithm
    for(const auto& crossroad : CROSSROADS) {
        crossroad.second->color = ColorType::WHITE;
    }

    // This for loop just to make sure we can process all nodes even when the graph is disconnected
    for(auto const& crossroad_pair : CROSSROADS){
        // Firstly, add any nodes to the queue
        if (crossroad_pair.second->color == ColorType::WHITE) {
            crossroad_pair.second->color = ColorType::GRAY;
            for ( const auto& pair_way : crossroad_pair.second->ways_from) {
                all_ways.push(WAYS.at(pair_way.first));
            }
        }

        while (!all_ways.empty()) {
            Way* way = all_ways.top();
            all_ways.pop();
            Crossroad* cross1 = CROSSROADS.at(way->way_coords.front());
            Crossroad* cross2 = CROSSROADS.at(way->way_coords.back());

            // Avoid cycle, by skipping edge that both of its ends have already visited
            if (cross1->color == ColorType::GRAY && cross2->color == ColorType::GRAY) {
                continue;
            }

            // Add edge to min spanning tree
            way->mst = true;
            total_distance += way->length;

            // Add all next edges of the current nodes to queue
            // Get the correct end of the edge
            Crossroad* visiting_cross = (cross1->color == ColorType::WHITE) ? cross1 : cross2;
            visiting_cross->color = ColorType::GRAY;

            for(auto const& next_way_pair : visiting_cross->ways_from) {
                // Only add edge that leading to the unvisited nodes
                if (CROSSROADS.at(next_way_pair.second)->color == ColorType::WHITE) {
                    all_ways.push(WAYS.at(next_way_pair.first));
                }
            }
        }
    }

    // Remove all ways that are not in MST
    for(const auto& way_pair : WAYS) {
        if(!way_pair.second->mst) {
            remove_way(way_pair.first);
        }
    }

    return total_distance;
}
