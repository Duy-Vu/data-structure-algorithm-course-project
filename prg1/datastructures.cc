// Datastructures.cc

#include "datastructures.hh"
#include <random>
#include <unordered_set>
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
}
Datastructures::~Datastructures()
{
    clear_all();
}

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
        qDebug() << "no axiest"<< id;
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
        AreaID child_id = stackID.back();   // O(1)
        stackID.pop_back();     // O(1)
        if (visited.find(child_id) == visited.end()) {  // O(1)
            // Label as visited
            visited.emplace(child_id);  // O(1)
            // Add to stack
            auto first_child = AREAS.at(child_id).childrenID.begin();
            auto last_child = AREAS.at(child_id).childrenID.end();
            stackID.reserve(stackID.size() + std::distance(first_child, last_child));   // O(n)
            stackID.insert(stackID.end(), first_child, last_child);     // O(n)
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
