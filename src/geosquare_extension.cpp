#define DUCKDB_EXTENSION_MAIN

#include "geosquare_extension.hpp"
#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension_util.hpp"
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>
#include "duckdb/function/scalar/string_functions.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include <geos/geom/Geometry.h>
#include <geos/io/WKTReader.h>
#include <geos/geom/Envelope.h>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <cmath>


namespace duckdb {

static const std::vector<std::vector<char>> CODE_ALPHABET = {
        {'2', '3', '4', '5', '6'},
        {'7', '8', '9', 'C', 'E'},
        {'F', 'G', 'H', 'J', 'L'},
        {'M', 'N', 'P', 'Q', 'R'},
        {'T', 'V', 'W', 'X', 'Y'}
    };
static const  std::unordered_map<char, std::pair<int, int>> CODE_ALPHABET_VALUE = {
        {'2', {0, 0}}, {'3', {0, 1}}, {'4', {0, 2}}, {'5', {0, 3}}, {'6', {0, 4}},
        {'7', {1, 0}}, {'8', {1, 1}}, {'9', {1, 2}}, {'C', {1, 3}}, {'E', {1, 4}},
        {'F', {2, 0}}, {'G', {2, 1}}, {'H', {2, 2}}, {'J', {2, 3}}, {'L', {2, 4}},
        {'M', {3, 0}}, {'N', {3, 1}}, {'P', {3, 2}}, {'Q', {3, 3}}, {'R', {3, 4}},
        {'T', {4, 0}}, {'V', {4, 1}}, {'W', {4, 2}}, {'X', {4, 3}}, {'Y', {4, 4}}
    };
static const std::unordered_map<int, std::vector<char>> CODE_ALPHABET_ = {
        {5, {'2', '3', '4', '5', '6', '7', '8', '9', 'C', 'E', 'F', 'G', 'H', 'J', 'L', 'M', 'N', 'P', 'Q', 'R', 'T', 'V', 'W', 'X', 'Y'}},
        {2, {'2', '3', '7', '8'}}
    };
static const std::vector<int> D = {5, 2, 5, 2, 5, 2, 5, 2, 5, 2, 5, 2, 5, 2, 5};
static const std::unordered_map<int, int> size_level = {
        {10000000, 1},
        {5000000, 2},
        {1000000, 3},
        {500000, 4},
        {100000, 5},
        {50000, 6},
        {10000, 7},
        {5000, 8},
        {1000, 9},
        {500, 10},
        {100, 11},
        {50, 12},
        {10, 13},
        {5, 14},
        {1, 15}
    };

std::string lonlatTOgid(double &longitude, double &latitude, int &level) {
    // Validate input parameters
    if (longitude < -180 || longitude > 180) {
        throw std::invalid_argument("Longitude must be between -180 and 180");
    }
    if (latitude < -90 || latitude > 90) {
        throw std::invalid_argument("Latitude must be between -90 and 90");
    }
    if (level < 1 || level > 15) {
        throw std::invalid_argument("Level must be between 1 and 15");
    }

    // Initialize gid with the specified level
    std::string gid(level, '0');
    double LAT_RANGED_MIN = LAT_RANGED_MIN_INIT;
    double LAT_RANGED_MAX = LAT_RANGED_MAX_INIT;
    double LON_RANGED_MIN = LON_RANGED_MIN_INIT;
    double LON_RANGED_MAX = LON_RANGED_MAX_INIT;

    // Calculate the gid based on the longitude, latitude, and level
    for (int i = 0; i < level; ++i) {
        double part_x = (LON_RANGED_MAX - LON_RANGED_MIN) / D[i];
        double part_y = (LAT_RANGED_MAX - LAT_RANGED_MIN) / D[i];

        int position_x = static_cast<int>((longitude - LON_RANGED_MIN) / part_x);
        int position_y = static_cast<int>((latitude - LAT_RANGED_MIN) / part_y);

        gid[i] = CODE_ALPHABET[position_y][position_x];

        LON_RANGED_MIN += part_x * position_x;
        LON_RANGED_MAX = LON_RANGED_MIN + part_x;
        LAT_RANGED_MIN += part_y * position_y;
        LAT_RANGED_MAX = LAT_RANGED_MIN + part_y;
    }

    return gid;
}

void gidTOlonlat(const string_t &gid, LatLng &lonlat) {
    // Validate input parameters
    if (gid.GetSize() == 0) {
        throw std::invalid_argument("GID cannot be empty");
    }

    double LAT_RANGED_MIN = LAT_RANGED_MIN_INIT;
    double LAT_RANGED_MAX = LAT_RANGED_MAX_INIT;
    double LON_RANGED_MIN = LON_RANGED_MIN_INIT;
    double LON_RANGED_MAX = LON_RANGED_MAX_INIT;

    for (size_t i = 0; i < gid.GetSize(); ++i) {
        char gid_char = gid.GetData()[i];
        if (CODE_ALPHABET_VALUE.find(gid_char) == CODE_ALPHABET_VALUE.end()) {
            throw std::invalid_argument("Invalid character in GID: " + std::string(1, gid_char));
        }

        double part_x = (LON_RANGED_MAX - LON_RANGED_MIN) / D[i];
        double part_y = (LAT_RANGED_MAX - LAT_RANGED_MIN) / D[i];

        double position_x = CODE_ALPHABET_VALUE.at(gid_char).second;
        double position_y = CODE_ALPHABET_VALUE.at(gid_char).first;

        LON_RANGED_MIN += part_x * position_x;
        LON_RANGED_MAX = LON_RANGED_MIN + part_x;
        LAT_RANGED_MIN += part_y * position_y;
        LAT_RANGED_MAX = LAT_RANGED_MIN + part_y;
    }

    lonlat.lng = LON_RANGED_MIN + (LON_RANGED_MAX - LON_RANGED_MIN) / 2;
    lonlat.lat = LAT_RANGED_MIN + (LAT_RANGED_MAX - LAT_RANGED_MIN) / 2;
}

void gidToBound(const string_t &gid, BBox &bound) {
    // Validate input parameters
    if (gid.GetSize() == 0) {
        throw std::invalid_argument("GID cannot be empty");
    }

    double LAT_RANGED_MIN = LAT_RANGED_MIN_INIT;
    double LAT_RANGED_MAX = LAT_RANGED_MAX_INIT;
    double LON_RANGED_MIN = LON_RANGED_MIN_INIT;
    double LON_RANGED_MAX = LON_RANGED_MAX_INIT;

    for (size_t i = 0; i < gid.GetSize(); ++i) {
        char gid_char = gid.GetData()[i];
        if (CODE_ALPHABET_VALUE.find(gid_char) == CODE_ALPHABET_VALUE.end()) {
            throw std::invalid_argument("Invalid character in GID: " + std::string(1, gid_char));
        }

        double part_x = (LON_RANGED_MAX - LON_RANGED_MIN) / D[i];
        double part_y = (LAT_RANGED_MAX - LAT_RANGED_MIN) / D[i];

        double position_x = CODE_ALPHABET_VALUE.at(gid_char).second;
        double position_y = CODE_ALPHABET_VALUE.at(gid_char).first;

        LON_RANGED_MIN += part_x * position_x;
        LON_RANGED_MAX = LON_RANGED_MIN + part_x;
        LAT_RANGED_MIN += part_y * position_y;
        LAT_RANGED_MAX = LAT_RANGED_MIN + part_y;
    }

    bound.west = LON_RANGED_MIN;
    bound.east = LON_RANGED_MAX;
    bound.south = LAT_RANGED_MIN;
    bound.north = LAT_RANGED_MAX;
}

void gidToBoundWkt(const string_t &gid, std::string &wkt) {
    // Validate input parameters
    if (gid.GetSize() == 0) {
        throw std::invalid_argument("GID cannot be empty");
    }

    // Get the bounding box for the GID
    BBox bound;
    gidToBound(gid, bound);

    // Create the WKT string
    std::ostringstream wkt_temp;
    wkt_temp.precision(16);
    wkt_temp << "POLYGON(("
             << bound.west << " " << bound.south << ", "
             << bound.east << " " << bound.south << ", "
             << bound.east << " " << bound.north << ", "
             << bound.west << " " << bound.north << ", "
             << bound.west << " " << bound.south << "))";

    // Assign the WKT string to the output parameter
    wkt = wkt_temp.str();
}

// Add overload for std::string parameter
void gidToBoundWkt(const std::string &gid, std::string &wkt) {
    // Validate input parameters
    if (gid.empty()) {
        throw std::invalid_argument("GID cannot be empty");
    }

    // Convert std::string to string_t for the existing function
    string_t gid_str(gid.c_str(), gid.size());
    gidToBoundWkt(gid_str, wkt);
}

void _area_ratio(string_t &a, string_t &b, double &area_ratio) {
    try {
        // Convert GIDs to geometries
        std::string geometry_a;
        gidToBoundWkt(a, geometry_a);
        std::string geometry_b;
        gidToBoundWkt(b, geometry_b);

        // Load geometries as GEOS objects
        geos::io::WKTReader reader;
        std::unique_ptr<geos::geom::Geometry> geom_a_ptr(reader.read(geometry_a));
        std::unique_ptr<geos::geom::Geometry> geom_b_ptr(reader.read(geometry_b));

        if (!geom_a_ptr || !geom_b_ptr) {
            throw std::runtime_error("Failed to read geometries from WKT");
        }

        // Calculate intersection area
        std::unique_ptr<geos::geom::Geometry> intersection_geom(geom_a_ptr->intersection(geom_b_ptr.get()));
        double intersection_area = intersection_geom->getArea();

        double area_a = geom_a_ptr->getArea();

        // Calculate area ratio
        if (area_a == 0) {
            throw std::runtime_error("Area of geometry A is zero, cannot calculate area ratio");
        }
        area_ratio = intersection_area / area_a;
    } catch (const std::exception &e) {
        throw std::runtime_error(std::string("Error in _area_ratio: ") + e.what());
    }
}

int _get_resolution(string_t &key) {
    return key.GetSize();
}

void _to_children(const string_t &key, std::vector<std::string> &children) {
    // Validate input parameters
    if (key.GetSize() == 0) {
        throw std::invalid_argument("Key cannot be empty");
    }

    size_t key_size = key.GetSize();
    if (key_size >= D.size()) {
        throw std::out_of_range("Key size is out of range");
    }

    // Get the alphabet for the current key size
    const std::vector<char>& alphabet = CODE_ALPHABET_.at(D[key_size]);

    // Generate children keys
    for (const char& c : alphabet) {
        children.push_back(key.GetString() + c);
    }
}

// Add overload for std::string parameter
void _to_children(const std::string &key, std::vector<std::string> &children) {
    // Validate input parameters
    if (key.empty()) {
        throw std::invalid_argument("Key cannot be empty");
    }

    size_t key_size = key.size();
    if (key_size >= D.size()) {
        throw std::out_of_range("Key size is out of range");
    }

    // Get the alphabet for the current key size
    const std::vector<char>& alphabet = CODE_ALPHABET_.at(D[key_size]);

    // Generate children keys
    for (const char& c : alphabet) {
        children.push_back(key + c);
    }
}

string_t _to_parent(string_t &key) {
    if (key.GetSize() > 1) {
        std::string parent_str(key.GetData(), key.GetSize() - 1);
        return string_t(parent_str);
    } else {
        return key;
    }
}

int count_max_gids(std::string &geometry, int &level) {
    try {
        // Read the geometry from the WKT string
        geos::io::WKTReader reader;
        std::unique_ptr<geos::geom::Geometry> geom_ptr(reader.read(geometry));

        // Ensure the geometry is valid
        if (!geom_ptr) {
            throw std::runtime_error("Failed to read geometry from WKT");
        }

        // Get the envelope of the geometry
        const geos::geom::Envelope* envelope = geom_ptr->getEnvelopeInternal();
        double lat_min = envelope->getMinY();
        double lat_max = envelope->getMaxY();
        double lon_min = envelope->getMinX();
        double lon_max = envelope->getMaxX();

        // Calculate the grid length
        double gridLength = LON_RANGED_MAX_INIT - LON_RANGED_MIN_INIT;
        for (int i = 0; i < level; ++i) {
            gridLength /= D[i];
        }

        // Calculate the number of grid cells in the x and y directions
        int x_count = std::ceil((lon_max - lon_min) / gridLength);
        int y_count = std::ceil((lat_max - lat_min) / gridLength);

        // Return the total number of grid cells
        return x_count * y_count;
    } catch (const std::exception &e) {
        throw std::runtime_error(std::string("Error in count_max_gids: ") + e.what());
    }
}

void get_intersect_key(std::string &geometry, const std::string &initial_key, int &resolution, bool parrent, int &max_gids, std::vector<std::string> &contained_keys) {
    // std::function<void(const std::string&, bool)> func = [&](std::string key, bool approved) {
    std::vector<std::vector<std::string>> queue;
    queue.push_back({initial_key, "false"});
    while (!queue.empty()) {
        std::vector<std::string> current = queue.back();
        queue.pop_back();
        std::string key = current[0];
        bool approved = current[1] == "true";
        if (approved) {
            if (key.length() == resolution) {
                contained_keys.push_back(key);
            } else {
                std::vector<std::string> childrens;
                _to_children(key, childrens);
                for (const std::string& child_key : childrens) {
                    // func(std::string(child_key.c_str(), child_key.size()), true);
                    queue.push_back({std::string(child_key.c_str(), child_key.size()), "true"});
                }
            }
        } else {
            std::string geometry_key;
            gidToBoundWkt(key, geometry_key);
            std::unique_ptr<geos::geom::Geometry> geom_key_ptr;
            std::unique_ptr<geos::geom::Geometry> geom_ptr;
            std::unique_ptr<geos::geom::Geometry> intersection_geom_ptr;

            try {
                geom_key_ptr = std::unique_ptr<geos::geom::Geometry>(geos::io::WKTReader().read(geometry_key));
                geom_ptr = std::unique_ptr<geos::geom::Geometry>(geos::io::WKTReader().read(geometry));
                intersection_geom_ptr = std::unique_ptr<geos::geom::Geometry>(geom_ptr->intersection(geom_key_ptr.get()));
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error reading geometries: ") + e.what());
            }
            double area_ratio = intersection_geom_ptr->getArea() / geom_key_ptr->getArea();
            if (area_ratio == 0) {
                if (key.length() == 1) {
                    const auto& last_idx = std::find(CODE_ALPHABET_.at(D[0]).begin(), CODE_ALPHABET_.at(D[0]).end(), key[key.length() - 1]) - CODE_ALPHABET_.at(D[0]).begin();
                    if (last_idx < 24 && key.length() == 1) {
                        std::string new_key(&CODE_ALPHABET_.at(D[0]).at(last_idx + 1), 1);
                        // func(new_key, false);
                        queue.push_back({new_key, "false"});
                    }
                }
            } else if (area_ratio == 1) {
                // func(key, true);
                queue.push_back({key, "true"});
            } else if (key.length() == resolution && parrent) {
                contained_keys.push_back(key);
            } else if (key.length() == resolution && area_ratio > 0.5 && !parrent) {
                contained_keys.push_back(key);
            } else if (key.length() == resolution) {
                // Do nothing
            } else {
                std::vector<std::string> childrens;
                _to_children(key, childrens);
                for (const std::string& child_key : childrens) {
                    // func(std::string(child_key.c_str(), child_key.size()), false);
                    queue.push_back({std::string(child_key.c_str(), child_key.size()), "false"});
                }
            }
        }
    };

    // func(initial_key, false);
}

void count_intersect_gids(std::string &geometry, const std::string &initial_key, int &resolution, bool parrent, int &counts) {
    std::function<void(const std::string&, bool)> func = [&](std::string key, bool approved) {
        if (approved) {
            if (key.length() == resolution) {
                counts++;
            } else {
                std::vector<std::string> childrens;
                _to_children(key, childrens);
                for (const std::string& child_key : childrens) {
                    func(std::string(child_key.c_str(), child_key.size()), true);
                }
            }
        } else {
            std::string geometry_key;
            gidToBoundWkt(key, geometry_key);
            std::unique_ptr<geos::geom::Geometry> geom_key_ptr;
            std::unique_ptr<geos::geom::Geometry> geom_ptr;
            std::unique_ptr<geos::geom::Geometry> intersection_geom_ptr;

            try {
                geom_key_ptr = std::unique_ptr<geos::geom::Geometry>(geos::io::WKTReader().read(geometry_key));
                geom_ptr = std::unique_ptr<geos::geom::Geometry>(geos::io::WKTReader().read(geometry));
                intersection_geom_ptr = std::unique_ptr<geos::geom::Geometry>(geom_ptr->intersection(geom_key_ptr.get()));
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error reading geometries: ") + e.what());
            }
            double area_ratio = intersection_geom_ptr->getArea() / geom_key_ptr->getArea();
            if (area_ratio == 0) {
                if (key.length() == 1) {
                    const auto& last_idx = std::find(CODE_ALPHABET_.at(D[0]).begin(), CODE_ALPHABET_.at(D[0]).end(), key[key.length() - 1]) - CODE_ALPHABET_.at(D[0]).begin();
                    if (last_idx < 24 && key.length() == 1) {
                        std::string new_key(&CODE_ALPHABET_.at(D[0]).at(last_idx + 1), 1);
                        func(new_key, false);
                    }
                }
            } else if (area_ratio == 1) {
                func(key, true);
            } else if (key.length() == resolution && parrent) {
                counts++;
            } else if (key.length() == resolution && area_ratio > 0.5 && !parrent) {
                counts++;
            } else if (key.length() == resolution) {
                // Do nothing
            } else {
                std::vector<std::string> childrens;
                _to_children(key, childrens);
                for (const std::string& child_key : childrens) {
                    func(std::string(child_key.c_str(), child_key.size()), false);
                }
            }
        }
    };

    func(initial_key, false);
}


inline void GeosquaregidTOlonlat(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    UnaryExecutor::Execute<string_t, list_entry_t>(inputs, result, args.size(),
        [&](string_t gid) {
            list_entry_t result_entry;
            LatLng lonlat;

            try {
                // Convert GID to latitude and longitude
                gidTOlonlat(gid, lonlat);
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error in gidTOlonlat: ") + e.what());
            }

            // Push the longitude and latitude to the result vector
            ListVector::PushBack(result, lonlat.lng);
            ListVector::PushBack(result, lonlat.lat);

            // Set the length of the result entry
            result_entry.length = 2;

            return result_entry;
        });
}

inline void GeosquaregidTOBound(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    UnaryExecutor::Execute<string_t, list_entry_t>(inputs, result, args.size(),
        [&](string_t gid) {
            list_entry_t result_entry;
            BBox bound;

            try {
                // Convert GID to bounding box
                gidToBound(gid, bound);
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error in gidToBound: ") + e.what());
            }

            // Push the bounding box coordinates to the result vector
            ListVector::PushBack(result, bound.west);
            ListVector::PushBack(result, bound.east);
            ListVector::PushBack(result, bound.south);
            ListVector::PushBack(result, bound.north);

            // Set the length of the result entry
            result_entry.length = 4;

            return result_entry;
        });
}

inline void GeosquaregidTOPointWkt(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    UnaryExecutor::Execute<string_t, string_t>(inputs, result, args.size(),
        [&](string_t gid) {
            LatLng lonlat;

            try {
                // Convert GID to latitude and longitude
                gidTOlonlat(gid, lonlat);
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error in gidTOlonlat: ") + e.what());
            }

            // Create WKT string for the point
            std::ostringstream wkt;
            wkt.precision(16);
            wkt << "POINT(" << lonlat.lng << " " << lonlat.lat << ")";

            // Return the WKT string
            return StringVector::AddString(result, wkt.str());
        });
}

inline void GeosquarelonlatTOgid(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    auto &inputs2 = args.data[1];
    auto &inputs3 = args.data[2];
    TernaryExecutor::Execute<double, double, int, string_t>(inputs, inputs2, inputs3, result, args.size(),
        [&](double longitude, double latitude, int level) {
            return StringVector::AddString(result, lonlatTOgid(longitude, latitude, level));
        });

}

inline void GeosquaregidTOBoundWkt(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    UnaryExecutor::Execute<string_t, string_t>(inputs, result, args.size(),
        [&](string_t gid) {
            BBox bound;

            try {
                // Convert GID to bounding box
                gidToBound(gid, bound);
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error in gidToBound: ") + e.what());
            }

            // Create WKT string for the bounding box
            std::ostringstream wkt;
            wkt.precision(16);
            wkt << "POLYGON((" << bound.west << " " << bound.south << ", "
                << bound.east << " " << bound.south << ", "
                << bound.east << " " << bound.north << ", "
                << bound.west << " " << bound.north << ", "
                << bound.west << " " << bound.south << "))";

            // Return the WKT string
            return StringVector::AddString(result, wkt.str());
        });
}

inline void Geosquarearea_ratio(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    auto &inputs2 = args.data[1];
    BinaryExecutor::Execute<string_t, string_t, double>(inputs, inputs2, result, args.size(),
        [&](string_t a, string_t b) {
            std::string geom_a;
            std::string geom_b;
            double area_ratio_temp;

            try {
                // Convert GIDs to WKT geometries
                gidToBoundWkt(a, geom_a);
                gidToBoundWkt(b, geom_b);

                // Calculate area ratio
                _area_ratio(a, b, area_ratio_temp);
            } catch (const std::exception &e) {
                throw std::runtime_error(std::string("Error in Geosquarearea_ratio: ") + e.what());
            }

            return area_ratio_temp;
        });
}

inline void Geosquaregid_children(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    UnaryExecutor::Execute<string_t, list_entry_t>(inputs, result, args.size(),
        [&](string_t key) {
            list_entry_t result_entry;
            std::vector<std::string> children;
            _to_children(key, children);
            for (auto &child : children) {
                ListVector::PushBack(result, child);
            }
            auto result_data = FlatVector::GetData<list_entry_t>(result);
            result_entry.length = children.size();
            return result_entry;
        });
}

inline void Geosquaregid_parent(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    UnaryExecutor::Execute<string_t, string_t>(inputs, result, args.size(),
        [&](string_t key) {
            return StringVector::AddString(result, _to_parent(key));
        });
}

inline void Geosquarecount_gids(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    auto &inputs2 = args.data[1];
    BinaryExecutor::Execute<string_t, int, int>(inputs, inputs2, result, args.size(),
        [&](string_t geometry, int resolution) {
            // Convert geometry to std::string
            std::string geometry_str(geometry.GetData(), geometry.GetSize());

            // Get the level from the resolution
            int level;
            try {
                level = size_level.at(resolution);
            } catch (const std::out_of_range &e) {
                throw std::runtime_error("Invalid resolution: " + std::to_string(resolution));
            }

            // Define the initial key
            std::string initial_key(1, '2');

            // Get the contained keys using polyfill
            int counts = 0;
            try {
                count_intersect_gids(geometry_str, initial_key, level, false, counts);
            } catch (const std::exception &e) {
                throw std::runtime_error("Error in polyfill: " + std::string(e.what()));
            }
            return counts;

        });
}

inline void Geosquarepolyfill(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    auto &inputs2 = args.data[1];
    BinaryExecutor::Execute<string_t, int, list_entry_t>(inputs, inputs2, result, args.size(),
        [&](string_t geometry, int resolution) {
            // Convert geometry to std::string
            std::string geometry_str(geometry.GetData(), geometry.GetSize());

            // Get the level from the resolution
            int level;
            try {
                level = size_level.at(resolution);
            } catch (const std::out_of_range &e) {
                throw std::runtime_error("Invalid resolution: " + std::to_string(resolution));
            }

            // Count the maximum number of GIDs
            int max_gids = count_max_gids(geometry_str, level);
            if (max_gids > 10000000) {
                throw std::runtime_error("The area is too big or the size is too small");
            }

            // Define the initial key
            std::string initial_key(1, '2');

            // Get the contained keys using polyfill
            std::vector<std::string> contained_keys;
            contained_keys.reserve(max_gids);
            try {
                get_intersect_key(geometry_str, initial_key, level, false, max_gids, contained_keys);
            } catch (const std::exception &e) {
                throw std::runtime_error("Error in polyfill: " + std::string(e.what()));
            }

            // Prepare the result entry
            list_entry_t result_entry;
            for (auto &key : contained_keys) {
                ListVector::PushBack(result, key);
            }

            // Set the length of the result entry
            result_entry.length = contained_keys.size();

            // Return the result entry
            return result_entry;
        });
}

inline void GeosquarepolyfillFull(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &inputs = args.data[0];
    auto &inputs2 = args.data[1];
    BinaryExecutor::Execute<string_t, int, list_entry_t>(inputs, inputs2, result, args.size(),
        [&](string_t geometry, int resolution) {
            // Convert geometry to std::string
            std::string geometry_str(geometry.GetData(), geometry.GetSize());

            // Get the level from the resolution
            int level;
            try {
                level = size_level.at(resolution);
            } catch (const std::out_of_range &e) {
                throw std::runtime_error("Invalid resolution: " + std::to_string(resolution));
            }

            // Count the maximum number of GIDs
            int max_gids = count_max_gids(geometry_str, level);
            if (max_gids > 10000000) {
                throw std::runtime_error("The area is too big or the size is too small");
            }

            // Define the initial key
            std::string initial_key(1, '2');

            // Get the contained keys using polyfill
            std::vector<std::string> contained_keys;
            contained_keys.reserve(max_gids);
            try {
                get_intersect_key(geometry_str, initial_key, level, true, max_gids, contained_keys);
            } catch (const std::exception &e) {
                throw std::runtime_error("Error in polyfill: " + std::string(e.what()));
            }

            // Prepare the result entry
            list_entry_t result_entry;
            for (auto &key : contained_keys) {
                ListVector::PushBack(result, key);
            }

            // Set the length of the result entry
            result_entry.length = contained_keys.size();

            // Return the result entry
            return result_entry;
        });
}

static void LoadInternal(DatabaseInstance &instance) {
    auto geosquare_lonlatTOgid_function = ScalarFunction("geosquare_lonlat_to_gid", {LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::INTEGER}, LogicalType::VARCHAR, GeosquarelonlatTOgid);
    ExtensionUtil::RegisterFunction(instance, geosquare_lonlatTOgid_function);

    auto geosquare_gidTOlonlat_function = ScalarFunction("geosquare_gid_to_lonlat", {LogicalType::VARCHAR}, LogicalType::LIST(LogicalType::DOUBLE), GeosquaregidTOlonlat);
    ExtensionUtil::RegisterFunction(instance, geosquare_gidTOlonlat_function);   

    auto geosquare_gidTOPointWkt_function = ScalarFunction("geosquare_gid_to_point_wkt", {LogicalType::VARCHAR}, LogicalType::VARCHAR, GeosquaregidTOPointWkt);
    ExtensionUtil::RegisterFunction(instance, geosquare_gidTOPointWkt_function);

    auto geosquare_gidTOBound_function = ScalarFunction("geosquare_gid_to_bound", {LogicalType::VARCHAR}, LogicalType::LIST(LogicalType::DOUBLE), GeosquaregidTOBound);
    ExtensionUtil::RegisterFunction(instance, geosquare_gidTOBound_function);

    auto geosquare_gidTOBoundWkt_function = ScalarFunction("geosquare_gid_to_bound_wkt", {LogicalType::VARCHAR}, LogicalType::VARCHAR, GeosquaregidTOBoundWkt);
    ExtensionUtil::RegisterFunction(instance, geosquare_gidTOBoundWkt_function);

    auto geosquare_to_children = ScalarFunction("geosquare_gid_children", {LogicalType::VARCHAR}, LogicalType::LIST(LogicalType::VARCHAR), Geosquaregid_children);
    ExtensionUtil::RegisterFunction(instance, geosquare_to_children);

    auto geosquare_to_parrent = ScalarFunction("geosquare_gid_parent", {LogicalType::VARCHAR}, LogicalType::VARCHAR, Geosquaregid_parent);
    ExtensionUtil::RegisterFunction(instance, geosquare_to_parrent);

    auto geosquare_polyfill = ScalarFunction("geosquare_polyfill", {LogicalType::VARCHAR, LogicalType::INTEGER}, LogicalType::LIST(LogicalType::VARCHAR), Geosquarepolyfill);
    ExtensionUtil::RegisterFunction(instance, geosquare_polyfill);

    auto geosquare_polyfill_full = ScalarFunction("geosquare_polyfill_full", {LogicalType::VARCHAR, LogicalType::INTEGER}, LogicalType::LIST(LogicalType::VARCHAR), GeosquarepolyfillFull);
    ExtensionUtil::RegisterFunction(instance, geosquare_polyfill_full);

    auto geosquare_count_gids = ScalarFunction("geosquare_count_gids", {LogicalType::VARCHAR, LogicalType::INTEGER}, LogicalType::INTEGER, Geosquarecount_gids);
    ExtensionUtil::RegisterFunction(instance, geosquare_count_gids);
}

void GeosquareExtension::Load(DuckDB &db) {
	LoadInternal(*db.instance);
}
std::string GeosquareExtension::Name() {
	return "geosquare";
}

std::string GeosquareExtension::Version() const {
#ifdef EXT_VERSION_geosquare
	return EXT_VERSION_geosquare;
#else
	return "";
#endif
}

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void geosquare_init(duckdb::DatabaseInstance &db) {
    duckdb::DuckDB db_wrapper(db);
    db_wrapper.LoadExtension<duckdb::GeosquareExtension>();
}

DUCKDB_EXTENSION_API const char *geosquare_version() {
	return duckdb::DuckDB::LibraryVersion();
}
}

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif