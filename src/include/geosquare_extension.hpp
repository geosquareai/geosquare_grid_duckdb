#pragma once

#include "duckdb.hpp"

namespace duckdb {

/** @struct LatLng
@brief latitude/longitude in radians
*/
typedef struct {
    double lat;  ///< latitude in radians
    double lng;  ///< longitude in radians
} LatLng;

/** @struct BBox
 *  @brief  Geographic bounding box with coordinates defined in radians
 */
typedef struct {
    double north;  ///< north latitude
    double south;  ///< south latitude
    double east;   ///< east longitude
    double west;   ///< west longitude
} BBox;

class GeosquareExtension : public Extension {
public:
	void Load(DuckDB &db) override;
	std::string Name() override;
        std::string Version() const override;
};

#define LAT_RANGED_MIN_INIT -216.0
#define LAT_RANGED_MAX_INIT 233.157642055036
#define LON_RANGED_MIN_INIT -217.0
#define LON_RANGED_MAX_INIT 232.157642055036

} // namespace duckdb