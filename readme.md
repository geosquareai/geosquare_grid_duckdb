# Geosquare DuckDB Extension

## About

`geosquare` is a fast grid-based spatial indexing extension for DuckDB. It allows you to convert spatial coordinates into discrete hierarchical geocodes (GIDs), fetch their corresponding bounding boxes, perform fast point-in-polygon queries, and efficiently intersect geometries across partitioned Parquet datasets directly within your DuckDB SQL queries. 

By avoiding complex intersections and converting spatial coordinates into a string-based grid index, Geosquare can significantly accelerate geospatial lookups and joins, especially on large datasets.

### Grid Size Levels

The `resolution` or `level` parameter in several functions corresponds to the grid's precision. You configure the target cell dimension using the meter-scale equivalent:

| Size (meter) | Level |
|----------|-------|
| 10000000 | 1     | 
| 1000000  | 3     | 
| 100000   | 5     | 
| 10000    | 7     | 
| 5000     | 8     | 
| 1000     | 9     | 
| 500      | 10    | 
| 100      | 11    | 
| 50       | 12    | 
| 10       | 13    | 
| 5        | 14    | 

## Installation

The extension is continuously built and published for versions generated via Git tags. To install and load the pre-built extension for a given published version, execute the following SQL commands in your DuckDB connection. 
Note: Since the extension is community-built, you must first enable `allow_unsigned_extensions`.

```sql
SET allow_unsigned_extensions=true;
INSTALL geosquare FROM 'https://geosquareai.github.io/geosquare_grid_duckdb';
LOAD geosquare;
```

For instance, within a Javascript/TypeScript environment using the DuckDB Node.js API:

```javascript
await db.run("SET allow_unsigned_extensions=true");
await db.run("INSTALL geosquare FROM 'https://geosquareai.github.io/geosquare_grid_duckdb'");
await db.run("LOAD geosquare");
```

### Local Build

If a pre-built binary is not available for your platform or DuckDB version, you must build the extension yourself from the source tree:

```bash
# Clone the repository
git clone https://github.com/geosquareai/geosquare_grid.duckdb_extension.git
cd geosquare_grid.duckdb_extension

# Build the extension using Makefile (Requires CMake)
make release
```

Once built, you can load your locally compiled binary into DuckDB directly by dropping the `.duckdb_extension` file into your DuckDB application, or loading it dynamically:

```sql
LOAD 'build/release/extension/geosquare/geosquare.duckdb_extension';
```

## How to Use & Examples

Below is a comprehensive list of the core functions provided by the extension.

### Core Functions

#### Geosquare Polyfill (`geosquare_polyfill`)
Retrieves a list of GIDs that cover a given WKT geometry at a specified resolution size (in meters). Set the boolean flag `fullcover` to `true` to include partially intersecting grid cells, or `false` to include only cells whose majority are covered.
```sql
SELECT geosquare_polyfill('POLYGON ((106.851 -6.288, 106.860 -6.288, 106.860 -6.279, 106.851 -6.279, 106.851 -6.288))', 100, true) AS gids;
```

#### Coordinate to Grid (`geosquare_lonlat_to_gid`)
Convert a longitude and latitude coordinate to a Geosquare GID at a specific level (1-15).
```sql
SELECT geosquare_lonlat_to_gid(106.8516, -6.2882, 11) AS gid;
```

#### Grid to Coordinate (`geosquare_gid_to_lonlat`)
Returns the centroid `[longitude, latitude]` for a given GID.
```sql
SELECT geosquare_gid_to_lonlat('2M74') AS lonlat;
```

#### Grid to Point WKT (`geosquare_gid_to_point_wkt`)
Returns the centroid as a WKT Point geometry string.
```sql
SELECT geosquare_gid_to_point_wkt('2M74') AS point_geom;
```

#### Grid to Bounding Box (`geosquare_gid_to_bound`)
Returns the `[west, east, south, north]` bounding box of a GID.
```sql
SELECT geosquare_gid_to_bound('2M74') AS bbox;
```

#### Grid to Bounding Box WKT (`geosquare_gid_to_bound_wkt`)
Returns the bounding box of the grid cell as a WKT Polygon string.
```sql
SELECT geosquare_gid_to_bound_wkt('2M74') AS polygon_geom;
```

#### Grid Hierarchy (`geosquare_gid_parent`, `geosquare_gid_children`)
Navigate the hierarchical grid by getting the parent or children of a given GID.
```sql
-- Get the parent cell (one level up)
SELECT geosquare_gid_parent('2M74') AS parent_gid;

-- Get all children cells (one level down)
SELECT geosquare_gid_children('2M74') AS child_gids;
```

#### Count GIDs (`geosquare_count_gids`)
Returns the total number of grid cells required to cover a given geometry at the specified resolution size.
```sql
SELECT geosquare_count_gids('POLYGON ((...))', 100) AS count;
```

### Table Macros

The extension also offers dynamic macros to simplify partitioned query operations.

#### Intersect Parquet (`geosquare_intersect_parquet`)
Dynamically reads from partitioned Parquet files (partitioned by `gid` folder structure e.g. `gid=2MH4/...`) and returns rows that intersect a given target geometry.

```sql
-- Using the default geometry column (`geometry`)
SELECT * FROM geosquare_intersect_parquet(
    ST_GeomFromText('POLYGON(...)'), 
    's3://my-bucket/dataset', 
    100
);

-- Specifying a custom geometry column (`geom_col`)
SELECT * FROM geosquare_intersect_parquet(
    ST_GeomFromText('POLYGON(...)'), 
    's3://my-bucket/dataset', 
    100,
    my_custom_geom_col
);
```