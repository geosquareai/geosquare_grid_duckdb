#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cctype>

namespace duckdb {

struct Point2D {
    double x, y;
};

struct LinearRing {
    std::vector<Point2D> points;
};

struct SimplePolygon {
    LinearRing exterior;
    std::vector<LinearRing> holes;
};

struct LightweightGeometry {
    std::vector<SimplePolygon> polygons;
    bool is_empty() const { return polygons.empty(); }
    
    // Envelope calculation
    void get_envelope(double& min_x, double& min_y, double& max_x, double& max_y) const {
        if (polygons.empty() || polygons[0].exterior.points.empty()) {
            min_x = min_y = max_x = max_y = 0;
            return;
        }
        min_x = max_x = polygons[0].exterior.points[0].x;
        min_y = max_y = polygons[0].exterior.points[0].y;
        
        for (const auto& poly : polygons) {
            for (const auto& pt : poly.exterior.points) {
                if (pt.x < min_x) min_x = pt.x;
                if (pt.x > max_x) max_x = pt.x;
                if (pt.y < min_y) min_y = pt.y;
                if (pt.y > max_y) max_y = pt.y;
            }
        }
    }
};

class LightweightWKTParser {
public:
    static LightweightGeometry parse(const std::string& wkt) {
        LightweightGeometry geom;
        std::string str = wkt;
        // make uppercase
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
        
        size_t pos = 0;
        skip_whitespace(str, pos);
        
        if (pos + 7 <= str.length() && str.substr(pos, 7) == "POLYGON") {
            pos += 7;
            geom.polygons.push_back(parse_polygon(str, pos));
        } else if (pos + 12 <= str.length() && str.substr(pos, 12) == "MULTIPOLYGON") {
            pos += 12;
            geom.polygons = parse_multi_polygon(str, pos);
        } else {
            throw std::runtime_error("Only POLYGON and MULTIPOLYGON are supported. wkt: " + str);
        }
        return geom;
    }

private:
    static void skip_whitespace(const std::string& str, size_t& pos) {
        while (pos < str.length() && std::isspace(static_cast<unsigned char>(str[pos]))) pos++;
    }
    
    static void expect_char(const std::string& str, size_t& pos, char c) {
        skip_whitespace(str, pos);
        if (pos >= str.length() || str[pos] != c) {
            throw std::runtime_error(std::string("Expected '") + c + "'");
        }
        pos++;
    }
    
    static Point2D parse_point(const std::string& str, size_t& pos) {
        skip_whitespace(str, pos);
        size_t next_pos;
        double x = std::stod(str.substr(pos), &next_pos);
        pos += next_pos;
        skip_whitespace(str, pos);
        double y = std::stod(str.substr(pos), &next_pos);
        pos += next_pos;
        return {x, y};
    }
    
    static LinearRing parse_ring(const std::string& str, size_t& pos) {
        expect_char(str, pos, '(');
        LinearRing ring;
        while (pos < str.length()) {
            ring.points.push_back(parse_point(str, pos));
            skip_whitespace(str, pos);
            if (pos < str.length() && str[pos] == ',') {
                pos++;
            } else if (pos < str.length() && str[pos] == ')') {
                pos++;
                break;
            } else {
                throw std::runtime_error("Expected ',' or ')' in ring");
            }
        }
        return ring;
    }
    
    static SimplePolygon parse_polygon(const std::string& str, size_t& pos) {
        expect_char(str, pos, '(');
        SimplePolygon poly;
        poly.exterior = parse_ring(str, pos);
        skip_whitespace(str, pos);
        while (pos < str.length() && str[pos] == ',') {
            pos++;
            poly.holes.push_back(parse_ring(str, pos));
            skip_whitespace(str, pos);
        }
        expect_char(str, pos, ')');
        return poly;
    }
    
    static std::vector<SimplePolygon> parse_multi_polygon(const std::string& str, size_t& pos) {
        expect_char(str, pos, '(');
        std::vector<SimplePolygon> polys;
        while (pos < str.length()) {
            polys.push_back(parse_polygon(str, pos));
            skip_whitespace(str, pos);
            if (pos < str.length() && str[pos] == ',') {
                pos++;
            } else if (pos < str.length() && str[pos] == ')') {
                pos++;
                break;
            } else {
                throw std::runtime_error("Expected ',' or ')' in multipolygon");
            }
        }
        return polys;
    }
};

class LightweightGeometryOperations {
public:
    static double get_area(const LinearRing& ring) {
        if (ring.points.size() < 3) return 0.0;
        double area = 0.0;
        for (size_t i = 0; i < ring.points.size(); i++) {
            size_t j = (i + 1) % ring.points.size();
            area += (ring.points[i].x * ring.points[j].y) - (ring.points[j].x * ring.points[i].y);
        }
        return std::abs(area) / 2.0;
    }

    static double get_area(const SimplePolygon& poly) {
        double area = get_area(poly.exterior);
        for (const auto& hole : poly.holes) {
            area -= get_area(hole);
        }
        return area;
    }

    static double get_area(const LightweightGeometry& geom) {
        double area = 0.0;
        for (const auto& poly : geom.polygons) {
            area += get_area(poly);
        }
        return area;
    }

    // Sutherland-Hodgman clipping against AABB
    static LinearRing clip_ring_to_aabb(const LinearRing& subject, double min_x, double min_y, double max_x, double max_y) {
        if (subject.points.empty()) return subject;
        
        auto clip_edge = [](const std::vector<Point2D>& input, int edge, double edge_val) {
            std::vector<Point2D> output;
            if (input.empty()) return output;
            
            auto inside = [edge, edge_val](const Point2D& p) {
                switch(edge) {
                    case 0: return p.x >= edge_val; // Left
                    case 1: return p.x <= edge_val; // Right
                    case 2: return p.y >= edge_val; // Bottom
                    case 3: return p.y <= edge_val; // Top
                }
                return false;
            };

            auto intersect = [edge, edge_val](const Point2D& p1, const Point2D& p2) {
                Point2D p;
                switch(edge) {
                    case 0: // Left
                    case 1: // Right
                        p.x = edge_val;
                        p.y = p1.y + (p2.y - p1.y) * (edge_val - p1.x) / (p2.x - p1.x);
                        break;
                    case 2: // Bottom
                    case 3: // Top
                        p.y = edge_val;
                        p.x = p1.x + (p2.x - p1.x) * (edge_val - p1.y) / (p2.y - p1.y);
                        break;
                }
                return p;
            };
            
            Point2D s = input.back();
            for (const auto& e : input) {
                if (inside(e)) {
                    if (!inside(s)) {
                        output.push_back(intersect(s, e));
                    }
                    output.push_back(e);
                } else if (inside(s)) {
                    output.push_back(intersect(s, e));
                }
                s = e;
            }
            return output;
        };

        std::vector<Point2D> pts = subject.points;
        pts = clip_edge(pts, 0, min_x); // Left
        pts = clip_edge(pts, 1, max_x); // Right
        pts = clip_edge(pts, 2, min_y); // Bottom
        pts = clip_edge(pts, 3, max_y); // Top
        
        LinearRing result;
        result.points = pts;
        return result;
    }

    static double get_intersection_area(const LightweightGeometry& geom, double min_x, double min_y, double max_x, double max_y) {
        double area = 0.0;
        for (const auto& poly : geom.polygons) {
            LinearRing clipped_ext = clip_ring_to_aabb(poly.exterior, min_x, min_y, max_x, max_y);
            double poly_area = get_area(clipped_ext);
            for (const auto& hole : poly.holes) {
                LinearRing clipped_hole = clip_ring_to_aabb(hole, min_x, min_y, max_x, max_y);
                poly_area -= get_area(clipped_hole);
            }
            area += poly_area;
        }
        return area;
    }
};

} // namespace duckdb
