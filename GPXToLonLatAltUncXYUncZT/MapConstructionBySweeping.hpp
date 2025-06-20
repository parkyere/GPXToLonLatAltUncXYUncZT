#pragma once
#include "Common.hpp"

namespace traj {

    /* ======== BASIC GEOMETRY PRIMITIVES ======== */

    struct Vertex {
        double x, y, z, t;
    };

    struct Edge {
        std::size_t a;          // index of first  vertex in vertices[]
        std::size_t b;          // index of second vertex
        double minx, miny;      // cached AABB
        double maxx, maxy;
    };

    // squared distance from point P to segment AB (2‑D)
    inline double pointSegmentDist2(double px, double py,
        double ax, double ay,
        double bx, double by) noexcept
    {
        const double vx = bx - ax, vy = by - ay;
        const double wx = px - ax, wy = py - ay;
        const double vv = vx * vx + vy * vy;

        double t = (vv > 0.0) ? (wx * vx + wy * vy) / vv : 0.0;
        t = std::clamp(t, 0.0, 1.0);

        const double dx = wx - t * vx;
        const double dy = wy - t * vy;
        return dx * dx + dy * dy;
    }

    /* ======== HASHED GRID KEY ======== */

    struct CellKey {
        int gx, gy;                         // integer grid coordinates
        bool operator==(const CellKey& k) const noexcept
        {
            return gx == k.gx && gy == k.gy;
        }
    };

    struct CellKeyHash {
        std::size_t operator()(const CellKey& k) const noexcept
        {
            // 64‑bit mix (Thomas Wang)
            uint64_t x = (uint64_t)k.gx;
            uint64_t y = (uint64_t)k.gy * 0x9E3779B97F4A7C15ULL; // φ
            uint64_t z = x ^ y;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            return (std::size_t)(z ^ (z >> 31));
        }
    };

    /* ======== MAIN CLASS ======== */

    class SpatialHashGrid
    {
    public:
        explicit SpatialHashGrid(double cellSize)
            : cell(cellSize), invCell(1.0 / cellSize),
            epsSqCache(-1.0), // invalid
            candVerts(), candEdges()
        {
        }

        /* ---- insertion ---- */

        std::size_t addVertex(double x, double y, double z, double t)
        {
            std::size_t id = vertices.size();
            vertices.push_back({ x,y,z,t });
            insertVertexToGrid(id, x, y);
            return id;
        }

        std::size_t addEdge(std::size_t a, std::size_t b)
        {
            const auto& A = vertices[a];
            const auto& B = vertices[b];

            Edge e{ a,b,
                   std::min(A.x,B.x), std::min(A.y,B.y),
                   std::max(A.x,B.x), std::max(A.y,B.y) };

            std::size_t id = edges.size();
            edges.push_back(e);
            insertEdgeToGrid(id, e);
            return id;
        }

        /* ---- query ---- */
        void query(double qx, double qy, double eps,
            std::vector<std::size_t>& outVerts,
            std::vector<std::size_t>& outEdges)
        {
            outVerts.clear();
            outEdges.clear();
            candVerts.clear();
            candEdges.clear();

            // gather candidates ---------------------------------------------------
            const int gx0 = gridCoord(qx - eps);
            const int gx1 = gridCoord(qx + eps);
            const int gy0 = gridCoord(qy - eps);
            const int gy1 = gridCoord(qy + eps);

            for (int gx = gx0; gx <= gx1; ++gx)
                for (int gy = gy0; gy <= gy1; ++gy)
                {
                    CellKey key{ gx,gy };
                    if (auto it = vIndex.find(key); it != vIndex.end())
                        candVerts.insert(candVerts.end(),
                            it->second.begin(), it->second.end());
                    if (auto it = eIndex.find(key); it != eIndex.end())
                        candEdges.insert(candEdges.end(),
                            it->second.begin(), it->second.end());
                }

            // eliminate duplicates (linear because typical candidate counts are low)
            std::sort(candVerts.begin(), candVerts.end());
            candVerts.erase(std::unique(candVerts.begin(), candVerts.end()),
                candVerts.end());

            std::sort(candEdges.begin(), candEdges.end());
            candEdges.erase(std::unique(candEdges.begin(), candEdges.end()),
                candEdges.end());

            // exact test ----------------------------------------------------------
            const double epsSq = eps * eps;
            for (std::size_t vid : candVerts)
            {
                const auto& V = vertices[vid];
                if ((V.x - qx) * (V.x - qx) + (V.y - qy) * (V.y - qy) <= epsSq)
                    outVerts.push_back(vid);
            }
            for (std::size_t eid : candEdges)
            {
                const auto& E = edges[eid];
                // quick AABB reject
                if (E.minx > qx + eps || E.maxx < qx - eps ||
                    E.miny > qy + eps || E.maxy < qy - eps)
                    continue;

                const auto& A = vertices[E.a];
                const auto& B = vertices[E.b];
                const double d2 = pointSegmentDist2(qx, qy, A.x, A.y, B.x, B.y);
                if (d2 <= epsSq)
                    outEdges.push_back(eid);
            }
        }

        /* ---- accessors (optional) ---- */
        const std::vector<Vertex>& getVertices() const noexcept { return vertices; }
        const std::vector<Edge>& getEdges()   const noexcept { return edges; }

    private:
        /* helpers ----------------------------------------------------------- */
        int gridCoord(double x) const noexcept
        {
            return static_cast<int>(std::floor(x * invCell));
        }

        CellKey makeKey(double x, double y) const noexcept
        {
            return { gridCoord(x), gridCoord(y) };
        }

        void insertVertexToGrid(std::size_t id, double x, double y)
        {
            vIndex[makeKey(x, y)].push_back(id);
        }

        void insertEdgeToGrid(std::size_t id, const Edge& e)
        {
            const int gx0 = gridCoord(e.minx);
            const int gx1 = gridCoord(e.maxx);
            const int gy0 = gridCoord(e.miny);
            const int gy1 = gridCoord(e.maxy);

            for (int gx = gx0; gx <= gx1; ++gx)
                for (int gy = gy0; gy <= gy1; ++gy)
                    eIndex[{gx, gy}].push_back(id);
        }

        /* data ----------------------------------------------------------------- */
        double cell, invCell;
        double epsSqCache;                   // for potential reuse

        std::vector<Vertex> vertices;
        std::vector<Edge>   edges;

        std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> vIndex;
        std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> eIndex;

        // re‑usable scratch buffers to avoid re‑allocating every query
        std::vector<std::size_t> candVerts;
        std::vector<std::size_t> candEdges;
    };

} // namespace traj
