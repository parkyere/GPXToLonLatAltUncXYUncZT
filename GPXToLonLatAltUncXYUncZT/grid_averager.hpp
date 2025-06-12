#pragma once
#include "Common.hpp"
#include "MapConstruction.hpp"
#include <barrier>
/*------------------------------------------------------------
 *  GridAverager – in‑place 10 m×10 m grid averaging of ECEF
 *-----------------------------------------------------------*/
class GridAverager
{
public:
    explicit GridAverager(BoundingBoxWithPoseFileInfo bBox, double cellSizeMetres = 15.0) noexcept
        : cellSize_(cellSizeMetres), bBox{bBox} {
    }
    void PerformAveragingOperation()
    {
        for (const auto& pF : bBox.pFiles) {
            if ((!(pF.curve)) || pF.curve->empty()) {
                continue;
            }
            if (pF.curve->size() < 2) {
                continue;
            }
            buildCells(*(pF.curve));
        }
        // 3.  Parallel pass: average each cell & write back
        parallelAverage();
        // 4.  Done – memory in map will be released on destruction.
    }

private:
    BoundingBoxWithPoseFileInfo bBox;
    std::map<std::pair<size_t, size_t>, std::vector<VertexPtr>> grid;
private:
    void buildCells(const std::vector<VertexPtr>& pts)
    {
        for (std::size_t i = 0; i < pts.size(); ++i)
        {
            const auto& p = pts[i];
            const size_t gx = static_cast<size_t>((std::floorl((long double)(p->getX()-bBox.minX))/cellSize_));
            const size_t gy = static_cast<size_t>((std::floorl((long double)(p->getY()-bBox.minY))/cellSize_));
            if (grid.contains({gx,gy})) grid[{gx, gy}].push_back(p);
			else grid[{gx, gy}] = { p };
        }
    }

    void parallelAverage()
    {
        // Collect raw pointers to each cell for cheap indexing
        const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
        const unsigned pool = std::max(1u, hw - 1);
		std::vector<std::vector<std::vector<VertexPtr>*>> work(pool,std::vector<std::vector<VertexPtr>*>()); // works for pools to divide work evenly
        std::size_t next{ 0 };
        for (auto& kv : grid) work[next++ % pool].push_back(&kv.second);
        std::vector<std::jthread> threads;
        threads.reserve(pool);
        std::barrier sync_point(pool);
        for (unsigned t = 0; t < pool; ++t)
        {
            threads.emplace_back([&, t]
                {
                    // Required: Each thread processes its own set of cells
                    for (auto& my : work[t])
                    {
						int n = static_cast<int>(my->size());
						if (n < 2) continue; // Skip cells with less than 2 points
						double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
                        for (auto& v : *my)
                        {
                            // Required: Accumulate values for averaging
                            sumX += v->getX();
                            sumY += v->getY();
                            sumZ += v->getZ();
                        }
                        for (auto& v : *my) 
                        {
                            // Required: Average each cell and write back to the original vector
                            v->setX(sumX / n);
                            v->setY(sumY / n);
							v->setZ(sumZ / n);
                        }
                    }
					sync_point.arrive_and_wait(); // Synchronize threads
                }); // jthread auto‑joins on destruction
        }
    }
    double cellSize_;      // metres (default 10)
};