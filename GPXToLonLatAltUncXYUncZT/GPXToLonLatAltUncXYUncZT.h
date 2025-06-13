#pragma once
#include "Common.hpp"
void IterateDirectoryAndParse(const std::filesystem::directory_entry& entry, std::filesystem::path& target, double& minX, double& minY, double& maxX, double& maxY, std::vector<PoseFile>& poseFiles);

inline bool NotValidArea(double minX, double minY, double maxX, double maxY) noexcept;
