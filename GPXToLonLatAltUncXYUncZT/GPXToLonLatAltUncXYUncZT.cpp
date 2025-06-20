#include "Common.hpp"
#include "tinyxml2.h"
#include "rts_smoother.h"
#include "MapConstruction.hpp"
#include "grid_averager.hpp"
#include "GPXToLonLatAltUncXYUncZT.h"

namespace fs = std::filesystem;
using namespace tinyxml2;


struct Record {
    double lon{0};
    double lat{0};
    double alt{0};
    double hacc{0};
    double vacc{0};
    long long timestamp{0};
};

static long long parseTimestamp(const std::string& ts) {
    std::tm tm{};
    std::istringstream ss(ts);
    ss >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%S");
    if (ss.fail())
        return std::numeric_limits<long long>::min();

    // Skip fractional seconds if present
    if (ss.peek() == '.') {
        ss.get();
        while (std::isdigit(ss.peek()))
            ss.get();
    }

    int sign = 0, h = 0, m = 0;
    if (ss.peek() == 'Z') {
        ss.get();
    } else if (ss.peek() == '+' || ss.peek() == '-') {
        char c = ss.get();
        sign = (c == '+') ? 1 : -1;
        ss >> std::setw(2) >> h;
        if (ss.peek() == ':')
            ss.get();
        ss >> std::setw(2) >> m;
    }

    tm.tm_isdst = 0;
#ifdef _WIN32
    time_t t = _mkgmtime(&tm);
#else
    time_t t = timegm(&tm);
#endif
    if (sign)
        t -= sign * (h * 3600 + m * 60);
    return static_cast<long long>(t);
}

static void processTrkpt(XMLElement* trkpt, std::vector<Record>& records) {
    try 
    {
        Record r;
        trkpt->QueryDoubleAttribute("lat", &r.lat);
        trkpt->QueryDoubleAttribute("lon", &r.lon);

        if (auto* ele = trkpt->FirstChildElement("ele"))
            ele->QueryDoubleText(&r.alt);

        if (auto* timeEl = trkpt->FirstChildElement("time")) {
            if (const char* ts = timeEl->GetText())
                r.timestamp = parseTimestamp(ts);
        }
        if (r.timestamp == std::numeric_limits<long long>::min()) throw std::runtime_error("Failed to parse timestamp in GPX file.");

        if (auto* ext = trkpt->FirstChildElement("extensions")) {
            if (auto* tpe = ext->FirstChildElement("rousen:TrackPointExtension")) {
                if (auto* h = tpe->FirstChildElement("rousen:horizontalAccuracy"))
                    h->QueryDoubleText(&r.hacc);
                if (auto* v = tpe->FirstChildElement("rousen:verticalAccuracy"))
                    v->QueryDoubleText(&r.vacc);
            }
        }
        records.push_back(r);
    }
    catch (std::exception& e) {
        std::cerr << "Error processing trackpoint: " << e.what() << "\n";
        return;
    }
    catch (...) {
        std::cerr << "Unknown error processing trackpoint.\n";
        return;
	}
}

static void traverse(XMLElement* node, std::vector<Record>& records) {
    for (auto* e = node; e; e = e->NextSiblingElement()) {
        if (std::string_view{e->Name()} == "trkpt")
            processTrkpt(e, records);
        if (auto* child = e->FirstChildElement())
            traverse(child, records);
    }
}

static bool parseGPX(const fs::path& path, std::vector<Record>& records) {
    XMLDocument doc;
    if (doc.LoadFile(path.string().c_str()) != XML_SUCCESS)
        return false;
    XMLElement* root = doc.RootElement();
    if (!root) return false;
    traverse(root, records);
    return true;
}

static bool writeSmoothedTrajectory(const fs::path& out, const std::vector<Record>& records) {
    std::ofstream ofs(out);
    if (!ofs) return false;
    //ofs << "Longitude,Latitude,Altitude,Timestamp\n";
    for (const auto& r : records) {
        ofs << r.lon << ' ' << r.lat << ' ' << r.alt << ' ' << r.timestamp << "\n";
    }
    return true;
}

void IterateDirectoryAndParse(const std::filesystem::directory_entry& entry, std::filesystem::path& target, double& minX, double& minY, double& maxX, double& maxY, std::vector<PoseFile>& poseFiles)
{
    std::vector<Record> records;
    if (!parseGPX(entry.path(), records)) throw std::runtime_error("Failed to parse " + entry.path().string() + "\n");
    if (records.size() < 3) throw std::runtime_error("Not enough trackpoints in " + entry.path().string() + " (need at least 3)\n");
	
    std::vector<std::vector<GPSMeasurement>> measurementsByTrack = std::vector<std::vector<GPSMeasurement>>(1, std::vector<GPSMeasurement>());
	int insertionIndex = 0;
	double prev_time = records[0].timestamp-1.0;
    for (const auto& r : records) 
    {
        if (r.timestamp <= prev_time) {
            std::cerr << "Skipping invalid timestamp in " << entry.path() << ": " << r.timestamp << "\n";
            continue;
		}
        if(r.timestamp - prev_time >25.0) {
            if (measurementsByTrack[insertionIndex].size() > 2) {
                insertionIndex++;
                measurementsByTrack.push_back(std::vector<GPSMeasurement>());
            }
            else
				measurementsByTrack[insertionIndex].clear();
		}
        prev_time = r.timestamp;
		GPSMeasurement m{ r.lat, r.lon, r.alt, r.hacc, r.vacc, r.timestamp, (r.hacc > 0 && r.vacc > 0) };
        measurementsByTrack[insertionIndex].push_back(m);
    }
    if (measurementsByTrack[insertionIndex].size() < 3) measurementsByTrack.pop_back();

    for (auto ii = 0; ii < measurementsByTrack.size(); ++ii) {
        auto& measurements = measurementsByTrack[ii];
        estimateUncertainties(measurements);
        std::vector<StepData> steps(measurements.size());
        State state{};
        state.pos = geodeticToECEF(measurements[0].lat, measurements[0].lon, measurements[0].alt);

        Matrix6 P = identity6();
        double accel_var = 1.0; // process noise acceleration variance
        double prev_time = measurements[0].timestamp;

        for (size_t i = 0;i < measurements.size();++i) {
            double t = measurements[i].timestamp;
            double dt = (i == 0) ? 0.0 : t - prev_time;
            prev_time = t;
            if (i > 0) predict(state, P, dt, accel_var);
            steps[i].predicted_state = state;
            steps[i].predicted_cov = P;
            steps[i].dt = dt;
            steps[i].timestamp = t;

            std::array<double, 3> z = geodeticToECEF(measurements[i].lat, measurements[i].lon, measurements[i].alt);
            Matrix3 R{};
            double hstd = measurements[i].has_std ? measurements[i].horiz_std : 5.0;
            double vstd = measurements[i].has_std ? measurements[i].vert_std : 5.0;
            R[0][0] = hstd * hstd;
            R[1][1] = hstd * hstd;
            R[2][2] = vstd * vstd;
            update(state, P, z, R);
            steps[i].filtered_state = state;
            steps[i].filtered_cov = P;
        }

        smooth(steps, accel_var);
        std::vector<Record> smoothedRecords;
        std::vector<StepData> smoothedSteps;

        for (const auto& step : steps) {
            auto llh = ecefToGeodetic(step.filtered_state.pos[0],
                step.filtered_state.pos[1],
                step.filtered_state.pos[2]);
            smoothedRecords.push_back({
                llh[1], // lon
                llh[0], // lat
                llh[2], // alt
                std::sqrt(step.filtered_cov[0][0]), // hacc
                std::sqrt(step.filtered_cov[2][2]), // vacc
                (long long)step.timestamp // timestamp
                });
            smoothedSteps.push_back(step);
            minX = std::min(minX, step.filtered_state.pos[0]);
            minY = std::min(minY, step.filtered_state.pos[1]);
            maxX = std::max(maxX, step.filtered_state.pos[0]);
            maxY = std::max(maxY, step.filtered_state.pos[1]);
        }
        poseFiles.push_back(PoseFile::readFile(entry.path().stem().string()+'_'+ std::to_string(ii)+entry.path().extension().string(), smoothedSteps));
        fs::path outPath = target / (entry.path().stem().string()+'_'+std::to_string(ii) + ".csv");
        if (!writeSmoothedTrajectory(outPath, smoothedRecords)) {
            std::cerr << "Failed to write " << outPath << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
	std::cin.tie(0)->sync_with_stdio(0);
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <SourceFolder> <TargetFolder>\n";
        return 1;
    }
    try {
        fs::path source = argv[1];
        if (!fs::exists(source) || !fs::is_directory(source)) {
            throw std::invalid_argument("Source folder does not exist or is not a directory.\n");
		}
        fs::path target = argv[2];
        if (!fs::exists(target)) {
            fs::create_directories(target);
        }
        std::vector<PoseFile> poseFiles;
//        std::vector<PoseFile> rawPoseFiles; //For debug;
        double minX = std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();
        for (const auto& entry : fs::directory_iterator(source)) {
            if (!entry.is_regular_file()) continue;
            std::string ext = entry.path().extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
            if (ext != ".gpx") continue;
            try {
                IterateDirectoryAndParse(entry, target, minX, minY, maxX, maxY, poseFiles);
            }
            catch (std::exception& e) {
                std::cerr << "Error processing " << entry.path() << ": " << e.what() << "\n";
                continue;
            }
            catch (...) {
                std::cerr << "Unknown error processing " << entry.path() << "\n";
                continue;
            }
        }
        if (NotValidArea(minX, minY, maxX, maxY)) {
            throw std::runtime_error("No valid GPX files found in the source directory.\n");
        }
        MapConstruction mapConstruction;
        std::vector<VertexPtr> constructedMap = mapConstruction.constructMapMain(poseFiles, 15, 15);
        MapConstruction::writeToFile(constructedMap, target.generic_string() + "/out");
        //std::vector<VertexPtr> constructedMap = 
    }
    catch (const std::exception& e) {
        std::cerr << "Error running the program: " << e.what() << "\n";
        return 1;
	}
    catch (...) {
        std::cerr << "Fatal error running the program.\n";
		return 1;
    }

}

inline bool NotValidArea(double minX, double minY, double maxX, double maxY) noexcept
{
    return minX == std::numeric_limits<double>::infinity() ||
        minY == std::numeric_limits<double>::infinity() ||
        maxX == -std::numeric_limits<double>::infinity() ||
        maxY == -std::numeric_limits<double>::infinity();
}

