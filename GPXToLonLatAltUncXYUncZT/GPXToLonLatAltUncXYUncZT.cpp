#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <cctype>
#include <algorithm>
#include "tinyxml2.h"
#include "rts_smoother.h"
#include "MapConstruction.hpp"

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
        return 0;

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
    Record r;
    trkpt->QueryDoubleAttribute("lat", &r.lat);
    trkpt->QueryDoubleAttribute("lon", &r.lon);

    if (auto* ele = trkpt->FirstChildElement("ele"))
        ele->QueryDoubleText(&r.alt);

    if (auto* timeEl = trkpt->FirstChildElement("time")) {
        if (const char* ts = timeEl->GetText())
            r.timestamp = parseTimestamp(ts);
    }

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

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <SourceFolder> <TargetFolder>\n";
        return 1;
    }
    fs::path source = argv[1];
    fs::path target = argv[2];
    if (!fs::exists(target)) {
        fs::create_directories(target);
    }
	std::vector<PoseFile> poseFiles;
    for (const auto& entry : fs::directory_iterator(source)) {
        if (!entry.is_regular_file()) continue;
        std::string ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
        if (ext != ".gpx") continue;

        std::vector<Record> records;
        if (!parseGPX(entry.path(), records)) {
            std::cerr << "Failed to parse " << entry.path() << "\n";
            continue;
        }
        if (records.empty()) {
            std::cerr << "No trackpoints in " << entry.path() << "\n";
            continue;
        }
        fs::path outPath = target / entry.path().filename().replace_extension(".csv");
        std::vector<GPSMeasurement> measurements;
		for (const auto& r : records) {
			GPSMeasurement m;
			m.lat = r.lat;
			m.lon = r.lon;
			m.alt = r.alt;
			m.horiz_std = r.hacc;
			m.vert_std = r.vacc;
			m.timestamp = r.timestamp;
			m.has_std = (r.hacc > 0 && r.vacc > 0);
			measurements.push_back(m);
		}
        estimateUncertainties(measurements);

        std::vector<StepData> steps(measurements.size());
        State state{};
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
        }
		poseFiles.push_back(PoseFile::readFile(entry.path().filename().string(), smoothedSteps));
        if (!writeSmoothedTrajectory(outPath, smoothedRecords)) {
            std::cerr << "Failed to write " << outPath << "\n";
        }
    }

}

