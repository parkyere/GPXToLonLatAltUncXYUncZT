#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include "tinyxml2.h"

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
    if (strptime(ts.c_str(), "%Y-%m-%dT%H:%M:%SZ", &tm)) {
        time_t t = timegm(&tm);
        return static_cast<long long>(t);
    }
    return 0;
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

static bool writeCSV(const fs::path& out, const std::vector<Record>& records) {
    std::ofstream ofs(out);
    if (!ofs) return false;
    ofs << "Longitude,Latitude,Altitude,XYUncertainty,ZUncertainty,Timestamp\n";
    for (const auto& r : records) {
        ofs << r.lon << ',' << r.lat << ',' << r.alt << ',' << r.hacc << ','
            << r.vacc << ',' << r.timestamp << "\n";
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

    for (const auto& entry : fs::directory_iterator(source)) {
        if (!entry.is_regular_file()) continue;
        auto ext = entry.path().extension();
        if (ext != ".gpx" && ext != ".GPX") continue;

        std::vector<Record> records;
        if (!parseGPX(entry.path(), records)) {
            std::cerr << "Failed to parse " << entry.path() << "\n";
            continue;
        }
        fs::path outPath = target / entry.path().filename().replace_extension(".csv");
        if (!writeCSV(outPath, records)) {
            std::cerr << "Failed to write " << outPath << "\n";
        }
    }
    return 0;
}

