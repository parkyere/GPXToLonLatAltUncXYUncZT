#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <libxml/parser.h>
#include <libxml/tree.h>

namespace fs = std::filesystem;

struct Record {
    double lon{0};
    double lat{0};
    double alt{0};
    double hacc{0};
    double vacc{0};
    long long timestamp{0};
};

static double getDoubleProp(xmlNodePtr node, const char* name) {
    xmlChar* prop = xmlGetProp(node, BAD_CAST name);
    if (!prop) return 0.0;
    double val = atof(reinterpret_cast<const char*>(prop));
    xmlFree(prop);
    return val;
}

static std::string getChildContent(xmlNodePtr parent, const char* name) {
    for (xmlNodePtr n = parent->children; n; n = n->next) {
        if (n->type == XML_ELEMENT_NODE && xmlStrEqual(n->name, BAD_CAST name)) {
            xmlChar* content = xmlNodeGetContent(n);
            if (!content) return "";
            std::string s(reinterpret_cast<const char*>(content));
            xmlFree(content);
            return s;
        }
    }
    return "";
}

static xmlNodePtr findChild(xmlNodePtr parent, const char* name) {
    for (xmlNodePtr n = parent->children; n; n = n->next) {
        if (n->type == XML_ELEMENT_NODE && xmlStrEqual(n->name, BAD_CAST name)) {
            return n;
        }
    }
    return nullptr;
}

static long long parseTimestamp(const std::string& ts) {
    struct tm tm{};
    if (strptime(ts.c_str(), "%Y-%m-%dT%H:%M:%SZ", &tm)) {
        time_t t = timegm(&tm);
        return static_cast<long long>(t);
    }
    return 0;
}

static void processTrkpt(xmlNodePtr trkpt, std::vector<Record>& records) {
    Record r;
    r.lat = getDoubleProp(trkpt, "lat");
    r.lon = getDoubleProp(trkpt, "lon");
    std::string ele = getChildContent(trkpt, "ele");
    if (!ele.empty()) r.alt = atof(ele.c_str());
    std::string timestr = getChildContent(trkpt, "time");
    if (!timestr.empty()) r.timestamp = parseTimestamp(timestr);

    xmlNodePtr extNode = findChild(trkpt, "extensions");
    if (extNode) {
        for (xmlNodePtr e = extNode->children; e; e = e->next) {
            if (e->type != XML_ELEMENT_NODE) continue;
            if (xmlStrEqual(e->name, BAD_CAST "TrackPointExtension")) {
                std::string hacc = getChildContent(e, "horizontalAccuracy");
                if (!hacc.empty()) r.hacc = atof(hacc.c_str());
                std::string vacc = getChildContent(e, "verticalAccuracy");
                if (!vacc.empty()) r.vacc = atof(vacc.c_str());
            }
        }
    }
    records.push_back(r);
}

static bool parseGPX(const fs::path& path, std::vector<Record>& records) {
    xmlDocPtr doc = xmlReadFile(path.string().c_str(), nullptr, 0);
    if (!doc) return false;
    xmlNodePtr root = xmlDocGetRootElement(doc);
    if (!root) { xmlFreeDoc(doc); return false; }

    // traverse to find trkpt
    std::vector<xmlNodePtr> stack{root};
    while (!stack.empty()) {
        xmlNodePtr node = stack.back();
        stack.pop_back();
        if (node->type == XML_ELEMENT_NODE && xmlStrEqual(node->name, BAD_CAST "trkpt")) {
            processTrkpt(node, records);
        }
        for (xmlNodePtr child = node->children; child; child = child->next) {
            if (child->type == XML_ELEMENT_NODE) stack.push_back(child);
        }
    }
    xmlFreeDoc(doc);
    return true;
}

static bool writeCSV(const fs::path& out, const std::vector<Record>& records) {
    std::ofstream ofs(out);
    if (!ofs) return false;
    ofs << "Longitude,Latitude,Altitude,XYUncertainty,ZUncertainty,Timestamp\n";
    for (const auto& r : records) {
        ofs << r.lon << ',' << r.lat << ',' << r.alt << ',' << r.hacc << ',' << r.vacc << ',' << r.timestamp << "\n";
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
        if (entry.path().extension() != ".gpx" && entry.path().extension() != ".GPX") continue;
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

