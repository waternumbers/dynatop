#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Function to split a CSV line into fields, handling quotes
std::vector<std::string> parseCSVLine(const std::string &line) {
    std::vector<std::string> result;
    std::string field;
    bool inQuotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];

        if (c == '"') {
            // Toggle inQuotes state or handle escaped quotes
            if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
                field += '"';
                ++i; // Skip escaped quote
            } else {
                inQuotes = !inQuotes;
            }
        } else if (c == ',' && !inQuotes) {
            result.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    result.push_back(field); // Add last field
    return result;
}

int main() {
    std::ifstream file("data.csv");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file.\n";
        return 1;
    }

    std::string line;
    bool isHeader = true;
    int wktColumnIndex = -1;

    while (std::getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines

        std::vector<std::string> fields = parseCSVLine(line);

        if (isHeader) {
            // Find WKT column index
            for (size_t i = 0; i < fields.size(); ++i) {
                if (fields[i] == "wkt") {
                    wktColumnIndex = static_cast<int>(i);
                    break;
                }
            }
            if (wktColumnIndex == -1) {
                std::cerr << "Error: No 'wkt' column found in header.\n";
                return 1;
            }
            isHeader = false;
        } else {
            if (wktColumnIndex >= 0 && wktColumnIndex < (int)fields.size()) {
                std::string wkt = fields[wktColumnIndex];
                std::cout << "WKT Geometry: " << wkt << "\n";
            }
        }
    }

    file.close();
    return 0;
}
