//Bradley Till, May 15, 2025

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>

struct SampleInfo {
    std::string name;
    char status; // '+', '-', or '0'
};

struct GenotypeCount {
    int homozygousRef = 0;    // 0/0
    int heterozygous = 0;     // 0/1
    int homozygousAlt = 0;    // 1/1
    int total = 0;            // Total valid genotypes (not ./.)
};

// Function to parse the samples file
std::vector<SampleInfo> parseSamplesFile(const std::string& filename) {
    std::vector<SampleInfo> samples;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error opening samples file: " << filename << std::endl;
        return samples;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        SampleInfo sample;
        
        if (iss >> sample.name >> sample.status) {
            samples.push_back(sample);
        }
    }
    
    return samples;
}

// Function to process the data file and calculate frequencies
void processDataFile(const std::string& dataFilename, const std::vector<SampleInfo>& samples) {
    std::ifstream dataFile(dataFilename);
    
    if (!dataFile.is_open()) {
        std::cerr << "Error opening data file: " << dataFilename << std::endl;
        return;
    }
    
    // Read the header line to get sample positions
    std::string headerLine;
    std::getline(dataFile, headerLine);
    
    std::istringstream headerStream(headerLine);
    std::vector<std::string> headers;
    std::string header;
    
    // Read all header columns
    while (headerStream >> header) {
        headers.push_back(header);
    }
    
    // Map sample names to their column index
    std::map<std::string, int> sampleColumns;
    for (size_t i = 2; i < headers.size(); ++i) {  // Start from column 3 (index 2)
        sampleColumns[headers[i]] = i;
    }
    
    // Print output header
    std::cout << "Chrom Pos 0/0_Percent_Case 0/0_Percent_Control 0/1_Percent_Case 0/1_Percent_Control "
              << "1/1_Percent_Case 1/1_Percent_Control" << std::endl;
    
    // Process each data line
    std::string line;
    while (std::getline(dataFile, line)) {
        std::istringstream lineStream(line);
        std::vector<std::string> fields;
        std::string field;
        
        while (lineStream >> field) {
            fields.push_back(field);
        }
        
        if (fields.size() < 3) continue;  // Skip lines with insufficient data
        
        // Initialize counters for cases and controls
        GenotypeCount caseCounts, controlCounts;
        
        // Process each sample's genotype
        for (const auto& sample : samples) {
            // Skip samples with status '0'
            if (sample.status == '0') continue;
            
            // Find the column for this sample
            auto it = sampleColumns.find(sample.name);
            if (it == sampleColumns.end()) continue;  // Sample not found in header
            
            int colIndex = it->second;
            if (colIndex >= fields.size()) continue;  // Column index out of bounds
            
            std::string genotype = fields[colIndex];
            
            // Skip missing genotypes
            if (genotype == "./.") continue;
            
            // Determine which counter to update
            GenotypeCount& counts = (sample.status == '+') ? caseCounts : controlCounts;
            
            // Count the genotype
            if (genotype == "0/0") {
                counts.homozygousRef++;
            } else if (genotype == "0/1" || genotype == "1/0") {
                counts.heterozygous++;
            } else if (genotype == "1/1") {
                counts.homozygousAlt++;
            }
            
            // Update total valid genotypes
            counts.total++;
        }
        
        // Calculate percentages
        auto calculatePercentage = [](int count, int total) -> int {
            if (total == 0) return 0;
            return (count * 100) / total;
        };
        
        int case_hom_ref_percent = calculatePercentage(caseCounts.homozygousRef, caseCounts.total);
        int control_hom_ref_percent = calculatePercentage(controlCounts.homozygousRef, controlCounts.total);
        int case_het_percent = calculatePercentage(caseCounts.heterozygous, caseCounts.total);
        int control_het_percent = calculatePercentage(controlCounts.heterozygous, controlCounts.total);
        int case_hom_alt_percent = calculatePercentage(caseCounts.homozygousAlt, caseCounts.total);
        int control_hom_alt_percent = calculatePercentage(controlCounts.homozygousAlt, controlCounts.total);
        
        // Print the results
        std::cout << fields[0] << " " << fields[1] << " "
                  << case_hom_ref_percent << " " << control_hom_ref_percent << " "
                  << case_het_percent << " " << control_het_percent << " "
                  << case_hom_alt_percent << " " << control_hom_alt_percent << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <data_file> <samples_file>" << std::endl;
        return 1;
    }
    
    std::string dataFilename = argv[1];
    std::string samplesFilename = argv[2];
    
    // Parse sample information
    std::vector<SampleInfo> samples = parseSamplesFile(samplesFilename);
    
    // Process data file and calculate frequencies
    processDataFile(dataFilename, samples);
    
    return 0;
}
