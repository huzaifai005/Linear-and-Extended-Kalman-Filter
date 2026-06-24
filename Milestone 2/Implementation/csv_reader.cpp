// csv_reader.cpp
// Team Thunder — Kalman Filter Milestone 2

#include "csv_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <cmath>

// Joint name table this is the header for all the joint names in order as were given in the true/noisy datasets
const std::string JOINT_NAMES[NUM_JOINTS] = {
    
    "pelvis",         //  0
    "L5",             //  1
    "L3",             //  2
    "T12",            //  3
    "T8",             //  4
    "neck",           //  5
    "head",           //  6
    "shoulderRight",  //  7
    "upperArmRight",  //  8
    "forearmRight",   //  9
    "handRight",      // 10
    "shoulderLeft",   // 11
    "upperArmLeft",   // 12
    "forearmLeft",    // 13
    "handLeft",       // 14
    "upperLegRight",  // 15
    "lowerLegRight",  // 16
    "footRight",      // 17
    "toeRight",       // 18
    "upperLegLeft",   // 19
    "lowerLegLeft",   // 20
    "footLeft",       // 21
    "toeLeft"         // 22
};

// GaitFrame helpers
std::array<double, 3> GaitFrame::get_xyz(int joint_idx) const {
    
    if (joint_idx < 0 || joint_idx >= NUM_JOINTS)
        
        throw std::out_of_range("get_xyz: joint index out of range");
    
    return { joints[joint_idx].x, joints[joint_idx].y, joints[joint_idx].z };
}

std::vector<double> GaitFrame::to_measurement_vector() const {

    std::vector<double> z;
    z.reserve(MEASUREMENT_SIZE);
    
    for (int j = 0; j < NUM_JOINTS; ++j) {
        
        z.push_back(joints[j].x);
        z.push_back(joints[j].y);
        z.push_back(joints[j].z);
    }
    return z;
}

// GaitDataset helpers
std::vector<JointPosition> GaitDataset::get_joint_trajectory(int joint_idx) const {
    
    if (joint_idx < 0 || joint_idx >= NUM_JOINTS)
        
        throw std::out_of_range("get_joint_trajectory: joint index out of range");
    
    std::vector<JointPosition> traj;
    traj.reserve(frames.size());
    
    for (const auto& f : frames)
        
        traj.push_back(f.joints[joint_idx]);
    
    return traj;
}

std::vector<double> GaitDataset::get_measurement(int frame_idx) const {
    
    if (frame_idx < 0 || frame_idx >= num_frames())
        
        throw std::out_of_range("get_measurement: frame index out of range");
    
    return frames[frame_idx].to_measurement_vector();
}

void GaitDataset::print_summary() const {
    
    std::cout << "================================================\n";
    std::cout << "  GaitDataset Summary\n";
    std::cout << "  File    : " << source_file << "\n";
    std::cout << "  Frames  : " << num_frames() << "\n";
    std::cout << "  dt      : " << dt << " s  ("
              << static_cast<int>(std::round(1.0 / dt)) << " Hz)\n";
    std::cout << "  Duration: " << num_frames() * dt << " s\n";
    std::cout << "  Joints  : " << NUM_JOINTS << "\n";
    std::cout << "  Meas/frame: " << MEASUREMENT_SIZE << "\n";

    if (!frames.empty()) {
        
        // Print position range for pelvis (joint 0) to confirm correction
        
        double min_x = 1e18, max_x = -1e18;
        double min_y = 1e18, max_y = -1e18;
        double min_z = 1e18, max_z = -1e18;
        
        for (const auto& f : frames) {
            
            min_x = std::min(min_x, f.joints[0].x);
            max_x = std::max(max_x, f.joints[0].x);
            min_y = std::min(min_y, f.joints[0].y);
            max_y = std::max(max_y, f.joints[0].y);
            min_z = std::min(min_z, f.joints[0].z);
            max_z = std::max(max_z, f.joints[0].z);
        }
        std::cout << "  Pelvis x: [" << std::setprecision(4) << min_x
                  << ", " << max_x << "] m\n";
        std::cout << "  Pelvis y: [" << min_y << ", " << max_y << "] m\n";
        std::cout << "  Pelvis z: [" << min_z << ", " << max_z << "] m\n";
    }
    std::cout << "================================================\n\n";
}

// CSV line parser
std::vector<double> parse_csv_line(const std::string& line) {
    
    std::vector<double> values;
    
    if (line.empty()) 
        
        return values;

    std::stringstream ss(line);
    std::string token;

    while (std::getline(ss, token, ',')) {
       
        // remove whitespace and carriage returns (Windows line endings)
        token.erase(0, token.find_first_not_of(" \t\r\n"));
        token.erase(token.find_last_not_of(" \t\r\n") + 1);

        if (token.empty()) {
            
            values.push_back(0.0);   // empty cell is 0
            continue;
        }

        try {
            
            values.push_back(std::stod(token));
        } catch (const std::invalid_argument&) {
            // Non-numeric token could be a header cell
            // Return empty vector so the caller can detect a header row
            return {};

        } catch (const std::out_of_range&) {
            
            values.push_back(0.0);   // overflow also treat as 0
        }
    }
    return values;
}

// get_joint_index
int get_joint_index(const std::string& joint_name) {
    
    for (int i = 0; i < NUM_JOINTS; ++i)
        
        if (JOINT_NAMES[i] == joint_name)
            
            return i;
    
    return -1;
}

// load_gait_csv
GaitDataset load_gait_csv(const std::string& filename,double dt, bool has_header) {
    
    std::ifstream file(filename);
    
    if (!file.is_open())
        
        throw std::runtime_error("load_gait_csv: cannot open file '" + filename + "'");

    GaitDataset dataset;
    dataset.source_file = filename;
    dataset.dt          = dt;
    dataset.has_header  = has_header;

    std::string line;
    int line_number = 0;
    int data_rows   = 0;

    while (std::getline(file, line)) {
       
        ++line_number;

        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
            
            continue;
        
        // for header
        if (line_number == 1) {
            
            std::vector<double> trial = parse_csv_line(line);
            
            if (trial.empty()) {

                std::cout << "[csv_reader] Detected header row in '"
                          << filename << "', skipping.\n";
                continue;
            }
        }

        // parse data row
        std::vector<double> values = parse_csv_line(line);

        if (values.empty()) {
            
            continue;
        }

        if (static_cast<int>(values.size()) != MEASUREMENT_SIZE) {
            
            std::cerr << "[csv_reader] WARNING: line " << line_number
                      << " has " << values.size() << " values (expected "
                      << MEASUREMENT_SIZE << "). Skipping.\n";
            continue;
        }

        // pack into GaitFrame
        
        GaitFrame frame;
        
        for (int j = 0; j < NUM_JOINTS; ++j) {
            
            int base        = j * DIMS_PER_JOINT;
            frame.joints[j].x = values[base + 0];
            frame.joints[j].y = values[base + 1];
            frame.joints[j].z = values[base + 2];
        }

        dataset.frames.push_back(frame);
        ++data_rows;
    }

    file.close();

    if (data_rows == 0)
        
        throw std::runtime_error("load_gait_csv: no valid data rows found in '" + filename + "'");

    std::cout << "[csv_reader] Loaded " << data_rows
              << " frames from '" << filename << "'.\n";
    return dataset;
}

// load_both_datasets
std::pair<GaitDataset, GaitDataset>

load_both_datasets(const std::string& true_file, const std::string& noisy_file, double dt) {
    
    std::cout << "\n[csv_reader] Loading true data...\n";
    GaitDataset true_data  = load_gait_csv(true_file,  dt);

    std::cout << "[csv_reader] Loading noisy data...\n";
    GaitDataset noisy_data = load_gait_csv(noisy_file, dt);

    // correction check: both files should have the same number of frames
    
    if (true_data.num_frames() != noisy_data.num_frames()) {
        
        std::cerr << "[csv_reader] WARNING: true data has " << true_data.num_frames() << " frames but noisy data has "
                  << noisy_data.num_frames() << " frames.\n" << "             Using the smaller count.\n";

        int n = std::min(true_data.num_frames(), noisy_data.num_frames());
        true_data.frames.resize(n);
        noisy_data.frames.resize(n);
    }

    return { true_data, noisy_data };
}
