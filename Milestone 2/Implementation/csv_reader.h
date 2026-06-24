#ifndef CSV_READER_H
#define CSV_READER_H

#include <vector>
#include <string>
#include <array>

// 23 joints × 3 axes = 69 columns

static constexpr int NUM_JOINTS = 23;
static constexpr int DIMS_PER_JOINT = 3;
static constexpr int MEASUREMENT_SIZE = NUM_JOINTS * DIMS_PER_JOINT;

// Joint names in CSV column order: pelvis, L5, L3, T12, T8, neck, head, and so on
extern const std::string JOINT_NAMES[NUM_JOINTS];

struct JointPosition {
    
    double x, y, z;
};

struct GaitFrame {
    
    std::array<JointPosition, NUM_JOINTS> joints;

    // Get position as [x, y, z] array
    std::array<double, 3> get_xyz(int joint_idx) const;
    
    // Flatten to 69-element measurement vector
    std::vector<double> to_measurement_vector() const;
};

struct GaitDataset {
    
    std::vector<GaitFrame> frames;
    double dt;              // sampling interval in seconds
    std::string source_file;
    bool has_header;

    int num_frames() const { 
        
        return static_cast<int>(frames.size()); 
    
    }
    std::vector<JointPosition> get_joint_trajectory(int joint_idx) const;
    std::vector<double> get_measurement(int frame_idx) const;
    
    void print_summary() const;
};

GaitDataset load_gait_csv(const std::string& filename,
    
    double dt = 0.01,
    bool has_header = true);

std::pair<GaitDataset, GaitDataset>

load_both_datasets(const std::string& true_file, const std::string& noisy_file, double dt = 0.01);

std::vector<double> parse_csv_line(const std::string& line);
int get_joint_index(const std::string& joint_name);

#endif // CSV_READER_H