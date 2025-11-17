#include "rclcpp/rclcpp.hpp"
#include "rules_cc/cc/runfiles/runfiles.h"
#include "operational-space-control/walter_sr_air/osc_node.h"
#include "rclcpp/executors/multi_threaded_executor.hpp"

int main(int argc, char** argv) {
    rclcpp::init(argc, argv);

    std::string error;
    std::unique_ptr<rules_cc::cc::runfiles::Runfiles> runfiles(
        rules_cc::cc::runfiles::Runfiles::Create(argv[0], "osc-floating-base", &error));

    if (!error.empty()) {
        std::cerr << "Failed to create runfiles: " << error << std::endl;
        return 1;
    }

    std::filesystem::path model_path = 
        runfiles->Rlocation("mujoco-models+/models/walter_sr/scene_walter_sr_updated_air.xml");

    // --- CONCURRENCY CHANGE ---
    // 1. Instantiate the MultiThreaded Executor
    rclcpp::executors::MultiThreadedExecutor executor;        

    try {
        auto node = std::make_shared<OSCNode>(model_path.string());
        // 2. Add the node to the MultiThreaded Executor
        executor.add_node(node);
        
        // 3. Spin the MultiThreaded Executor
        RCLCPP_INFO(rclcpp::get_logger("main"), "Starting MultiThreaded Executor to fix thread starvation...");
        executor.spin(); // 👈 Use executor.spin() instead of rclcpp::spin(node)
    } catch (const std::exception& e) {
        RCLCPP_FATAL(rclcpp::get_logger("main"), "Exception caught: %s", e.what());
        // Ensure the executor is stopped if an exception occurs
        executor.cancel(); 
        rclcpp::shutdown();        
        return 1;
    }

    rclcpp::shutdown();
    return 0;
}