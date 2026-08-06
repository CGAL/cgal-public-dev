#ifndef EXECUTION_TIMING_VISUALIZER_H
#define EXECUTION_TIMING_VISUALIZER_H

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <cmath>

// Polyscope & ImGui Includes
#include "polyscope/polyscope.h"
#include "imgui.h"

// Execution Stats Definition
#include "../testBVH/ExecutionStats.h"

class ExecutionTimingVisualizer {
public:
    struct TimingSegment {
        std::string name;
        double ms = 0.0;
        ImVec4 color;
        std::string description;
    };

private:
    bool m_isVisible = false;
    double m_totalTimeMs = 0.0;
    std::vector<TimingSegment> m_segments;

public:
    ExecutionTimingVisualizer() = default;

    // Registers user callback hook into Polyscope
    void init() {
        // Appends or sets polyscope callback loop
        auto existingCallback = polyscope::state::userCallback;
        polyscope::state::userCallback = [this, existingCallback]() {
            if (existingCallback) existingCallback();
            this->renderUI();
        };
    }

    // Call this inside runComputeLogic or cmdTestBVH to update data and display window
    void updateAndShow(const ExecutionStats& stats) {
        m_segments.clear();
        m_totalTimeMs = stats.GPUTotalTime;

        // 1. Pipeline Outer Setup Phase
        addSegment("GPU Cross-Check Engine",     stats.gpuCrossCheckEngineMs,            ImVec4(0.2f, 0.5f, 0.85f, 1.0f), "Initial BVH criss-cross intersection pass");
        addSegment("Dual Tree Expansion",       stats.dualTreeStepMs,                   ImVec4(0.1f, 0.8f, 0.6f,  1.0f), "Dual-tree traversal depth step");
        addSegment("Parallel DFS Descent (A&B)", stats.parallelDfsDescentBMs,            ImVec4(0.2f, 0.7f, 0.75f, 1.0f), "BFS/DFS primitive descent and reverse map creation");
        
        // 2. Loop Initialization
        addSegment("Preallocation Phase",        stats.loopTracker.preallocateTimeMs,    ImVec4(0.6f, 0.4f, 0.85f, 1.0f), "Batch chunk workspace allocation");
        
        // 3. Batched Chunk Loop Phases (In order of per-chunk execution)
        addSegment("Assembly Phase",            stats.loopTracker.assemblyPhaseMs,       ImVec4(0.9f, 0.55f, 0.2f, 1.0f), "Chunk buffer staging and batch index kernel");
        addSegment("AABB Candidate Execution",  stats.loopTracker.executionPhaseMs,      ImVec4(0.85f, 0.85f, 0.2f, 1.0f), "Fine AABB overlap counting and candidate generation");
        addSegment("GPU Float Predicates",      stats.loopTracker.fineEvaluationPhaseMs, ImVec4(0.3f, 0.85f, 0.3f, 1.0f), "GPU float-assisted exact predicates evaluation");
        addSegment("GPU Double Predicates",     stats.loopTracker.gpuDoublePredicatesMs, ImVec4(0.1f, 0.6f, 0.35f, 1.0f), "GPU double-precision predicates pass (Yellow)");
        addSegment("Thrust Compaction & D2H",   stats.loopTracker.DownloadAndClean,      ImVec4(0.7f, 0.4f, 0.7f, 1.0f), "Thrust copy_if compaction, D2H transfers & chunk frees");
        addSegment("CPU Narrow-Phase Compute",   stats.loopTracker.CPUPredicates,         ImVec4(0.85f, 0.3f, 0.55f, 1.0f), "Exact CPU CGAL TBB filtering (Orange)");

        // 4. Pipeline Finalization & Teardown
        addSegment("Loop Workspace Cleanup",    stats.loopTracker.cleanupTimeMs,         ImVec4(0.5f, 0.3f, 0.3f, 1.0f), "End-of-loop workspace deallocation and stream sync");
        addSegment("Explicit Cleanup Sync",      stats.finalCleanupSyncMs,               ImVec4(0.85f, 0.35f, 0.35f, 1.0f), "Device memory synchronization and driver teardown");

        // 5. Compute residual unaccounted time ("Grey Area")
        double trackedSum = 0.0;
        for (const auto& seg : m_segments) {
            trackedSum += seg.ms;
        }

        double unaccountedMs = std::max(0.0, m_totalTimeMs - trackedSum);
        if (unaccountedMs > 0.001) {
            m_segments.push_back({
                "Unaccounted Overhead (Grey Area)",
                unaccountedMs,
                ImVec4(0.5f, 0.5f, 0.5f, 1.0f), // Neutral Grey
                "Unmeasured driver overhead, device vector allocations, stream syncs, or kernel launch latencies"
            });
        }

        m_isVisible = true;
    }

    void show() { m_isVisible = true; }
    void hide() { m_isVisible = false; }
    bool isVisible() const { return m_isVisible; }

private:
    void addSegment(const std::string& name, double ms, ImVec4 color, const std::string& desc) {
        if (ms > 0.0001) {
            m_segments.push_back({name, ms, color, desc});
        }
    }

    void renderUI() {
        if (!m_isVisible) return;

        ImGui::SetNextWindowSize(ImVec2(680, 420), ImGuiCond_FirstUseEver);

        if (ImGui::Begin("Pipeline Execution Profiler", &m_isVisible, ImGuiWindowFlags_NoCollapse)) {
            
            // Pipeline Summary Header
            ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 1.0f), "Pipeline Timeline Overview");
            ImGui::Text("Total Execution Time (GPUTotalTime): %.3f ms", m_totalTimeMs);
            ImGui::Separator();

            if (m_totalTimeMs <= 0.0 || m_segments.empty()) {
                ImGui::TextDisabled("No timing profile data available. Execute 'compute' or 'testBVH' first.");
                ImGui::End();
                return;
            }

            // ----------------------------------------------------------------
            // 1. STACKED HORIZONTAL BAR CHART
            // ----------------------------------------------------------------
            ImGui::Text("Sub-Phase Distribution Bar:");
            ImVec2 canvasPos = ImGui::GetCursorScreenPos();
            ImVec2 canvasSize(ImGui::GetContentRegionAvail().x, 35.0f);

            ImDrawList* drawList = ImGui::GetWindowDrawList();

            // Background frame bounding rect
            drawList->AddRectFilled(canvasPos, ImVec2(canvasPos.x + canvasSize.x, canvasPos.y + canvasSize.y),
                                    IM_COL32(30, 30, 30, 255), 4.0f);

            float currentX = canvasPos.x;
            ImVec2 mousePos = ImGui::GetMousePos();
            const TimingSegment* hoveredSeg = nullptr;

            for (const auto& seg : m_segments) {
                float segWidth = static_cast<float>((seg.ms / m_totalTimeMs) * canvasSize.x);
                if (segWidth < 1.0f) segWidth = 1.0f; // Minimum 1px visual visibility

                ImVec2 pMin(currentX, canvasPos.y);
                ImVec2 pMax(std::min(currentX + segWidth, canvasPos.x + canvasSize.x), canvasPos.y + canvasSize.y);

                ImU32 col32 = ImGui::ColorConvertFloat4ToU32(seg.color);
                drawList->AddRectFilled(pMin, pMax, col32, 0.0f);
                drawList->AddRect(pMin, pMax, IM_COL32(20, 20, 20, 180), 0.0f); // Border line

                // Check mouse hover for tooltip
                if (mousePos.x >= pMin.x && mousePos.x <= pMax.x &&
                    mousePos.y >= pMin.y && mousePos.y <= pMax.y) {
                    hoveredSeg = &seg;
                }

                currentX += segWidth;
                if (currentX >= canvasPos.x + canvasSize.x) break;
            }

            // Dummy item to reserve space in ImGui layout
            ImGui::Dummy(canvasSize);

            // Render Hover Tooltip
            if (hoveredSeg) {
                double pct = (hoveredSeg->ms / m_totalTimeMs) * 100.0;
                ImGui::BeginTooltip();
                ImGui::TextColored(hoveredSeg->color, "%s", hoveredSeg->name.c_str());
                ImGui::Separator();
                ImGui::Text("Duration   : %.3f ms", hoveredSeg->ms);
                ImGui::Text("Percentage : %.2f%%", pct);
                ImGui::TextDisabled("%s", hoveredSeg->description.c_str());
                ImGui::EndTooltip();
            }

            ImGui::Spacing();
            ImGui::Separator();

            // ----------------------------------------------------------------
            // 2. DETAILED BREAKDOWN TABLE & LEGEND
            // ----------------------------------------------------------------
            ImGui::Text("Detailed Timing Breakdown:");

            static ImGuiTableFlags flags = ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg | 
                                           ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV;

            if (ImGui::BeginTable("TimingBreakdownTable", 4, flags)) {
                ImGui::TableSetupColumn("Color", ImGuiTableColumnFlags_WidthFixed, 45.0f);
                ImGui::TableSetupColumn("Sub-Phase Name", ImGuiTableColumnFlags_WidthStretch);
                ImGui::TableSetupColumn("Time (ms)", ImGuiTableColumnFlags_WidthFixed, 90.0f);
                ImGui::TableSetupColumn("Share (%)", ImGuiTableColumnFlags_WidthFixed, 75.0f);
                ImGui::TableHeadersRow();

                for (const auto& seg : m_segments) {
                    double pct = (seg.ms / m_totalTimeMs) * 100.0;

                    ImGui::TableNextRow();
                    
                    // Column 0: Color Box
                    ImGui::TableSetColumnIndex(0);
                    ImGui::ColorButton(("##col_" + seg.name).c_str(), seg.color, 
                                        ImGuiColorEditFlags_NoTooltip, ImVec2(30, 15));

                    // Column 1: Name + Tooltip
                    ImGui::TableSetColumnIndex(1);
                    ImGui::Text("%s", seg.name.c_str());
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("%s", seg.description.c_str());
                    }

                    // Column 2: Milliseconds
                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text("%.3f ms", seg.ms);

                    // Column 3: Percentage
                    ImGui::TableSetColumnIndex(3);
                    ImGui::Text("%.2f%%", pct);
                }

                ImGui::EndTable();
            }
        }
        ImGui::End();
    }
};

#endif // EXECUTION_TIMING_VISUALIZER_H