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
    std::string m_windowTitle = "Pipeline Execution Profiler";
    double m_totalTimeMs = 0.0;
    std::vector<TimingSegment> m_segments;

public:
    ExecutionTimingVisualizer() = default;

    void init() {
        auto existingCallback = polyscope::state::userCallback;
        polyscope::state::userCallback = [this, existingCallback]() {
            if (existingCallback) existingCallback();
            this->renderUI();
        };
    }

    // --------------------------------------------------------------------
    // 1. DEDICATED BREAKDOWN FOR 'COMPUTE' (Runtime Traversal & Predicates)
    // --------------------------------------------------------------------
    void updateAndShowCompute(const ExecutionStats& stats) {
        m_segments.clear();
        m_windowTitle = "Pipeline Profiler: COMPUTE (Query Only)";
        m_totalTimeMs = stats.GPUTotalTime;

        // Query & Traversal Setup
        addSegment("GPU Cross-Check Engine",     stats.gpuCrossCheckEngineMs,            ImVec4(0.20f, 0.50f, 0.85f, 1.0f), "Initial BVH criss-cross intersection pass");
        addSegment("Dual Tree Expansion",       stats.dualTreeStepMs,                   ImVec4(0.10f, 0.80f, 0.60f, 1.0f), "Dual-tree traversal depth step");
        addSegment("Parallel DFS Descent (A&B)", stats.parallelDfsDescentBMs,            ImVec4(0.20f, 0.70f, 0.75f, 1.0f), "BFS/DFS primitive descent and reverse map creation");
        
        // Chunk Processing Loop
        addSegment("Preallocation Phase",        stats.loopTracker.preallocateTimeMs,    ImVec4(0.60f, 0.40f, 0.85f, 1.0f), "Batch chunk workspace allocation");
        addSegment("Assembly Phase",            stats.loopTracker.assemblyPhaseMs,       ImVec4(0.90f, 0.55f, 0.20f, 1.0f), "Chunk buffer staging and batch index kernel");
        addSegment("AABB Candidate Execution",  stats.loopTracker.executionPhaseMs,      ImVec4(0.85f, 0.85f, 0.20f, 1.0f), "Fine AABB overlap counting and candidate generation");
        addSegment("GPU Float Predicates",      stats.loopTracker.fineEvaluationPhaseMs, ImVec4(0.30f, 0.85f, 0.30f, 1.0f), "GPU float-assisted exact predicates evaluation");
        addSegment("GPU Double Predicates",     stats.loopTracker.gpuDoublePredicatesMs, ImVec4(0.10f, 0.60f, 0.35f, 1.0f), "GPU double-precision predicates pass (Yellow)");
        addSegment("Thrust Compaction & D2H",   stats.loopTracker.DownloadAndClean,      ImVec4(0.70f, 0.40f, 0.70f, 1.0f), "Thrust copy_if compaction, D2H transfers & chunk frees");
        addSegment("CPU Narrow-Phase Compute",   stats.loopTracker.CPUPredicates,         ImVec4(0.85f, 0.30f, 0.55f, 1.0f), "Exact CPU CGAL TBB filtering (Orange)");

        // Finalization
        addSegment("Loop Workspace Cleanup",    stats.loopTracker.cleanupTimeMs,         ImVec4(0.50f, 0.30f, 0.30f, 1.0f), "End-of-loop workspace deallocation and stream sync");
        addSegment("Explicit Cleanup Sync",      stats.finalCleanupSyncMs,               ImVec4(0.85f, 0.35f, 0.35f, 1.0f), "Device memory synchronization and driver teardown");

        computeUnaccountedOverhead();
        m_isVisible = true;
    }

    // --------------------------------------------------------------------
    // 2. DEDICATED BREAKDOWN FOR 'TESTBVH' (Allocations + Trees + Compute)
    // --------------------------------------------------------------------
    void updateAndShowTestBVH(const ExecutionStats& stats) {
        m_segments.clear();
        m_windowTitle = "Pipeline Profiler: TESTBVH (End-To-End)";
        m_totalTimeMs = stats.GPUTotalTime;

        // Phase 1: Allocations & Transfers
        addSegment("1. Raw Copy & Box Gen",      stats.initialAllocAndCopyMs,           ImVec4(0.95f, 0.40f, 0.20f, 1.0f), "Host-to-Device buffer copies and box generation");
        addSegment("2. Thrust Workspaces Init",  stats.thrustInitOverheadMs,            ImVec4(0.85f, 0.50f, 0.30f, 1.0f), "Thrust device vector allocations");

        // Phase 2: BVH Construction
        addSegment("3. Build/Refit Tree Mesh A", stats.buildRefitMeshAMs,               ImVec4(0.35f, 0.65f, 0.95f, 1.0f), "Linear init, Forest BVH expansion & refitting (Mesh A)");
        addSegment("4. Build/Refit Tree Mesh B", stats.buildRefitMeshBMs,               ImVec4(0.25f, 0.75f, 0.85f, 1.0f), "Linear init, Forest BVH expansion & refitting (Mesh B)");

        // Phase 3: Traversal & Cross Check
        addSegment("5. GPU Cross-Check Engine",  stats.gpuCrossCheckEngineMs,            ImVec4(0.20f, 0.50f, 0.85f, 1.0f), "Initial BVH criss-cross intersection pass");
        addSegment("6. Dual Tree Expansion",    stats.dualTreeStepMs,                   ImVec4(0.10f, 0.80f, 0.60f, 1.0f), "Dual-tree traversal depth step");
        addSegment("7. Parallel DFS Descent",   stats.parallelDfsDescentBMs,            ImVec4(0.20f, 0.70f, 0.75f, 1.0f), "BFS/DFS primitive descent and reverse map creation");

        // Phase 4: Chunk Execution Loop
        addSegment("8. Preallocation Phase",     stats.loopTracker.preallocateTimeMs,    ImVec4(0.60f, 0.40f, 0.85f, 1.0f), "Batch chunk workspace allocation");
        addSegment("9. Assembly Phase",         stats.loopTracker.assemblyPhaseMs,       ImVec4(0.90f, 0.55f, 0.20f, 1.0f), "Chunk buffer staging and batch index kernel");
        addSegment("10. AABB Candidate Exec",   stats.loopTracker.executionPhaseMs,      ImVec4(0.85f, 0.85f, 0.20f, 1.0f), "Fine AABB overlap counting and candidate generation");
        addSegment("11. GPU Float Predicates",   stats.loopTracker.fineEvaluationPhaseMs, ImVec4(0.30f, 0.85f, 0.30f, 1.0f), "GPU float-assisted exact predicates evaluation");
        addSegment("12. GPU Double Predicates",  stats.loopTracker.gpuDoublePredicatesMs, ImVec4(0.10f, 0.60f, 0.35f, 1.0f), "GPU double-precision predicates pass (Yellow)");
        addSegment("13. Thrust Compaction & D2H",stats.loopTracker.DownloadAndClean,      ImVec4(0.70f, 0.40f, 0.70f, 1.0f), "Thrust copy_if compaction, D2H transfers & chunk frees");
        addSegment("14. CPU Narrow-Phase",       stats.loopTracker.CPUPredicates,         ImVec4(0.85f, 0.30f, 0.55f, 1.0f), "Exact CPU CGAL TBB filtering (Orange)");

        // Phase 5: Cleanup
        addSegment("15. Loop Workspace Cleanup", stats.loopTracker.cleanupTimeMs,         ImVec4(0.50f, 0.30f, 0.30f, 1.0f), "End-of-loop workspace deallocation and stream sync");
        addSegment("16. Explicit Cleanup Sync",  stats.finalCleanupSyncMs,               ImVec4(0.85f, 0.35f, 0.35f, 1.0f), "Device memory synchronization and driver teardown");

        computeUnaccountedOverhead();
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

    void computeUnaccountedOverhead() {
        double trackedSum = 0.0;
        for (const auto& seg : m_segments) {
            trackedSum += seg.ms;
        }

        double unaccountedMs = std::max(0.0, m_totalTimeMs - trackedSum);
        if (unaccountedMs > 0.001) {
            m_segments.push_back({
                "Unaccounted Overhead (Grey Area)",
                unaccountedMs,
                ImVec4(0.5f, 0.5f, 0.5f, 1.0f),
                "Unmeasured driver overhead, memory allocations, or launch latencies"
            });
        }
    }

    void renderUI() {
        if (!m_isVisible) return;

        ImGui::SetNextWindowSize(ImVec2(700, 440), ImGuiCond_FirstUseEver);

        if (ImGui::Begin(m_windowTitle.c_str(), &m_isVisible, ImGuiWindowFlags_NoCollapse)) {
            ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 1.0f), "Pipeline Timeline Overview");
            ImGui::Text("Total Execution Time (GPUTotalTime): %.3f ms", m_totalTimeMs);
            ImGui::Separator();

            if (m_totalTimeMs <= 0.0 || m_segments.empty()) {
                ImGui::TextDisabled("No timing profile data available.");
                ImGui::End();
                return;
            }

            // Stacked Bar
            ImGui::Text("Sub-Phase Distribution Bar:");
            ImVec2 canvasPos = ImGui::GetCursorScreenPos();
            ImVec2 canvasSize(ImGui::GetContentRegionAvail().x, 35.0f);
            ImDrawList* drawList = ImGui::GetWindowDrawList();

            drawList->AddRectFilled(canvasPos, ImVec2(canvasPos.x + canvasSize.x, canvasPos.y + canvasSize.y),
                                    IM_COL32(30, 30, 30, 255), 4.0f);

            float currentX = canvasPos.x;
            ImVec2 mousePos = ImGui::GetMousePos();
            const TimingSegment* hoveredSeg = nullptr;

            for (const auto& seg : m_segments) {
                float segWidth = static_cast<float>((seg.ms / m_totalTimeMs) * canvasSize.x);
                if (segWidth < 1.0f) segWidth = 1.0f;

                ImVec2 pMin(currentX, canvasPos.y);
                ImVec2 pMax(std::min(currentX + segWidth, canvasPos.x + canvasSize.x), canvasPos.y + canvasSize.y);

                ImU32 col32 = ImGui::ColorConvertFloat4ToU32(seg.color);
                drawList->AddRectFilled(pMin, pMax, col32, 0.0f);
                drawList->AddRect(pMin, pMax, IM_COL32(20, 20, 20, 180), 0.0f);

                if (mousePos.x >= pMin.x && mousePos.x <= pMax.x &&
                    mousePos.y >= pMin.y && mousePos.y <= pMax.y) {
                    hoveredSeg = &seg;
                }

                currentX += segWidth;
                if (currentX >= canvasPos.x + canvasSize.x) break;
            }

            ImGui::Dummy(canvasSize);

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

            // Table Breakdown
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
                    ImGui::TableSetColumnIndex(0);
                    ImGui::ColorButton(("##col_" + seg.name).c_str(), seg.color, 
                                        ImGuiColorEditFlags_NoTooltip, ImVec2(30, 15));

                    ImGui::TableSetColumnIndex(1);
                    ImGui::Text("%s", seg.name.c_str());
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("%s", seg.description.c_str());
                    }

                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text("%.3f ms", seg.ms);

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