#ifndef MESH_TRIANGLE_DEGENERACY_VISUALIZER_H
#define MESH_TRIANGLE_DEGENERACY_VISUALIZER_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <iomanip>
#include <sstream>

// TBB Includes
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

// Polyscope & ImGui Includes
#include "polyscope/polyscope.h"
#include "imgui.h"

// Your Mesh Definition Header
#include "CgalDefinitions.h"

class MeshTriangleDegeneracyVisualizer {
public:
    struct HistogramData {
        std::string meshName;
        std::vector<float> binCounts;
        std::vector<double> binEdges;
        double minLength = 0.0; // Represents Min Quality Q
        double maxLength = 0.0; // Represents Max Quality Q
        double avgLength = 0.0; // Represents Avg Quality Q
        double totalLength = 0.0;
        size_t totalEdges = 0;  // Represents total face count
        bool isValid = false;
    };

    struct LocalStats {
        double minLen = std::numeric_limits<double>::infinity();
        double maxLen = -1.0;
        double totalLen = 0.0;
        size_t count = 0;
        bool valid = false;
    };

private:
    HistogramData m_dataA;
    HistogramData m_dataB;
    bool m_isVisible = false;
    int m_numBins = 25;      // Default 25 bins
    bool m_useLogScale = false;

    // Shared global range across both meshes for consistent X-axis comparison
    double m_globalMinLength = 0.0;
    double m_globalMaxLength = 0.0;

    // Cached references for dynamic UI re-binning
    const Mesh* m_cachedMeshA = nullptr;
    const Mesh* m_cachedMeshB = nullptr;

    // Helper to evaluate triangle quality Q in [0, 1]
    static inline double computeTriangleQuality(double dx0, double dy0, double dz0,
                                                 double dx1, double dy1, double dz1) {
        // Cross product for Area vector
        double cx = dy0 * dz1 - dz0 * dy1;
        double cy = dz0 * dx1 - dx0 * dz1;
        double cz = dx0 * dy1 - dy0 * dx1;
        double doubleArea = std::sqrt(cx * cx + cy * cy + cz * cz);

        // Edge 2 vector (p1 - p0)
        double dx2 = dx0 - dx1;
        double dy2 = dy0 - dy1;
        double dz2 = dz0 - dz1;

        double l0_sq = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;
        double l1_sq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
        double l2_sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

        double sumSqEdges = l0_sq + l1_sq + l2_sq;
        if (sumSqEdges <= 1e-15) return 0.0;

        // Q = (4 * sqrt(3) * Area) / (a^2 + b^2 + c^2)
        // Note: doubleArea = 2 * Area, so 2 * sqrt(3) * doubleArea
        constexpr double TWO_SQRT_3 = 3.464101615137754587;
        double q = (TWO_SQRT_3 * doubleArea) / sumSqEdges;

        return std::clamp(q, 0.0, 1.0);
    }

public:
    MeshTriangleDegeneracyVisualizer() = default;

    // Registers the render hook into Polyscope's main loop
    void init() {
        polyscope::state::userCallback = [this]() {
            this->renderUI();
        };
    }

    // Call this inside your CLI command handler
    void computeAndShow(const Mesh& meshA, const Mesh& meshB) {
        m_cachedMeshA = &meshA;
        m_cachedMeshB = &meshB;
        recompute();
        m_isVisible = true; // Opens/unhides the window
    }

    void recompute() {
        if (!m_cachedMeshA && !m_cachedMeshB) return;

        // Pre-pass lambda to compute degeneracy stats across a mesh
        auto computeStats = [](const Mesh* mesh) -> LocalStats {
            LocalStats s;
            if (!mesh) return s;
            size_t nFaces = num_faces(*mesh);
            if (nFaces == 0) return s;

            auto coords = mesh->points();

            struct ReductionPass {
                double minLength = std::numeric_limits<double>::infinity();
                double maxLength = -1.0;
                double totalLength = 0.0;
                size_t validEdges = 0;
            };

            ReductionPass pass = tbb::parallel_reduce(
                tbb::blocked_range<size_t>(0, nFaces), ReductionPass(),
                [&](const tbb::blocked_range<size_t>& r, ReductionPass local) {
                    for (size_t i = r.begin(); i != r.end(); ++i) {
                        typename Mesh::Face_index f(static_cast<typename Mesh::size_type>(i));
                        if (mesh->is_removed(f)) continue;

                        auto h0 = mesh->halfedge(f);
                        auto h1 = mesh->next(h0);
                        auto h2 = mesh->next(h1);

                        auto p0 = coords[mesh->target(h0)];
                        auto p1 = coords[mesh->target(h1)];
                        auto p2 = coords[mesh->target(h2)];

                        double dx0 = p1.x() - p0.x(), dy0 = p1.y() - p0.y(), dz0 = p1.z() - p0.z();
                        double dx1 = p2.x() - p0.x(), dy1 = p2.y() - p0.y(), dz1 = p2.z() - p0.z();

                        double q = computeTriangleQuality(dx0, dy0, dz0, dx1, dy1, dz1);

                        local.minLength = std::min(local.minLength, q);
                        local.maxLength = std::max(local.maxLength, q);
                        local.totalLength += q;
                        local.validEdges++;
                    }
                    return local;
                },
                [](ReductionPass a, const ReductionPass& b) {
                    a.minLength = std::min(a.minLength, b.minLength);
                    a.maxLength = std::max(a.maxLength, b.maxLength);
                    a.totalLength += b.totalLength;
                    a.validEdges += b.validEdges;
                    return a;
                }
            );

            s.minLen = pass.minLength;
            s.maxLen = pass.maxLength;
            s.totalLen = pass.totalLength;
            s.count = pass.validEdges;
            s.valid = (pass.validEdges > 0);
            return s;
        };

        LocalStats statsA = computeStats(m_cachedMeshA);
        LocalStats statsB = computeStats(m_cachedMeshB);

        m_globalMinLength = std::min(statsA.minLen, statsB.minLen);
        m_globalMaxLength = std::max(statsA.maxLen, statsB.maxLen);

        // Fallback bounds clamp
        if (m_globalMinLength >= m_globalMaxLength) {
            m_globalMaxLength = m_globalMinLength + 1e-4;
        }

        // Bin both meshes using exact same global range
        if (m_cachedMeshA) m_dataA = computeHistogram(*m_cachedMeshA, "Mesh A", statsA);
        if (m_cachedMeshB) m_dataB = computeHistogram(*m_cachedMeshB, "Mesh B", statsB);
    }

    void show() { m_isVisible = true; }
    void hide() { m_isVisible = false; }
    bool isVisible() const { return m_isVisible; }

private:
    HistogramData computeHistogram(const Mesh& mesh, const std::string& name, const LocalStats& stats) {
        HistogramData data;
        data.meshName = name;
        if (!stats.valid) return data;

        size_t bins = static_cast<size_t>(std::max(1, m_numBins));

        data.binCounts.assign(bins, 0.0f);
        data.binEdges.assign(bins + 1, 0.0);

        data.minLength = stats.minLen;
        data.maxLength = stats.maxLen;
        data.totalLength = stats.totalLen;
        data.avgLength = stats.totalLen / static_cast<double>(stats.count);
        data.totalEdges = stats.count;

        double globalDiff = m_globalMaxLength - m_globalMinLength;

        // Bounding floor for log scales when metrics approach zero
        double safeMinLen = std::max(m_globalMinLength, 1e-8);
        double logMin = std::log10(safeMinLen);
        double logMax = std::log10(std::max(m_globalMaxLength, safeMinLen + 1e-8));
        double logRange = logMax - logMin;

        if (m_useLogScale) {
            double logStep = logRange / static_cast<double>(bins);
            for (size_t i = 0; i <= bins; ++i) {
                data.binEdges[i] = std::pow(10.0, logMin + static_cast<double>(i) * logStep);
            }
        } else {
            double binWidth = globalDiff / static_cast<double>(bins);
            for (size_t i = 0; i <= bins; ++i) {
                data.binEdges[i] = m_globalMinLength + static_cast<double>(i) * binWidth;
            }
        }

        size_t nFaces = num_faces(mesh);
        auto coords = mesh.points();

        // Parallel Face Quality Bin Classification
        tbb::concurrent_vector<size_t> localBins(bins, 0);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, nFaces), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                typename Mesh::Face_index f(static_cast<typename Mesh::size_type>(i));
                if (mesh.is_removed(f)) continue;

                auto h0 = mesh.halfedge(f);
                auto h1 = mesh.next(h0);
                auto h2 = mesh.next(h1);

                auto p0 = coords[mesh.target(h0)];
                auto p1 = coords[mesh.target(h1)];
                auto p2 = coords[mesh.target(h2)];

                double dx0 = p1.x() - p0.x(), dy0 = p1.y() - p0.y(), dz0 = p1.z() - p0.z();
                double dx1 = p2.x() - p0.x(), dy1 = p2.y() - p0.y(), dz1 = p2.z() - p0.z();

                double q = computeTriangleQuality(dx0, dy0, dz0, dx1, dy1, dz1);

                int idx = 0;
                if (m_useLogScale) {
                    double val = std::log10(std::max(q, 1e-8));
                    idx = static_cast<int>((val - logMin) / (logRange / static_cast<double>(bins)));
                } else {
                    idx = static_cast<int>((q - m_globalMinLength) / (globalDiff / static_cast<double>(bins)));
                }

                idx = std::clamp(idx, 0, static_cast<int>(bins - 1));
                localBins[idx]++;
            }
        });

        for (size_t i = 0; i < bins; ++i) {
            data.binCounts[i] = static_cast<float>(localBins[i]);
        }

        data.isValid = true;
        return data;
    }

    void renderUI() {
        if (!m_isVisible) return;

        ImGui::SetNextWindowSize(ImVec2(580, 560), ImGuiCond_FirstUseEver);

        if (ImGui::Begin("Triangle Degeneracy Distribution", &m_isVisible, ImGuiWindowFlags_NoCollapse)) {

            // Shared Axis Overview
            ImGui::Text("Shared Quality (Q) Bounds (0=Degenerate, 1=Equilateral):");
            ImGui::TextDisabled("Global Min Q: %.4e | Global Max Q: %.4e", m_globalMinLength, m_globalMaxLength);
            
            ImGui::Separator();

            // Controls Toolbar
            ImGui::Text("Histogram Settings:");
            ImGui::SetNextItemWidth(160.0f);
            if (ImGui::SliderInt("Bins", &m_numBins, 5, 100)) {
                recompute();
            }
            ImGui::SameLine();
            if (ImGui::Checkbox("Log Scale", &m_useLogScale)) {
                recompute();
            }
            ImGui::SameLine();
            if (ImGui::Button("Reset")) {
                m_numBins = 25;
                m_useLogScale = false;
                recompute();
            }

            ImGui::Separator();

            renderMeshSection(m_dataA);
            ImGui::Separator();
            renderMeshSection(m_dataB);
        }
        ImGui::End();
    }

    void renderMeshSection(const HistogramData& d) {
        if (!d.isValid) {
            ImGui::Text("%s: No geometry data available.", d.meshName.c_str());
            return;
        }

        // Metrics Section Header
        ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 1.0f), "%s Metrics", d.meshName.c_str());
        ImGui::Text("Total Faces: %zu | Cumulative Q: %.4e", d.totalEdges, d.totalLength);
        ImGui::Text("Mesh Local Quality: Min Q: %.4e | Avg Q: %.4f | Max Q: %.4f", d.minLength, d.avgLength, d.maxLength);

        float maxCount = *std::max_element(d.binCounts.begin(), d.binCounts.end());

        ImGui::TextDisabled("Y-Axis Max: %.0f faces", maxCount);

        std::string label = "##Hist_Degeneracy_" + d.meshName;
        ImVec2 canvasSize(ImGui::GetContentRegionAvail().x, 100.0f);

        // Store current cursor position to overlay interaction rect
        ImVec2 startPos = ImGui::GetCursorScreenPos();

        // Render Histogram Visual
        ImGui::PlotHistogram(
            label.c_str(),
            d.binCounts.data(),
            static_cast<int>(d.binCounts.size()),
            0,
            nullptr,
            0.0f,
            maxCount * 1.05f,
            canvasSize
        );

        // Overlay Invisible Button for reliable hover & mouse calculation
        ImGui::SetCursorScreenPos(startPos);
        std::string btnId = "##HistBtn_" + d.meshName;
        ImGui::InvisibleButton(btnId.c_str(), canvasSize);

        // Interactive Hover Tooltip detailing exact bin range
        if (ImGui::IsItemHovered() && d.binCounts.size() > 0 && d.binEdges.size() > d.binCounts.size()) {
            ImVec2 mousePos = ImGui::GetMousePos();

            float hoverRelX = (mousePos.x - startPos.x) / canvasSize.x;
            hoverRelX = std::clamp(hoverRelX, 0.0f, 0.9999f);

            int binIdx = static_cast<int>(hoverRelX * static_cast<float>(d.binCounts.size()));
            binIdx = std::clamp(binIdx, 0, static_cast<int>(d.binCounts.size()) - 1);

            size_t count = static_cast<size_t>(d.binCounts[binIdx]);
            double pct = d.totalEdges > 0 ? (100.0 * count / static_cast<double>(d.totalEdges)) : 0.0;

            double binStart = d.binEdges[binIdx];
            double binEnd = d.binEdges[binIdx + 1];
            double binMid = 0.5 * (binStart + binEnd);

            ImGui::BeginTooltip();
            ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f), "[ %s ] Bin #%d of %zu", d.meshName.c_str(), binIdx + 1, d.binCounts.size());
            ImGui::Separator();
            ImGui::Text("Quality Metric (Q) Range:");
            ImGui::Text("  Min Bound : %.6e (%.4f)", binStart, binStart);
            ImGui::Text("  Max Bound : %.6e (%.4f)", binEnd, binEnd);
            ImGui::Text("  Midpoint  : %.6e (%.4f)", binMid, binMid);
            ImGui::Separator();
            ImGui::Text("Face Count  : %zu (%.2f%% of mesh)", count, pct);
            ImGui::EndTooltip();
        }

        // Axis Bounds Display below plot
        ImGui::Text("X-Axis Range: %.2e", m_globalMinLength);
        ImGui::SameLine(ImGui::GetContentRegionAvail().x - 110);
        ImGui::Text("%.2e (Quality)", m_globalMaxLength);
    }
};

#endif // MESH_TRIANGLE_DEGENERACY_VISUALIZER_H