#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace cuBQL {
  namespace samples {

    std::vector<Triangle> loadOFF(const std::string &offFile)
    {
        std::ifstream in(offFile);
        if (!in.is_open()) {
            throw std::runtime_error("Could not open OFF file: " + offFile);
        }

        std::string header;
        in >> header;
        if (header != "OFF") {
            throw std::runtime_error("Not a valid OFF file: " + offFile);
        }

        size_t numVertices, numFaces, numEdges;
        in >> numVertices >> numFaces >> numEdges;

        std::vector<vec3f> vertices(numVertices);
        for (size_t i = 0; i < numVertices; ++i) {
            in >> vertices[i].x >> vertices[i].y >> vertices[i].z;
        }

        std::vector<Triangle> triangles;
        for (size_t i = 0; i < numFaces; ++i) {
            int numV;
            in >> numV;
            
            // OFF supports polygons; this assumes triangulation
            if (numV == 3) {
                int i0, i1, i2;
                in >> i0 >> i1 >> i2;
                triangles.push_back({vertices[i0], vertices[i1], vertices[i2]});
            } else if (numV > 3) {
                // Basic fan triangulation for polygons
                std::vector<int> polyIndices(numV);
                for (int j = 0; j < numV; ++j) in >> polyIndices[j];
                for (int j = 1; j < numV - 1; ++j) {
                    triangles.push_back({vertices[polyIndices[0]], 
                                         vertices[polyIndices[j]], 
                                         vertices[polyIndices[j+1]]});
                }
            }
        }

        return triangles;
    }

  } // ::cuBQL::samples
} // ::cuBQL
