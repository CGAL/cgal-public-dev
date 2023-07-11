// Copyright (c) 2023 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Tolga Talha Yildiz  <tolgayildiz@tamu.edu>)

#ifndef POLYHEDRON_RENDERER
#define POLYHEDRON_RENDERER

#include <CGAL/Surface_mesh.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#include <CGAL/draw_face_graph.h>
#include <CGAL/Random.h>
#include <QVulkanWindow>
#include <QVulkanFunctions>
#include <QFile>


#include <CGAL/Dynamic_property_map.h>
#include <CGAL/boost/graph/helpers.h>

namespace glm {
    struct vec3 {
        float x;
        float y;
        float z;
    };
    struct vec2 {
        float x;
        float y;
    };
    struct mat4 {
        float a_11;
        float a_12;
        float a_13;
        float a_14;
        float a_21;
        float a_22;
        float a_23;
        float a_24;
        float a_31;
        float a_32;
        float a_33;
        float a_34;
        float a_41;
        float a_42;
        float a_43;
        float a_44;
    };
}

struct Vertex {
    glm::vec3 pos;
    glm::vec3 color;

    bool operator==(const Vertex& other) const {
        return pos.x == other.pos.x && pos.y == other.pos.y && pos.z == other.pos.z;
    }
};

inline void hash_combine(size_t& seed, size_t hash)
{
    hash += 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash;
}

struct UniformBufferObject {
    alignas(16) glm::mat4 model;
};

namespace std {
    template<> struct hash<Vertex> {
        size_t operator()(Vertex const& vertex) const {
            size_t seed1 = 0;
            hash_combine(seed1, vertex.pos.x);
            hash_combine(seed1, vertex.pos.y);
            hash_combine(seed1, vertex.pos.z);
            return seed1;
        }
    };
}

namespace CGAL {
#define CGAL_POLY_TYPE Surface_mesh<K>

    static const int UNIFORM_DATA_SIZE = 16 * sizeof(float);

    static inline VkDeviceSize aligned(VkDeviceSize v, VkDeviceSize byteAlign)
    {
        return (v + byteAlign - 1) & ~(byteAlign - 1);
    }

    template <typename K>
        class Polyhedron_renderer : public QVulkanWindowRenderer {
        public:
            Polyhedron_renderer(QVulkanWindow* w, CGAL_POLY_TYPE& poly) : m_window(w), m_buffVertex(nullptr), m_buffVertexMem(nullptr), m_descPool(nullptr), m_descSetLayout(nullptr), m_pipeline(nullptr), m_devFuncs(nullptr), m_pipelineCache(nullptr), m_pipelineLayout(nullptr), m_proj(QMatrix4x4()), m_rotation(0.0f), m_descSet(), m_uniformBufferInfo(), m_poly(poly) {}

            void initResources() override {
                m_window->availablePhysicalDevices();
                std::vector<float> vData{};
                std::unordered_map<Vertex, uint32_t> uniqueVertices{};
                for (typename CGAL_POLY_TYPE::Face_iterator i = m_poly.facets_begin(); i != m_poly.facets_end(); i++) {
                    typename CGAL_POLY_TYPE::Halfedge_around_facet_circulator j = i->facet_begin();
                    typename CGAL_POLY_TYPE::Halfedge_around_facet_circulator k = j;
                    k++;
                    do {
                        Vertex v1{ {j->vertex()->point().x(),j->vertex()->point().y(),j->vertex()->point().z()},{0.5f,0.5f,1.0f } };
                        if (uniqueVertices.count(v1) == 0) {
                            uniqueVertices[v1] = static_cast<uint32_t>(vertices.size());
                            vertices.push_back(v1);
                        }
                        indices.push_back(uniqueVertices[v1]);

                        Vertex v2{ {k->vertex()->point().x(),k->vertex()->point().y(),k->vertex()->point().z()},{0.5f,0.5f,1.0f } };
                        if (uniqueVertices.count(v2) == 0) {
                            uniqueVertices[v2] = static_cast<uint32_t>(vertices.size());
                            vertices.push_back(v2);
                        }
                        indices.push_back(uniqueVertices[v2]);
                        k++;
                        Vertex v3{ {k->vertex()->point().x(),k->vertex()->point().y(),k->vertex()->point().z()},{0.5f,0.5f,1.0f } };
                        if (uniqueVertices.count(v3) == 0) {
                            uniqueVertices[v3] = static_cast<uint32_t>(vertices.size());
                            vertices.push_back(v3);
                        }
                        indices.push_back(uniqueVertices[v3]);
                    } while (k != i->facet_begin());
                }
                for (typename CGAL_POLY_TYPE::Face_iterator i = m_poly.facets_begin(); i != m_poly.facets_end(); i++) {
                    typename CGAL_POLY_TYPE::Halfedge_around_facet_circulator j = i->facet_begin();
                    do {
                        Vertex v1{ {j->vertex()->point().x(),j->vertex()->point().y(),j->vertex()->point().z()},{0.0f,0.0f,0.0f} };
                        indicesGraph.push_back(uniqueVertices[v1]);
                        j++;
                        Vertex v2{ {j->vertex()->point().x(),j->vertex()->point().y(),j->vertex()->point().z()},{0.0f,0.0f,0.0f} };
                        indicesGraph.push_back(uniqueVertices[v2]);
                    } while (j != i->facet_begin());
                }
                for (int i = 0; i < indices.size(); i += 3) {
                    glm::vec3 v1 = vertices[indices[i]].pos;
                    glm::vec3 v2 = vertices[indices[i + 1]].pos;
                    glm::vec3 v3 = vertices[indices[i + 2]].pos;
                    //printf("Triangle %i: \n\t v1: %f, %f, %f\n\t v2: %f, %f, %f\n\t v3: %f, %f, %f\n",i/3, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z);
                }
                /*
                for (typename CGAL_POLY_TYPE::Face_iterator i = m_poly.facets_begin(); i != m_poly.facets_end(); i++) {
                    typename CGAL_POLY_TYPE::Halfedge_around_facet_circulator j = i->facet_begin();
                    typename CGAL_POLY_TYPE::Halfedge_around_facet_circulator k = j;
                    k++;
                    do {

                        //printf("%f, %f, %f, 1.0f, 1.0f, 1.0f,\n", j->vertex()->point().x(), j->vertex()->point().y(), j->vertex()->point().z());
                        vData.push_back(j->vertex()->point().x());
                        vData.push_back(j->vertex()->point().y());
                        vData.push_back(j->vertex()->point().z());
                        vData.push_back(1.0f);
                        vData.push_back(1.0f);
                        vData.push_back(1.0f);
                        printf("%f, %f, %f, 1.0f, 1.0f, 1.0f,\n", k->vertex()->point().x(), k->vertex()->point().y(), k->vertex()->point().z());
                        vData.push_back(k->vertex()->point().x());
                        vData.push_back(k->vertex()->point().y());
                        vData.push_back(k->vertex()->point().z());
                        vData.push_back(1.0f);
                        vData.push_back(1.0f);
                        vData.push_back(1.0f);
                        ++k;
                        printf("%f, %f, %f, 1.0f, 1.0f, 1.0f,\n", k->vertex()->point().x(), k->vertex()->point().y(), k->vertex()->point().z());
                        vData.push_back(k->vertex()->point().x());
                        vData.push_back(k->vertex()->point().y());
                        vData.push_back(k->vertex()->point().z());
                        vData.push_back(1.0f);
                        vData.push_back(1.0f);
                        vData.push_back(1.0f);
                    } while (k != i->facet_begin());
                }
                */
                //for (auto f : faces(m_poly)) {
                //    I_HalfedgeDS_facet_circ<Iterator_from_circulator,  v_start = f->facet_begin();
                //    //Halfedge_around_face_circulator<PolyhedronTraits_3> v_start = f->facet_begin();
                //    I_HalfedgeDS_facet_circ v = v_start;
                //    do
                //    {
                //        printf("a possible vertex in a face : %f, %f, %f \n", v->point().x(), v->point().y(), v->point().z());
                //    } while (++v != v_start);
                //}
                //for (auto v : vertices(m_poly))
                //{
                //    printf("a possible vertex position: %f, %f, %f \n", v->point().x(), v->point().y(), v->point().z());
                //    vData.push_back(v->point().x() / 5.0f);
                //    vData.push_back(v->point().y() / 5.0f);
                //    vData.push_back(v->point().z() / 5.0f);
                //    vData.push_back(1.0f);
                //    vData.push_back(1.0f);
                //    vData.push_back(1.0f);
                //}

                qDebug("initResources call");

                VkDevice dev = m_window->device();
                m_devFuncs = m_window->vulkanInstance()->deviceFunctions(dev);

                const int concFrameCount = m_window->concurrentFrameCount();
                const VkPhysicalDeviceLimits* pDevLimits = &m_window->physicalDeviceProperties()->limits;
                const VkDeviceSize uniAlign = pDevLimits->minUniformBufferOffsetAlignment;
                qDebug("uniform buffer offset alignment is %u", (uint)uniAlign);

                VkBufferCreateInfo buffInfo{};
                buffInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;

                const VkDeviceSize vertexAllocSize = aligned((size_t)sizeof(vertices) * vertices.size(), uniAlign);
                const VkDeviceSize uniformAllocSize = aligned(UNIFORM_DATA_SIZE, uniAlign);
                buffInfo.size = vertexAllocSize + concFrameCount * uniformAllocSize;
                buffInfo.usage = VK_BUFFER_USAGE_VERTEX_BUFFER_BIT | VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT;

                VkResult err = m_devFuncs->vkCreateBuffer(dev, &buffInfo, nullptr, &m_buffVertex);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create buffer: %d", err);

                VkMemoryRequirements memReq;
                m_devFuncs->vkGetBufferMemoryRequirements(dev, m_buffVertex, &memReq);

                VkMemoryAllocateInfo memAllocInfo{};
                memAllocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
                memAllocInfo.allocationSize = memReq.size;
                memAllocInfo.memoryTypeIndex = m_window->hostVisibleMemoryIndex();

                err = m_devFuncs->vkAllocateMemory(dev, &memAllocInfo, nullptr, &m_buffVertexMem);
                if (err != VK_SUCCESS)
                    qFatal("Failed to allocate memory: %d", err);

                err = m_devFuncs->vkBindBufferMemory(dev, m_buffVertex, m_buffVertexMem, 0);
                if (err != VK_SUCCESS)
                    qFatal("Failed to bind buffer memory: %d", err);

                quint8* p;
                err = m_devFuncs->vkMapMemory(dev, m_buffVertexMem, 0, memReq.size, 0, reinterpret_cast<void**>(&p));
                if (err != VK_SUCCESS)
                    qFatal("Failed to map memory: %d", err);

                memcpy(p, vertices.data(), (size_t)(sizeof(vertices[0]) * vertices.size()));
                QMatrix4x4 ident;
                memset(m_uniformBufferInfo, 0, sizeof(m_uniformBufferInfo));
                for (int i = 0; i < concFrameCount; i++) {
                    const VkDeviceSize offset = vertexAllocSize + i * uniformAllocSize;
                    memcpy(p + offset, ident.constData(), 16 * sizeof(float));
                    m_uniformBufferInfo[i].buffer = m_buffVertex;
                    m_uniformBufferInfo[i].offset = offset;
                    m_uniformBufferInfo[i].range = uniformAllocSize;
                }
                m_devFuncs->vkUnmapMemory(dev, m_buffVertexMem);
                createIndexBuffers(m_buffIndex, m_buffIndexMem, indices);
                createIndexBuffers(m_buffIndexGraph, m_buffIndexGraphMem, indicesGraph);

                VkVertexInputBindingDescription vertexBindingDesc{};
                vertexBindingDesc.binding = 0;
                vertexBindingDesc.stride = 6 * sizeof(float);
                vertexBindingDesc.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

                VkVertexInputAttributeDescription vertexAttrDesc[] = {
                    {//pos
                        0,//location
                        0,//binding
                        VK_FORMAT_R32G32B32_SFLOAT,
                        0
                    },
                    {//col
                        1,
                        0,
                        VK_FORMAT_R32G32B32_SFLOAT,
                        3 * sizeof(float)
                    }
                };

                VkPipelineVertexInputStateCreateInfo vertexInputInfo{};
                vertexInputInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;
                vertexInputInfo.pNext = nullptr;
                vertexInputInfo.flags = 0;
                vertexInputInfo.vertexBindingDescriptionCount = 1;
                vertexInputInfo.pVertexBindingDescriptions = &vertexBindingDesc;
                vertexInputInfo.vertexAttributeDescriptionCount = 2;
                vertexInputInfo.pVertexAttributeDescriptions = vertexAttrDesc;

                VkDescriptorPoolSize descPoolSizes{};
                descPoolSizes.type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
                descPoolSizes.descriptorCount = uint32_t(concFrameCount);

                VkDescriptorPoolCreateInfo descPoolInfo{};
                descPoolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
                descPoolInfo.maxSets = concFrameCount;
                descPoolInfo.poolSizeCount = 1;
                descPoolInfo.pPoolSizes = &descPoolSizes;

                err = m_devFuncs->vkCreateDescriptorPool(dev, &descPoolInfo, nullptr, &m_descPool);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create descriptor pool: %d", err);

                VkDescriptorSetLayoutBinding layoutBinding{};
                layoutBinding.binding = 0;
                layoutBinding.descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
                layoutBinding.descriptorCount = 1;
                layoutBinding.stageFlags = VK_SHADER_STAGE_VERTEX_BIT;

                VkDescriptorSetLayoutCreateInfo descSetLayoutInfo{};
                descSetLayoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
                descSetLayoutInfo.bindingCount = 1;
                descSetLayoutInfo.pBindings = &layoutBinding;

                err = m_devFuncs->vkCreateDescriptorSetLayout(dev, &descSetLayoutInfo, nullptr, &m_descSetLayout);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create descriptor set layout: %d", err);

                for (int i = 0; i < concFrameCount; i++) {
                    VkDescriptorSetAllocateInfo descSetAllocateInfo{};
                    descSetAllocateInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
                    descSetAllocateInfo.descriptorPool = m_descPool;
                    descSetAllocateInfo.pSetLayouts = &m_descSetLayout;
                    descSetAllocateInfo.descriptorSetCount = 1;

                    err = m_devFuncs->vkAllocateDescriptorSets(dev, &descSetAllocateInfo, &m_descSet[i]);
                    if (err != VK_SUCCESS)
                        qFatal("Failed to allocate descriptor set %i: %d, %i", i, err);

                    VkWriteDescriptorSet descWrite{};
                    descWrite.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
                    descWrite.dstSet = m_descSet[i];
                    descWrite.descriptorCount = 1;
                    descWrite.descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
                    descWrite.pBufferInfo = &m_uniformBufferInfo[i];
                    m_devFuncs->vkUpdateDescriptorSets(dev, 1, &descWrite, 0, nullptr);
                }

                VkPipelineCacheCreateInfo pipelineCacheInfo{};
                pipelineCacheInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_CACHE_CREATE_INFO;
                err = m_devFuncs->vkCreatePipelineCache(dev, &pipelineCacheInfo, nullptr, &m_pipelineCache);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create pipeline cache: %d", err);

                VkPipelineLayoutCreateInfo pipelineLayoutInfo{};
                pipelineLayoutInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
                pipelineLayoutInfo.setLayoutCount = 1;
                pipelineLayoutInfo.pSetLayouts = &m_descSetLayout;
                err = m_devFuncs->vkCreatePipelineLayout(dev, &pipelineLayoutInfo, nullptr, &m_pipelineLayout);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create pipeline layout: %d", err);

                VkShaderModule vertShader = createShader("../resources/color_vert.spv");
                VkShaderModule fragShader = createShader("../resources/color_frag.spv");
                VkShaderModule vertShaderWire = createShader("../resources/color_vert_wire.spv");
                VkShaderModule fragShaderWire = createShader("../resources/color_frag_wire.spv");
                VkShaderModule vertShaderPoint = createShader("../resources/color_vert_point.spv");
                VkShaderModule fragShaderPoint = createShader("../resources/color_frag_point.spv");

                VkPipelineShaderStageCreateInfo shaderStages[2] = {
                    {
                        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
                        nullptr,
                        0,
                        VK_SHADER_STAGE_VERTEX_BIT,
                        vertShader,
                        "main",
                        nullptr
                    },
                    {
                        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
                        nullptr,
                        0,
                        VK_SHADER_STAGE_FRAGMENT_BIT,
                        fragShader,
                        "main",
                        nullptr
                    }
                };

                VkPipelineShaderStageCreateInfo wireShaderStages[2] = {
                    {
                        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
                        nullptr,
                        0,
                        VK_SHADER_STAGE_VERTEX_BIT,
                        vertShaderWire,
                        "main",
                        nullptr
                    },
                    {
                        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
                        nullptr,
                        0,
                        VK_SHADER_STAGE_FRAGMENT_BIT,
                        fragShaderWire,
                        "main",
                        nullptr
                    }
                };

                VkPipelineShaderStageCreateInfo pointShaderStages[2] = {
                    {
                        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
                        nullptr,
                        0,
                        VK_SHADER_STAGE_VERTEX_BIT,
                        vertShaderPoint,
                        "main",
                        nullptr
                    },
                    {
                        VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
                        nullptr,
                        0,
                        VK_SHADER_STAGE_FRAGMENT_BIT,
                        fragShaderPoint,
                        "main",
                        nullptr
                    }
                };

                VkPipelineInputAssemblyStateCreateInfo iaStateInfo{};
                iaStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
                iaStateInfo.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

                VkPipelineInputAssemblyStateCreateInfo iaStateInfoEdges{};
                iaStateInfoEdges.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
                iaStateInfoEdges.topology = VK_PRIMITIVE_TOPOLOGY_LINE_LIST;

                VkPipelineViewportStateCreateInfo vpStateInfo{};
                vpStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
                vpStateInfo.viewportCount = 1;
                vpStateInfo.scissorCount = 1;

                VkPipelineRasterizationStateCreateInfo rastStateInfo{};
                rastStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
                rastStateInfo.polygonMode = VK_POLYGON_MODE_FILL;
                rastStateInfo.cullMode = VK_CULL_MODE_NONE;
                rastStateInfo.frontFace = VK_FRONT_FACE_CLOCKWISE;
                rastStateInfo.lineWidth = 1.0f;

                VkPipelineRasterizationStateCreateInfo rastStateWireInfo{};
                rastStateWireInfo = rastStateInfo;
                rastStateWireInfo.polygonMode = VK_POLYGON_MODE_FILL;
                rastStateWireInfo.lineWidth = 3.0f;

                VkPipelineRasterizationStateCreateInfo rastStatePointInfo{};
                rastStatePointInfo = rastStateInfo;
                rastStatePointInfo.polygonMode = VK_POLYGON_MODE_POINT;
                rastStatePointInfo.lineWidth = 10.0f;

                VkPipelineMultisampleStateCreateInfo msStateInfo{};
                msStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
                msStateInfo.rasterizationSamples = m_window->sampleCountFlagBits();

                VkPipelineDepthStencilStateCreateInfo dsStateInfo{};
                dsStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO;
                dsStateInfo.depthTestEnable = VK_TRUE;
                dsStateInfo.depthWriteEnable = VK_TRUE;
                dsStateInfo.depthCompareOp = VK_COMPARE_OP_LESS_OR_EQUAL;

                VkPipelineColorBlendAttachmentState cbAttachState{};
                cbAttachState.colorWriteMask = 0xF;

                VkPipelineColorBlendStateCreateInfo cbStateInfo{};
                cbStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
                cbStateInfo.attachmentCount = 1;
                cbStateInfo.pAttachments = &cbAttachState;

                VkDynamicState dyn[] = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR, VK_DYNAMIC_STATE_LINE_WIDTH };
                VkPipelineDynamicStateCreateInfo dynStateInfo{};
                dynStateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO;
                dynStateInfo.dynamicStateCount = 3;
                dynStateInfo.pDynamicStates = dyn;

                VkGraphicsPipelineCreateInfo pipelineInfo{};
                pipelineInfo.sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
                pipelineInfo.stageCount = 2;
                pipelineInfo.pStages = shaderStages;
                pipelineInfo.pVertexInputState = &vertexInputInfo;
                pipelineInfo.pInputAssemblyState = &iaStateInfo;
                pipelineInfo.pViewportState = &vpStateInfo;
                pipelineInfo.pRasterizationState = &rastStateInfo;
                pipelineInfo.pMultisampleState = &msStateInfo;
                pipelineInfo.pDepthStencilState = &dsStateInfo;
                pipelineInfo.pColorBlendState = &cbStateInfo;
                pipelineInfo.pDynamicState = &dynStateInfo;
                pipelineInfo.layout = m_pipelineLayout;
                pipelineInfo.renderPass = m_window->defaultRenderPass();

                err = m_devFuncs->vkCreateGraphicsPipelines(dev, m_pipelineCache, 1, &pipelineInfo, nullptr, &m_pipeline);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create graphics pipeline: %d", err);

                pipelineInfo.pInputAssemblyState = &iaStateInfoEdges;
                pipelineInfo.pRasterizationState = &rastStateWireInfo;
                pipelineInfo.pStages = wireShaderStages;

                err = m_devFuncs->vkCreateGraphicsPipelines(dev, m_pipelineCache, 1, &pipelineInfo, nullptr, &m_pipeline_wire);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create graphics pipeline: %d", err);

                pipelineInfo.pInputAssemblyState = &iaStateInfo;
                pipelineInfo.pRasterizationState = &rastStatePointInfo;
                pipelineInfo.pStages = pointShaderStages;

                err = m_devFuncs->vkCreateGraphicsPipelines(dev, m_pipelineCache, 1, &pipelineInfo, nullptr, &m_pipeline_point);
                if (err != VK_SUCCESS)
                    qFatal("Failed to create graphics pipeline: %d", err);

                if (vertShader)
                    m_devFuncs->vkDestroyShaderModule(dev, vertShader, nullptr);
                if (fragShader)
                    m_devFuncs->vkDestroyShaderModule(dev, fragShader, nullptr);
                if (vertShaderWire)
                    m_devFuncs->vkDestroyShaderModule(dev, vertShaderWire, nullptr);
                if (fragShaderWire)
                    m_devFuncs->vkDestroyShaderModule(dev, fragShaderWire, nullptr);
                if (vertShaderPoint)
                    m_devFuncs->vkDestroyShaderModule(dev, vertShaderPoint, nullptr);
                if (fragShaderPoint)
                    m_devFuncs->vkDestroyShaderModule(dev, fragShaderPoint, nullptr);
            }

            void initSwapChainResources() override {
                qDebug("initSwapChainResources");

                m_proj = m_window->clipCorrectionMatrix();
                const QSize sz = m_window->swapChainImageSize();
                m_proj.perspective(45.0f, sz.width() / (float)sz.height(), 0.01f, 100.0f);
                m_proj.translate(0, 0, -20);
            }

            void releaseSwapChainResources() override {
                qDebug("releaseSwapChainResources");
            }

            void releaseResources() override {
                qDebug("releaseResources");

                VkDevice dev = m_window->device();

                if (m_pipeline) {
                    m_devFuncs->vkDestroyPipeline(dev, m_pipeline, nullptr);
                    m_pipeline = VK_NULL_HANDLE;
                }

                if (m_pipeline_wire) {
                    m_devFuncs->vkDestroyPipeline(dev, m_pipeline_wire, nullptr);
                    m_pipeline_wire = VK_NULL_HANDLE;
                }

                if (m_pipeline_point) {
                    m_devFuncs->vkDestroyPipeline(dev, m_pipeline_point, nullptr);
                    m_pipeline_point = VK_NULL_HANDLE;
                }

                if (m_pipelineLayout) {
                    m_devFuncs->vkDestroyPipelineLayout(dev, m_pipelineLayout, nullptr);
                    m_pipelineLayout = VK_NULL_HANDLE;
                }

                if (m_pipelineCache) {
                    m_devFuncs->vkDestroyPipelineCache(dev, m_pipelineCache, nullptr);
                    m_pipelineCache = VK_NULL_HANDLE;
                }

                if (m_descSetLayout) {
                    m_devFuncs->vkDestroyDescriptorSetLayout(dev, m_descSetLayout, nullptr);
                    m_descSetLayout = VK_NULL_HANDLE;
                }

                if (m_descPool) {
                    m_devFuncs->vkDestroyDescriptorPool(dev, m_descPool, nullptr);
                    m_descPool = VK_NULL_HANDLE;
                }

                if (m_buffVertex) {
                    m_devFuncs->vkDestroyBuffer(dev, m_buffVertex, nullptr);
                    m_buffVertex = VK_NULL_HANDLE;
                }

                if (m_buffVertexMem) {
                    m_devFuncs->vkFreeMemory(dev, m_buffVertexMem, nullptr);
                    m_buffVertexMem = VK_NULL_HANDLE;
                }

                if (m_buffIndex) {
                    m_devFuncs->vkDestroyBuffer(dev, m_buffIndex, nullptr);
                    m_buffIndex = VK_NULL_HANDLE;
                }

                if (m_buffIndexMem) {
                    m_devFuncs->vkFreeMemory(dev, m_buffIndexMem, nullptr);
                    m_buffIndexMem = VK_NULL_HANDLE;
                }

                if (m_buffIndexGraph) {
                    m_devFuncs->vkDestroyBuffer(dev, m_buffIndexGraph, nullptr);
                    m_buffIndexGraph = VK_NULL_HANDLE;
                }

                if (m_buffIndexGraphMem) {
                    m_devFuncs->vkFreeMemory(dev, m_buffIndexGraphMem, nullptr);
                    m_buffIndexGraphMem = VK_NULL_HANDLE;
                }
            }

            void createIndexBuffers(VkBuffer& indexBuffer, VkDeviceMemory& indexBufferMemory, std::vector<uint32_t>& indices) {
                VkDeviceSize bufferSize = sizeof(indices[0]) * indices.size();

                VkBuffer stagingBuffer;
                VkDeviceMemory stagingBufferMemory;

                createBuffer(bufferSize, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, stagingBuffer, stagingBufferMemory);

                void* data;
                m_devFuncs->vkMapMemory(m_window->device(), stagingBufferMemory, 0, bufferSize, 0, &data);
                memcpy(data, indices.data(), (size_t)bufferSize);
                m_devFuncs->vkUnmapMemory(m_window->device(), stagingBufferMemory);

                createBuffer(bufferSize, VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_INDEX_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, indexBuffer, indexBufferMemory);
                copyBuffer(stagingBuffer, indexBuffer, bufferSize);

                m_devFuncs->vkDestroyBuffer(m_window->device(), stagingBuffer, nullptr);
                m_devFuncs->vkFreeMemory(m_window->device(), stagingBufferMemory, nullptr);

            }

            void createVertexBuffers(VkBuffer& vertexBuffer, VkDeviceMemory& vertexBufferMemory, std::vector<Vertex>& vertices) {
                VkDeviceSize bufferSize = sizeof(vertices[0]) * vertices.size();

                VkBuffer stagingBuffer;
                VkDeviceMemory stagingBufferMemory;

                createBuffer(bufferSize, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, stagingBuffer, stagingBufferMemory);

                void* data;
                m_devFuncs->vkMapMemory(m_window->device(), stagingBufferMemory, 0, bufferSize, 0, &data);
                memcpy(data, vertices.data(), (size_t)bufferSize);
                m_devFuncs->vkUnmapMemory(m_window->device(), stagingBufferMemory);

                createBuffer(bufferSize, VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, vertexBuffer, vertexBufferMemory);
                copyBuffer(stagingBuffer, vertexBuffer, bufferSize);

                m_devFuncs->vkDestroyBuffer(m_window->device(), stagingBuffer, nullptr);
                m_devFuncs->vkFreeMemory(m_window->device(), stagingBufferMemory, nullptr);
            }

            void copyBuffer(VkBuffer src, VkBuffer dst, VkDeviceSize size) {
                VkCommandBuffer commandBuffer = beginSingleTimeCommands();

                VkBufferCopy copyRegion{};
                copyRegion.srcOffset = 0;
                copyRegion.dstOffset = 0;
                copyRegion.size = size;
                m_devFuncs->vkCmdCopyBuffer(commandBuffer, src, dst, 1, &copyRegion);

                endSingleTimeCommands(commandBuffer);
            }

            VkCommandBuffer beginSingleTimeCommands() {
                VkCommandBufferAllocateInfo allocInfo{};
                allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
                allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
                allocInfo.commandPool = m_window->graphicsCommandPool();
                allocInfo.commandBufferCount = 1;

                VkCommandBuffer commandBuffer;
                m_devFuncs->vkAllocateCommandBuffers(m_window->device(), &allocInfo, &commandBuffer);

                VkCommandBufferBeginInfo beginInfo{};
                beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
                beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

                m_devFuncs->vkBeginCommandBuffer(commandBuffer, &beginInfo);

                return commandBuffer;
            }

            void endSingleTimeCommands(VkCommandBuffer commandBuffer) {
                m_devFuncs->vkEndCommandBuffer(commandBuffer);
                VkSubmitInfo submitInfo{};
                submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
                submitInfo.commandBufferCount = 1;
                submitInfo.pCommandBuffers = &commandBuffer;

                m_devFuncs->vkQueueSubmit(m_window->graphicsQueue(), 1, &submitInfo, VK_NULL_HANDLE);
                m_devFuncs->vkQueueWaitIdle(m_window->graphicsQueue());

                m_devFuncs->vkFreeCommandBuffers(m_window->device(), m_window->graphicsCommandPool(), 1, &commandBuffer);
            }

            void createUniformBuffers() {
                VkDeviceSize bufferSize = sizeof(UniformBufferObject);
                uniformBuffers.resize(m_window->MAX_CONCURRENT_FRAME_COUNT);
                uniformBuffersMemory.resize(m_window->MAX_CONCURRENT_FRAME_COUNT);

                for (size_t i = 0; i < m_window->MAX_CONCURRENT_FRAME_COUNT; i++) {
                    createBuffer(bufferSize, VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, uniformBuffers[i], uniformBuffersMemory[i]);
                }
            }

            void createBuffer(VkDeviceSize size, VkBufferUsageFlags usage, VkMemoryPropertyFlags properties, VkBuffer& buffer, VkDeviceMemory& bufferMemory) {
                VkBufferCreateInfo bufferInfo{};
                bufferInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
                bufferInfo.size = size;
                bufferInfo.usage = usage;
                bufferInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

                if (m_devFuncs->vkCreateBuffer(m_window->device(), &bufferInfo, nullptr, &buffer) != VK_SUCCESS)
                    throw std::runtime_error("failed to create buffer!");

                VkMemoryRequirements memReq;
                m_devFuncs->vkGetBufferMemoryRequirements(m_window->device(), buffer, &memReq);

                VkMemoryAllocateInfo allocInfo{};
                allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
                allocInfo.allocationSize = memReq.size;
                allocInfo.memoryTypeIndex = m_window->hostVisibleMemoryIndex();

                if (m_devFuncs->vkAllocateMemory(m_window->device(), &allocInfo, nullptr, &bufferMemory) != VK_SUCCESS)
                    throw std::runtime_error("failed to allocate memory!");

                m_devFuncs->vkBindBufferMemory(m_window->device(), buffer, bufferMemory, 0);
            }

            VkShaderModule createShader(const QString& path) {
                QFile file(path);
                if (!file.open(QIODevice::ReadOnly)) {
                    qWarning("failed to read shader at %s", qPrintable(path)); return VK_NULL_HANDLE;
                }

                QByteArray blob = file.readAll();
                file.close();

                VkShaderModuleCreateInfo shaderCreateInfo{};
                shaderCreateInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
                shaderCreateInfo.pCode = reinterpret_cast<const uint32_t*>(blob.constData());
                shaderCreateInfo.codeSize = blob.size();

                VkShaderModule shaderModule;
                VkResult err;

                err = m_devFuncs->vkCreateShaderModule(m_window->device(), &shaderCreateInfo, nullptr, &shaderModule);
                if (err != VK_SUCCESS) {
                    qFatal("Failed to create shader from %s: %d", path, err); return VK_NULL_HANDLE;
                }

                return shaderModule;
            }

            void startNextFrame() override {
                m_rotation += 0.1f;
                VkDevice dev = m_window->device();
                const QSize sz = m_window->swapChainImageSize();

                VkClearColorValue clearColor = { {.5f, .5f, .5f, 1.0f} };
                VkClearDepthStencilValue depthStencilClearColor = { 1.0f, 0 };
                VkClearValue clearValues[3];
                memset(clearValues, 0, sizeof(clearValues));
                clearValues[0].color = clearValues[2].color = clearColor;
                clearValues[1].depthStencil = depthStencilClearColor;

                VkRenderPassBeginInfo renderPassBeginInfo{};
                renderPassBeginInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
                renderPassBeginInfo.renderPass = m_window->defaultRenderPass();
                renderPassBeginInfo.framebuffer = m_window->currentFramebuffer();
                renderPassBeginInfo.renderArea.extent.width = sz.width();
                renderPassBeginInfo.renderArea.extent.height = sz.height();
                renderPassBeginInfo.clearValueCount = m_window->sampleCountFlagBits() > VK_SAMPLE_COUNT_1_BIT ? 3 : 2;
                renderPassBeginInfo.pClearValues = clearValues;
                VkCommandBuffer cmdBuf = m_window->currentCommandBuffer();
                m_devFuncs->vkCmdBeginRenderPass(cmdBuf, &renderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

                quint8* p;
                VkResult err = m_devFuncs->vkMapMemory(dev, m_buffVertexMem, m_uniformBufferInfo[m_window->currentFrame()].offset, UNIFORM_DATA_SIZE, 0, reinterpret_cast<void**>(&p));
                if (err != VK_SUCCESS)
                    qFatal("Failed to map memory: %d", err);
                QMatrix4x4 m = m_proj;
                m.rotate(m_rotation, 0, 1, 0);
                memcpy(p, m.constData(), 16 * sizeof(float));
                m_devFuncs->vkUnmapMemory(dev, m_buffVertexMem);

                m_devFuncs->vkCmdBindPipeline(cmdBuf, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
                m_devFuncs->vkCmdBindDescriptorSets(cmdBuf, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipelineLayout, 0, 1, &m_descSet[m_window->currentFrame()], 0, nullptr);

                VkDeviceSize vbOffset = 0;
                m_devFuncs->vkCmdBindVertexBuffers(cmdBuf, 0, 1, &m_buffVertex, &vbOffset);
                m_devFuncs->vkCmdBindIndexBuffer(cmdBuf, m_buffIndex, 0, VK_INDEX_TYPE_UINT32);

                VkViewport vp;
                vp.x = vp.y = 0;
                vp.width = sz.width();
                vp.height = sz.height();
                vp.minDepth = 0;
                vp.maxDepth = 1;
                m_devFuncs->vkCmdSetViewport(cmdBuf, 0, 1, &vp);

                VkRect2D sc;
                sc.offset.x = sc.offset.y = 0;
                sc.extent.width = vp.width;
                sc.extent.height = vp.height;
                m_devFuncs->vkCmdSetScissor(cmdBuf, 0, 1, &sc);

                if (renderWire) {
                    m_devFuncs->vkCmdBindIndexBuffer(cmdBuf, m_buffIndexGraph, 0, VK_INDEX_TYPE_UINT32);
                    m_devFuncs->vkCmdBindPipeline(cmdBuf, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline_wire);
                    m_devFuncs->vkCmdSetLineWidth(cmdBuf, 1.0f);
                    //m_devFuncs->vkCmdDraw(cb, 342, 1, 0, 0);
                    m_devFuncs->vkCmdDrawIndexed(cmdBuf, static_cast<uint32_t>(indicesGraph.size()), 1, 0, 0, 0);
                }
                if (renderFace) {
                    m_devFuncs->vkCmdBindIndexBuffer(cmdBuf, m_buffIndex, 0, VK_INDEX_TYPE_UINT32);
                    m_devFuncs->vkCmdBindPipeline(cmdBuf, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline);
                    m_devFuncs->vkCmdSetLineWidth(cmdBuf, 3.0f);
                    //m_devFuncs->vkCmdDraw(cb, 342, 1, 0, 0);
                    m_devFuncs->vkCmdDrawIndexed(cmdBuf, static_cast<uint32_t>(indices.size()), 1, 0, 0, 0);
                    //m_devFuncs->vkCmdDraw(cmdBuf, vertices.size(), 1, 0, 0);
                }
                if (renderPoints) {
                    m_devFuncs->vkCmdBindIndexBuffer(cmdBuf, m_buffIndex, 0, VK_INDEX_TYPE_UINT32);
                    m_devFuncs->vkCmdBindPipeline(cmdBuf, VK_PIPELINE_BIND_POINT_GRAPHICS, m_pipeline_point);
                    m_devFuncs->vkCmdSetLineWidth(cmdBuf, 10.0f);
                    //m_devFuncs->vkCmdDraw(cb, 342, 1, 0, 0);
                    m_devFuncs->vkCmdDrawIndexed(cmdBuf, static_cast<uint32_t>(indices.size()), 1, 0, 0, 0);
                }


                m_devFuncs->vkCmdEndRenderPass(cmdBuf);

                m_window->frameReady();
                m_window->requestUpdate();
            }

            void toggleWireframe() {
                renderWire = !renderWire;
            }
            void toggleFaceRender() {
                renderFace = !renderFace;
            }
            void togglePointRender() {
                renderPoints = !renderPoints;
            }


        protected:
            QVulkanWindow* m_window;
            QVulkanDeviceFunctions* m_devFuncs;

            VkBuffer m_buffVertex;
            VkBuffer m_buffIndex;
            VkBuffer m_buffIndexGraph;
            VkDeviceMemory m_buffIndexMem;
            VkDeviceMemory m_buffIndexGraphMem;
            VkDeviceMemory m_buffVertexMem;
            VkDescriptorBufferInfo m_uniformBufferInfo[QVulkanWindow::MAX_CONCURRENT_FRAME_COUNT];

            VkDescriptorPool m_descPool;
            VkDescriptorSetLayout m_descSetLayout;
            VkDescriptorSet m_descSet[QVulkanWindow::MAX_CONCURRENT_FRAME_COUNT];
            VkPipelineCache m_pipelineCache;
            VkPipelineLayout m_pipelineLayout;
            VkPipeline m_pipeline;
            VkPipeline m_pipeline_wire;
            VkPipeline m_pipeline_point;

            std::vector<VkBuffer> uniformBuffers;
            std::vector<VkDeviceMemory> uniformBuffersMemory;

            bool renderWire = true;
            bool renderPoints = true;
            bool renderFace = true;
            QMatrix4x4 m_proj;
            float m_rotation = 0.4f;
            CGAL_POLY_TYPE m_poly;
            std::vector<Vertex> vertices;
            std::vector<uint32_t> indices;
            std::vector<uint32_t> indicesGraph;

    };
}

#endif // !POLYHEDRON_RENDERER
