// Copyright (c) 2023 Texas A&M University (United States).
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

#include <QVulkanWindow>
#include <QVulkanFunctions>
#include <QFile>

static float vertexData[] = { // Y up, front = CCW
     0.0f,   0.5f,   1.0f, 0.0f, 0.0f,
    -0.5f,  -0.5f,   0.0f, 1.0f, 0.0f,
     0.5f,  -0.5f,   0.0f, 0.0f, 1.0f
};

static const int UNIFORM_DATA_SIZE = 16 * sizeof(float);

static inline VkDeviceSize aligned(VkDeviceSize v, VkDeviceSize byteAlign)
{
    return (v + byteAlign - 1) & ~(byteAlign - 1);
}


class Polyhedron_renderer : public QVulkanWindowRenderer {
public:
	Polyhedron_renderer(QVulkanWindow* w) : m_window(w){}
	void initResources() override {
		qDebug("initResources call");

		VkDevice dev = m_window->device();
        m_devFuncs = m_window->vulkanInstance()->deviceFunctions(dev);

        const int concFrameCount = m_window->concurrentFrameCount();
        const VkPhysicalDeviceLimits* pDevLimits = &m_window->physicalDeviceProperties()->limits;
        const VkDeviceSize uniAlign = pDevLimits->minUniformBufferOffsetAlignment;
        qDebug("uniform buffer offset alignment is %u", (uint)uniAlign);

        VkBufferCreateInfo buffInfo{};
        buffInfo.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;

        const VkDeviceSize vertexAllocSize = aligned(sizeof(vertexData), uniAlign);

	}
	void initSwapChainResources() override;
	void releaseSwapChainResources() override;
	void releaseResources() override;

	void startNextFrame() override {
        m_green += 0.005f;
        if (m_green > 1.0f)
            m_green = 0.0f;
        VkClearColorValue clearColor = { {.0f, m_green, .0f, 1.0f} };
        VkClearDepthStencilValue depthStencilClearColor = { 1.0f, 0 };
        VkClearValue clearValues[2];
        memset(clearValues, 0, sizeof(clearValues));
        clearValues[0].color = clearColor;
        clearValues[1].depthStencil = depthStencilClearColor;

        VkRenderPassBeginInfo renderPassBeginInfo{};
        renderPassBeginInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
        renderPassBeginInfo.renderPass = m_window->defaultRenderPass();
        renderPassBeginInfo.framebuffer = m_window->currentFramebuffer();
        const QSize size = m_window->swapChainImageSize();
        renderPassBeginInfo.renderArea.extent.width = size.width();
        renderPassBeginInfo.renderArea.extent.height = size.width();
        renderPassBeginInfo.clearValueCount = 2;
        renderPassBeginInfo.pClearValues = clearValues;
        VkCommandBuffer cmdBuf = m_window->currentCommandBuffer();
        m_devFuncs->vkCmdBeginRenderPass(cmdBuf, &renderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);
        m_devFuncs->vkCmdEndRenderPass(cmdBuf);

        m_window->frameReady();
        m_window->requestUpdate();
    }

protected:
	QVulkanWindow* m_window;
	QVulkanDeviceFunctions* m_devFuncs;
    float m_green;

};
#endif // !POLYHEDRON_RENDERER

