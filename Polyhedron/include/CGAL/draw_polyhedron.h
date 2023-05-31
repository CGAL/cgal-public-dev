// Copyright (c) 2018-2020  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POLYHEDRON_H
#define CGAL_DRAW_POLYHEDRON_H

#include <CGAL/license/Polyhedron.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Random.h>

#include <QVulkanWindow>
#include <QVulkanInstance>
#include <QVulkanFunctions>
namespace CGAL
{

// Specialization of draw function.
#define CGAL_POLY_TYPE CGAL::Polyhedron_3 \
  <PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>

class VulkanRenderer : public QVulkanWindowRenderer {
public: 
    VulkanRenderer(QVulkanWindow* w) : m_window(w), m_devFuncs(nullptr), m_green(0.0f) {};

    void initResources() override {
        qDebug("initResources");
        m_devFuncs = m_window->vulkanInstance()->deviceFunctions(m_window->device());
    }
    void initSwapChainResources() override {}
    void releaseSwapChainResources() override {}
    void releaseResources() override {}

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

private:
    QVulkanWindow* m_window;
    QVulkanDeviceFunctions* m_devFuncs;
    float m_green = 0;
};

class VulkanWindow : public QVulkanWindow {
public:
    QVulkanWindowRenderer* createRenderer() override {
        return new VulkanRenderer(this);
    }
};

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc>
void draw(const CGAL_POLY_TYPE& apoly,
          const char* title="Polyhedron Basic Viewer",
          bool nofill=false)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    //CGAL::Qt::init_ogl_context(4,3);
    int argc=1;
    const char* argv[2]={"polyhedron_viewer", nullptr};
    QApplication app(argc,const_cast<char**>(argv));
    QVulkanInstance inst;
    inst.setLayers(QByteArrayList() << "VK_LAYER_LUNARG_standard_validation");
    if (!inst.create()) {
        qFatal("Failed to create a Vulkan instance: %d", inst.errorCode());
    }

    VulkanWindow w;
    w.setVulkanInstance(&inst);
    w.showMaximized();
    //VulkanRenderer* renderer =  w.createRenderer();
    //renderer->startNextFrame();
    //w.requestUpdate();
    //SimpleFaceGraphViewerQt  mainwindow(app.activeWindow(), apoly, title, nofill);
    //mainwindow.show();
    app.exec();
  }
}

#undef CGAL_POLY_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYHEDRON_H
