#ifndef CGAL_BASIC_VIEWER_GLFW_IMPL_H
#define CGAL_BASIC_VIEWER_GLFW_IMPL_H

#include "Basic_viewer.h"

namespace CGAL 
{
namespace GLFW 
{

  inline 
  void Basic_viewer::error_callback(int error, const char *description)
  {
    std::cerr << "GLFW returned an error:\n\t" << description << "(" << error << ")\n";
  }

  inline
  GLFWwindow* Basic_viewer::create_window(int width, int height, const char* title, bool hidden)
  {
    // Initialise GLFW
    if (!glfwInit())
    {
      std::cerr << "Could not start GLFW\n";
      exit(EXIT_FAILURE);
    }

    // OpenGL 2.1 with compatibilty
    // if (hidden)
    // {
    //   glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    // }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Enable the GLFW runtime error callback function defined previously.
    glfwSetErrorCallback(error_callback);

    // Set additional window options
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
    glfwWindowHint(GLFW_SAMPLES, WINDOW_SAMPLES); // MSAA

    // Create window using GLFW
    GLFWwindow* window = glfwCreateWindow(width, height, title, nullptr, nullptr);

    // Ensure the window is set up correctly
    if (!window)
    {
      std::cerr << "Could not open GLFW window\n";
      glfwTerminate();
      exit(EXIT_FAILURE);
    }

    // Let the window be the current OpenGL context and initialise glad
    glfwMakeContextCurrent(window);

    // Initialized GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
      std::cerr << "Failed to initialized GLAD!";
      glfwDestroyWindow(window);
      glfwTerminate();
      exit(EXIT_FAILURE);
    }

    // Print various OpenGL information to stdout
    std::cout << glGetString(GL_VENDOR) << ": " << glGetString(GL_RENDERER) << '\n';
    std::cout << "GLFW\t " << glfwGetVersionString() << '\n';
    std::cout << "OpenGL\t " << glGetString(GL_VERSION) << '\n';
    std::cout << "GLSL\t " << glGetString(GL_SHADING_LANGUAGE_VERSION) << "\n\n";

    return window;
  }

  inline
  Basic_viewer::Basic_viewer(
    const Graphics_scene* graphicScene,
    const char* title,
    bool drawVertices,
    bool drawEdges,
    bool drawFaces,
    bool drawRays,
    bool drawLines, 
    bool useMonoColor,
    bool inverseNormal,
    bool flatShading
  ) : m_scene(graphicScene),
    m_title(title),
    m_drawVertices(drawVertices),
    m_drawEdges(drawEdges),
    m_drawFaces(drawFaces),
    m_drawRays(drawRays),
    m_drawLines(drawLines),
    m_useMonoColor(useMonoColor),
    m_inverseNormal(inverseNormal),
    m_flatShading(flatShading)
  {
    initialize();
    GLint maxGeometryOutputVertices = 0;
    glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES, &maxGeometryOutputVertices);
    std::cout << "Maximum geometry output vertices: " << maxGeometryOutputVertices << "\n";
    GLint maxGeometryOutputComponents = 0;
    glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_COMPONENTS, &maxGeometryOutputComponents);
    std::cout << "Maximum geometry output Component: " << maxGeometryOutputComponents << "\n";
  }

  inline
  void Basic_viewer::initialize(bool screenshotOnly)
  {
    m_window = create_window(m_windowSize.x(), m_windowSize.y(), m_title, screenshotOnly);

    m_aspectRatio = static_cast<float>(m_windowSize.x())/m_windowSize.y();

    if (!screenshotOnly)
    {
      // Set event callbacks 
      glfwSetWindowUserPointer(m_window, this);
      glfwSetKeyCallback(m_window, key_callback);
      glfwSetCursorPosCallback(m_window, cursor_callback);
      glfwSetMouseButtonCallback(m_window, mouse_btn_callback);
      glfwSetScrollCallback(m_window, scroll_callback);
      glfwSetFramebufferSizeCallback(m_window, window_size_callback);
    }

    GLint openglMajorVersion, openglMinorVersion;
    glGetIntegerv(GL_MAJOR_VERSION, &openglMajorVersion);
    glGetIntegerv(GL_MINOR_VERSION, &openglMinorVersion);

    if (openglMajorVersion > 4 || openglMajorVersion == 4 && openglMinorVersion >= 3)
    {
      m_isOpengl4_3 = true;
    }

    compile_shaders();
    initialize_camera();
    initialize_buffers();
    initialize_and_load_world_axis();
    initialize_and_load_clipping_plane();
    if (!screenshotOnly)
    {
      initialize_keys_actions();
    }
  }

  inline
  void Basic_viewer::show()
  {
    float elapsedTime = 0.0f;
    float lastFrame = 0.0;
    while (!glfwWindowShouldClose(m_window))
    {
      float currentFrame = static_cast<float>(glfwGetTime());
      m_deltaTime = currentFrame - lastFrame;
      lastFrame = currentFrame;
      if (m_deltaTime < 1e-3) m_deltaTime = 1e-3;

      handle_events(m_deltaTime);
      if (need_update()) 
      {
        render_scene(m_deltaTime);
      }
      print_application_state(elapsedTime, m_deltaTime);
    }

    clear_application();
  }

  void Basic_viewer::clear_application()
  {
    m_shaderPl.destroy();
    m_shaderFace.destroy();
    m_shaderLine.destroy();
    m_shaderArrow.destroy();
    m_shaderPlane.destroy();
    glDeleteBuffers(NB_GL_BUFFERS, m_vbo);
    glDeleteVertexArrays(NB_VAO_BUFFERS, m_vao);
    glfwDestroyWindow(m_window);
    glfwTerminate();
  }

  inline
  void Basic_viewer::make_screenshot(const std::string& filePath)
  {
    draw(0);
    glfwSwapBuffers(m_window);
    capture_screenshot(filePath);
    clear_application();
  }

  void generate_grid(Line_renderer& renderer, const vec3f& color, float size, int nbSubdivisions=10)
  {
    for (unsigned int i = 0; i <= nbSubdivisions; ++i)
    {
      float pos = float(size * (2.0 * i / nbSubdivisions - 1.0));
      renderer.add_line(vec3f(pos, -size, 0.f), vec3f(pos, size, 0.f), color);
      renderer.add_line(vec3f(-size, pos, 0.f), vec3f(size, pos, 0.f), color);
    }
  }

  inline
  void Basic_viewer::compile_shaders()
  {
    const char* FACE_VERTEX = m_isOpengl4_3 ? VERTEX_SOURCE_COLOR : VERTEX_SOURCE_COLOR_COMP;
    const char* FACE_FRAGMENT = m_isOpengl4_3 ? FRAGMENT_SOURCE_COLOR : FRAGMENT_SOURCE_COLOR_COMP;
    const char* PL_VERTEX = m_isOpengl4_3 ? VERTEX_SOURCE_P_L : VERTEX_SOURCE_P_L_COMP;
    const char* PL_FRAGMENT = m_isOpengl4_3 ? FRAGMENT_SOURCE_P_L : FRAGMENT_SOURCE_P_L_COMP;
    const char* PLANE_VERTEX = VERTEX_SOURCE_CLIPPING_PLANE;
    const char* PLANE_FRAGMENT = FRAGMENT_SOURCE_CLIPPING_PLANE;

    m_shaderPl = Shader::create_shader(PL_VERTEX, PL_FRAGMENT);
    m_shaderFace = Shader::create_shader(FACE_VERTEX, FACE_FRAGMENT);
    m_shaderPlane = Shader::create_shader(PLANE_VERTEX, PLANE_FRAGMENT);

    const char* POINT_GEOMETRY = GEOMETRY_SOURCE_SPHERE;
    m_shaderSphere = Shader::create_shader(PL_VERTEX, PL_FRAGMENT, POINT_GEOMETRY);

    const char* EDGE_GEOMETRY = GEOMETRY_SOURCE_CYLINDER;
    m_shaderCylinder = Shader::create_shader(PL_VERTEX, PL_FRAGMENT, EDGE_GEOMETRY);

    // For world axis and grid 
    const char* LINE_VERTEX = VERTEX_SOURCE_LINE;
    const char* LINE_GEOMETRY = GEOMETRY_SOURCE_LINE;
    const char* LINE_FRAGMENT = FRAGMENT_SOURCE_LINE;
    m_shaderLine = Shader::create_shader(LINE_VERTEX, LINE_FRAGMENT, LINE_GEOMETRY);

    const char* ARROW_GEOMETRY = GEOMETRY_SOURCE_ARROW;
    m_shaderArrow = Shader::create_shader(LINE_VERTEX, LINE_FRAGMENT, ARROW_GEOMETRY);

    const char* NORMAL_VERTEX = VERTEX_SOURCE_NORMAL;
    const char* NORMAL_GEOMETRY = GEOMETRY_SOURCE_NORMAL;
    const char* NORMAL_FRAGMENT = FRAGMENT_SOURCE_NORMAL;
    m_shaderNormal = Shader::create_shader(NORMAL_VERTEX, NORMAL_FRAGMENT, NORMAL_GEOMETRY);
  }

  inline 
  void Basic_viewer::initialize_camera()
  {
    vec3f pmin(
      m_scene->bounding_box().xmin(),
      m_scene->bounding_box().ymin(),
      m_scene->bounding_box().zmin());

    vec3f pmax(
      m_scene->bounding_box().xmax(),
      m_scene->bounding_box().ymax(),
      m_scene->bounding_box().zmax());

    m_boundingBox = { pmin, pmax };

    m_camera.lookat(pmin, pmax);  
  }

  /// @brief 
  inline 
  void Basic_viewer::initialize_and_load_world_axis()
  {
    // World axis initialization
    m_worldAxisRenderer.initialize_buffers();
    m_worldAxisRenderer.set_width(3.f);
    m_worldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitX(), vec3f(1, 0, 0)); // x-axis
    m_worldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitY(), vec3f(0, 1, 0)); // y-axis
    m_worldAxisRenderer.add_line(vec3f::Zero(), .1f*vec3f::UnitZ(), vec3f(0, 0, 1)); // z-axis
    m_worldAxisRenderer.load_buffers();

    float cameraSize = m_camera.get_size() * 0.5;
    // XY grid axis initialization
    m_XYAxisRenderer.initialize_buffers();
    m_XYAxisRenderer.set_width(5.f);
    m_XYAxisRenderer.add_line(vec3f::Zero(), cameraSize*vec3f::UnitX(), vec3f(1, 0, 0)); // x-axis
    m_XYAxisRenderer.add_line(vec3f::Zero(), cameraSize*vec3f::UnitY(), vec3f(0, 1, 0)); // y-axis
    m_XYAxisRenderer.load_buffers();

    // XY grid initialization 
    m_XYGridRenderer.initialize_buffers();
    m_XYGridRenderer.set_width(2.f);
    m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitX(), vec3f(.8f, .8f, .8f)); // -x-axis
    m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitY(), vec3f(.8f, .8f, .8f)); // -y-axis
    m_XYGridRenderer.add_line(vec3f::Zero(), -2.f*cameraSize*vec3f::UnitZ(), vec3f(.8f, .8f, .8f)); // -z-axis

    vec3f color(.8f, .8f, .8f);

    generate_grid(m_XYGridRenderer, color, cameraSize);
    
    m_XYGridRenderer.load_buffers();
  }

  inline
  void Basic_viewer::load_buffer(int i, int location, const std::vector<float>& vector, int dataCount)
  {
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo[i]);

    glBufferData(GL_ARRAY_BUFFER, vector.size() * sizeof(float), vector.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(location, dataCount, GL_FLOAT, GL_FALSE, dataCount * sizeof(float), nullptr);

    glEnableVertexAttribArray(location);
  }

  inline
  void Basic_viewer::load_buffer(int i, int location, int gsEnum, int dataCount)
  {
    const auto& vector = m_scene->get_array_of_index(gsEnum);
    load_buffer(i, location, vector, dataCount);
  }

  inline
  void Basic_viewer::initialize_buffers()
  {
    glGenBuffers(NB_GL_BUFFERS, m_vbo);
    glGenVertexArrays(NB_VAO_BUFFERS, m_vao);
  }

  inline
  void Basic_viewer::load_scene()
  {
    unsigned int bufn = 0;

    // 1) POINT SHADER

    // 1.1) Mono points
    glBindVertexArray(m_vao[VAO_MONO_POINTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_POINTS, 3);

    // 1.2) Color points
    glBindVertexArray(m_vao[VAO_COLORED_POINTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_POINTS, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_POINTS, 3);

    // 2) SEGMENT SHADER

    // 2.1) Mono segments
    glBindVertexArray(m_vao[VAO_MONO_SEGMENTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_SEGMENTS, 3);

    // 2.2) Colored segments
    glBindVertexArray(m_vao[VAO_COLORED_SEGMENTS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_SEGMENTS, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_SEGMENTS, 3);

    // 3) RAYS SHADER

    // 2.1) Mono segments
    glBindVertexArray(m_vao[VAO_MONO_RAYS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_RAYS, 3);

    // 2.2) Colored segments
    glBindVertexArray(m_vao[VAO_COLORED_RAYS]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_RAYS, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_RAYS, 3);

    // 4) LINES SHADER

    // 2.1) Mono lines
    glBindVertexArray(m_vao[VAO_MONO_LINES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_LINES, 3);

    // 2.2) Colored lines
    glBindVertexArray(m_vao[VAO_COLORED_LINES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_LINES, 3);
    load_buffer(bufn++, 1, Graphics_scene::COLOR_LINES, 3);

    // 5) FACE SHADER

    // 5.1) Mono faces
    glBindVertexArray(m_vao[VAO_MONO_FACES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_MONO_FACES, 3);
    if (m_flatShading)
    {
      load_buffer(bufn++, 1, Graphics_scene::FLAT_NORMAL_MONO_FACES, 3);
    }
    else
    {
      load_buffer(bufn++, 1, Graphics_scene::SMOOTH_NORMAL_MONO_FACES, 3);
    }

    // 5.2) Colored faces
    glBindVertexArray(m_vao[VAO_COLORED_FACES]);
    load_buffer(bufn++, 0, Graphics_scene::POS_COLORED_FACES, 3);
    if (m_flatShading)
    {
      load_buffer(bufn++, 1, Graphics_scene::FLAT_NORMAL_COLORED_FACES, 3);
    }
    else
    {
      load_buffer(bufn++, 1, Graphics_scene::SMOOTH_NORMAL_COLORED_FACES, 3);
    }
    load_buffer(bufn++, 2, Graphics_scene::COLOR_FACES, 3);

    m_areBuffersInitialized = true;
  }

  inline 
  CGAL::Plane_3<Basic_viewer::Local_kernel> Basic_viewer::clipping_plane() const
  {
    mat4f CPM = m_clippingPlane.get_matrix();
    CGAL::Aff_transformation_3<Basic_viewer::Local_kernel> aff(
      CPM(0, 0), CPM(0, 1), CPM(0, 2), CPM(0, 3),
      CPM(1, 0), CPM(1, 1), CPM(1, 2), CPM(1, 3),
      CPM(2, 0), CPM(2, 1), CPM(2, 2), CPM(2, 3)
    );

    CGAL::Plane_3<Local_kernel> p3(0, 0, 1, 0);
    return p3.transform(aff);
  }

  inline 
  void Basic_viewer::compute_model_view_projection_matrix(const float deltaTime)
  {
    m_camera.update(deltaTime);
    m_clippingPlane.update(deltaTime);

    if (m_animationController.is_running())
    {
      AnimationKeyFrame animationFrame = m_animationController.run();
      m_camera.set_orientation(animationFrame.orientation);
      m_camera.set_position(animationFrame.position);
    }

    m_viewMatrix = m_camera.view();
    m_projectionMatrix = m_camera.projection(m_windowSize.x(), m_windowSize.y());

    m_viewProjectionMatrix = m_projectionMatrix * m_viewMatrix;
  }

  inline
  void Basic_viewer::update_uniforms(const float deltaTime)
  {
    compute_model_view_projection_matrix(deltaTime);

    // ================================================================

    update_face_uniforms();
    update_pl_uniforms();
    update_clipping_uniforms();
  }

  inline
  void Basic_viewer::update_face_uniforms()
  {
    m_shaderFace.use();

    m_shaderFace.set_mat4f("u_Mvp", m_viewProjectionMatrix.data());
    m_shaderFace.set_mat4f("u_Mv",  m_viewMatrix.data());

    m_shaderFace.set_vec4f("u_LightPos",  m_lightPosition.data());
    m_shaderFace.set_vec4f("u_LightDiff", m_diffuseColor.data());
    m_shaderFace.set_vec4f("u_LightSpec", m_specularColor.data());
    m_shaderFace.set_vec4f("u_LightAmb",  m_ambientColor.data());
    m_shaderFace.set_float("u_SpecPower", m_shininess);

    m_shaderFace.set_vec4f("u_ClipPlane",  m_clipPlane.data());
    m_shaderFace.set_vec4f("u_PointPlane", m_pointPlane.data());
    m_shaderFace.set_float("u_RenderingTransparency", m_clippingPlane.get_transparency());
  }

  inline 
  void Basic_viewer::update_point_uniforms()
  {
    m_shaderSphere.use();

    bool half = m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    m_shaderSphere.set_float("u_RenderingMode", static_cast<float>(mode));

    m_shaderSphere.set_mat4f("u_Mvp", m_viewProjectionMatrix.data());
    m_shaderSphere.set_vec4f("u_ClipPlane",  m_clipPlane.data());
    m_shaderSphere.set_vec4f("u_PointPlane", m_pointPlane.data());
    m_shaderSphere.set_float("u_UseGeometryShader", 1.0);
    m_shaderSphere.set_float("u_PointSize", m_sizeVertices);
    m_shaderSphere.set_float("u_Radius", m_camera.get_radius()*m_sizeVertices*0.001);
  }

  inline 
  void Basic_viewer::update_edge_uniforms()
  {
    m_shaderCylinder.use();

    bool half = m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    m_shaderCylinder.set_float("u_RenderingMode", static_cast<float>(mode));

    m_shaderCylinder.set_mat4f("u_Mvp", m_viewProjectionMatrix.data());
    m_shaderCylinder.set_vec4f("u_ClipPlane",  m_clipPlane.data());
    m_shaderCylinder.set_vec4f("u_PointPlane", m_pointPlane.data());
    m_shaderCylinder.set_float("u_PointSize", m_sizeEdges);
    m_shaderCylinder.set_float("u_UseGeometryShader", 1.0);
    m_shaderCylinder.set_float("u_Radius", m_camera.get_radius()*m_sizeEdges*0.001);
  }

  inline
  void Basic_viewer::update_pl_uniforms()
  {
    m_shaderPl.use();

    bool half = m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY;
    auto mode = half ? RenderingMode::DRAW_INSIDE_ONLY : RenderingMode::DRAW_ALL;

    m_shaderPl.set_float("u_RenderingMode", static_cast<float>(mode));

    m_shaderPl.set_mat4f("u_Mvp", m_viewProjectionMatrix.data());
    m_shaderPl.set_vec4f("u_ClipPlane",  m_clipPlane.data());
    m_shaderPl.set_vec4f("u_PointPlane", m_pointPlane.data());
    m_shaderPl.set_float("u_PointSize", m_sizeVertices);
    m_shaderPl.set_float("u_UseGeometryShader", 0.0);
  }

  inline
  void Basic_viewer::update_clipping_uniforms()
  {
    mat4f clippingModelMatrix = m_clippingPlane.get_matrix();

    m_pointPlane = clippingModelMatrix * vec4f(0, 0, 0, 1);
    m_clipPlane = clippingModelMatrix * vec4f(0, 0, 1, 0);

    m_shaderPlane.use();
    m_shaderPlane.set_mat4f("u_Vp", m_viewProjectionMatrix.data());
    m_shaderPlane.set_mat4f("u_M",  clippingModelMatrix.data());
  }

  inline
  void Basic_viewer::update_world_axis_uniforms()
  {
    mat4f view = m_viewMatrix;

    // we only want the rotation part of the view matrix  
    mat3f rotation = view.block<3,3>(0,0);

    mat4f rotation4x4 = mat4f::Identity();
    rotation4x4.block<3,3>(0,0) = rotation;

    float halfWidth = m_aspectRatio * 0.1f;
    float halfHeight = 0.1f;
    mat4f projection = utils::ortho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0f, 1.0f);

    mat4f translation = transform::translation(vec3f(halfWidth - 0.1f*m_aspectRatio, halfHeight - 0.1f, 0.0f));

    mat4f mvp = projection * rotation4x4 * translation;

    m_shaderArrow.use();
    m_shaderArrow.set_mat4f("u_Mvp", mvp.data()); 
    m_shaderArrow.set_float("u_SceneRadius", 1.0f); 
  }

  inline
  void Basic_viewer::update_XY_axis_uniforms()
  {
    m_shaderArrow.use();
    m_shaderArrow.set_mat4f("u_Mvp", m_viewProjectionMatrix.data());
    m_shaderArrow.set_float("u_SceneRadius", m_camera.get_radius()); 
  }

  inline
  void Basic_viewer::update_XY_grid_uniforms()
  {
    m_shaderLine.use();
    m_shaderLine.set_mat4f("u_Mvp", m_viewProjectionMatrix.data());
  }

  inline 
  void Basic_viewer::update_normals_uniforms()
  {
    m_shaderNormal.use();
   
    vec4f color = color_to_normalized_vec4(m_normalsMonoColor);
    m_shaderNormal.set_mat4f("u_Mv", m_viewMatrix.data());
    m_shaderNormal.set_vec4f("u_Color", color.data());
    if (m_useNormalMonoColor)
    {
      m_shaderNormal.set_float("u_UseMonoColor", 1.0);
    }
    else
    {
      m_shaderNormal.set_float("u_UseMonoColor", 0.0);
    }
    m_shaderNormal.set_mat4f("u_Projection", m_projectionMatrix.data());
    m_shaderNormal.set_float("u_Factor", m_normalHeightFactor);
    m_shaderNormal.set_float("u_SceneRadius", m_camera.get_radius());

  }

  inline 
  bool Basic_viewer::need_update() const
  {
    return m_camera.need_update() 
        || m_clippingPlane.need_update() 
        || m_animationController.is_running()
        || has_active_actions() 
        ; 
  }

  inline 
  void Basic_viewer::render_scene(const float deltaTime)
  {
    draw(deltaTime);
    glfwSwapBuffers(m_window);
  }

  inline
  void Basic_viewer::draw(const float deltaTime)
  {
    if (!m_areBuffersInitialized)
    {
      load_scene();
    }

    glClearColor(1.0f, 1.0f, 1.0f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_LINE_SMOOTH);

    // update_uniforms(deltaTime);

    compute_model_view_projection_matrix(deltaTime);

    update_pl_uniforms();
    if (m_drawRays)
    {
      draw_rays();
    }
    if (m_drawLines)
    {
      draw_lines();
    }
    if (m_drawEdges)
    {
      update_pl_uniforms();
      if (m_drawCylinderEdge)
      {
        update_edge_uniforms();
      }
      draw_edges();
    }
    if (m_drawVertices)
    {
      update_pl_uniforms();
      if (m_drawSphereVertex)
      {
        update_point_uniforms();
      }
      draw_vertices();
    }

    if (clipping_plane_enable())
    {
      update_clipping_uniforms();
      render_clipping_plane();
    }

    if (m_drawNormals)
    {
      update_normals_uniforms();
      draw_normals();
    }  
        
    if (m_drawFaces)
    {
      update_face_uniforms();
      draw_faces();
    }
    
    if (m_drawWorldAxis)  
    {
      update_world_axis_uniforms();
      draw_world_axis();
    }

    if (m_drawXYGrid) 
    { 
      draw_xy_grid();
    }
  }

  inline
  vec4f Basic_viewer::color_to_normalized_vec4(const CGAL::IO::Color& c) const
  {
    return { static_cast<float>(c.red()) / 255, static_cast<float>(c.green()) / 255, static_cast<float>(c.blue()) / 255, 1.0f };
  }

  inline
  void Basic_viewer::draw_faces()
  {
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(2.0, 1.0); 
    glDepthFunc(GL_LESS);
    if (m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF)
    {
      // The z-buffer will prevent transparent objects from being displayed behind other transparent objects.
      // Before rendering all transparent objects, disable z-testing first.

      // 1. draw solid first
      draw_faces_bis(RenderingMode::DRAW_INSIDE_ONLY);

      // 2. draw transparent layer second with back face culling to avoid messy triangles
      glDepthMask(false); // disable z-testing
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_CULL_FACE);
      glCullFace(GL_BACK);
      glFrontFace(GL_CW);
      draw_faces_bis(RenderingMode::DRAW_OUTSIDE_ONLY);

      // 3. draw solid again without culling and blend to make sure the solid mesh is visible
      glDepthMask(true); // enable z-testing
      glDisable(GL_CULL_FACE);
      glDisable(GL_BLEND);
      draw_faces_bis(RenderingMode::DRAW_INSIDE_ONLY);

      // 4. render clipping plane here
      // render_clipping_plane();
    } 
    else // Not CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF
    {
      if (m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF ||
          m_displayMode == DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY)
      {
        // 1. draw solid HALF
        draw_faces_bis(RenderingMode::DRAW_INSIDE_ONLY);

        // 2. render clipping plane here
        // render_clipping_plane();
      } 
      else
      {
        // 1. draw solid FOR ALL
        draw_faces_bis(RenderingMode::DRAW_ALL);
      } 
    }
    glDisable(GL_POLYGON_OFFSET_FILL); 
  }

  inline
  void Basic_viewer::draw_faces_bis(RenderingMode mode)
  {
    m_shaderFace.set_float("u_RenderingMode", static_cast<float>(mode));

    vec4f color = color_to_normalized_vec4(m_facesMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_FACES]);
    glVertexAttrib4fv(2, color.data());
    glDrawArrays(GL_TRIANGLES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_FACES));

    glBindVertexArray(m_vao[VAO_COLORED_FACES]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(2);
    }
    else
    {
      glEnableVertexAttribArray(2);
    }
    glDrawArrays(GL_TRIANGLES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_FACES));

  }

  inline
  void Basic_viewer::draw_rays()
  {
    m_shaderPl.set_float("u_RenderingMode", static_cast<float>(RenderingMode::DRAW_ALL));

    vec4f color = color_to_normalized_vec4(m_raysMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_RAYS]);
    glVertexAttrib4fv(1, color.data());

    glLineWidth(m_sizeRays);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_RAYS));

    glBindVertexArray(m_vao[VAO_COLORED_RAYS]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_RAYS));
  }

  inline
  void Basic_viewer::draw_vertices()
  {
    vec4f color = color_to_normalized_vec4(m_verticeMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_POINTS]);
    glVertexAttrib4fv(1, color.data());
    glDrawArrays(GL_POINTS, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_POINTS));

    glBindVertexArray(m_vao[VAO_COLORED_POINTS]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_POINTS, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_POINTS));
  }

  inline
  void Basic_viewer::draw_lines()
  {
    m_shaderPl.set_float("u_RenderingMode", static_cast<float>(RenderingMode::DRAW_ALL));

    vec4f color = color_to_normalized_vec4(m_linesMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_LINES]);
    glVertexAttrib4fv(1, color.data());
    glLineWidth(m_sizeLines);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_LINES));

    glBindVertexArray(m_vao[VAO_COLORED_LINES]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_LINES));
  }

  inline
  void Basic_viewer::draw_edges()
  {
    glDepthFunc(GL_LEQUAL);

    vec4f color = color_to_normalized_vec4(m_edgesMonoColor);

    glBindVertexArray(m_vao[VAO_MONO_SEGMENTS]);
    glVertexAttrib4fv(1, color.data());
    glLineWidth(m_sizeEdges);
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_MONO_SEGMENTS));

    glBindVertexArray(m_vao[VAO_COLORED_SEGMENTS]);
    if (m_useMonoColor)
    {
      glDisableVertexAttribArray(1);
    }
    else
    {
      glEnableVertexAttribArray(1);
    }
    glDrawArrays(GL_LINES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_SEGMENTS));
  }


  void Basic_viewer::draw_world_axis() 
  {
    glDepthFunc(GL_LEQUAL);

    int &w = m_windowSize.x();
    int &h = m_windowSize.y();

    glViewport(w - w / 5, h - h / 5, w / 5, h / 5);
    m_worldAxisRenderer.draw();

    // Restore the main viewport
    glViewport(0, 0, w, h);
  }

  void Basic_viewer::draw_xy_grid() 
  { 
    glDepthFunc(GL_LEQUAL);

    update_XY_grid_uniforms();
    m_XYGridRenderer.draw();
    update_XY_axis_uniforms();
    m_XYAxisRenderer.draw();
  }

  void Basic_viewer::draw_normals()
  {
    glDepthFunc(GL_LEQUAL);

    glBindVertexArray(m_vao[VAO_COLORED_FACES]);
    glLineWidth(m_sizeNormals);
    glDrawArrays(GL_TRIANGLES, 0, m_scene->number_of_elements(Graphics_scene::POS_COLORED_FACES));
  }

  inline
  void Basic_viewer::initialize_and_load_clipping_plane()
  {
    float size = ((m_scene->bounding_box().xmax() - m_scene->bounding_box().xmin()) +
                  (m_scene->bounding_box().ymax() - m_scene->bounding_box().ymin()) +
                  (m_scene->bounding_box().zmax() - m_scene->bounding_box().zmin()));

    const unsigned int NB_SUBDIVISIONS = 30;

    vec3f color(0,0,0);
    m_clippingPlane.initialize_buffers();
    m_clippingPlane.set_width(0.1f);
    generate_grid(m_clippingPlane, color, size, NB_SUBDIVISIONS);
    m_clippingPlane.load_buffers();

    m_clippingPlane.set_size(m_camera.get_size());  
  }

  inline
  void Basic_viewer::render_clipping_plane()
  {
    if (!m_isOpengl4_3)
    {
      return;
    }

    update_clipping_uniforms();
    if (m_drawClippingPlane)
    {
      m_clippingPlane.draw();
    }
  }

  inline
  void Basic_viewer::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_key_event(key, scancode, action, mods);
  }

  inline
  void Basic_viewer::cursor_callback(GLFWwindow* window, double xpos, double ypo)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    int windowWidth, windowHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    viewer->on_cursor_event(xpos, ypo, windowWidth, windowHeight);
  }

  inline
  void Basic_viewer::mouse_btn_callback(GLFWwindow* window, int button, int action, int mods)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_mouse_button_event(button, action, mods);
  }

  inline
  void Basic_viewer::window_size_callback(GLFWwindow* window, int width, int height)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->m_windowSize = {width, height};

    viewer->m_aspectRatio = static_cast<float>(width) / height;
    glViewport(0, 0, width, height);

    viewer->render_scene(viewer->m_deltaTime);
  }

  inline
  void Basic_viewer::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
  {
    auto viewer = static_cast<Basic_viewer*>(glfwGetWindowUserPointer(window));
    viewer->on_scroll_event(xoffset, yoffset);
  }

  inline 
  void Basic_viewer::print_application_state(float& elapsedTime, const float deltaTime)
  {
    elapsedTime += deltaTime;
    if (elapsedTime * 1000 > 100) // update terminal display each 100ms
    {
      elapsedTime = 0.0f;
      if (m_printApplicationState)
      {
        std::cout << "\33[2K"  
                  << "FPS: "                        << std::round(1 / deltaTime)                               << "\n\33[2K" 
                  << "Camera translation speed: "   << m_camera.get_translation_speed()                        << "    " 
                  << "Camera rotation speed: "      << std::round(m_camera.get_rotation_speed())               << "    "
                  << "Camera constraint axis: "     << m_camera.get_constraint_axis_str()                      << "\n\33[2K"     
                  << "CP translation speed: "       << m_clippingPlane.get_translation_speed()                 << "    "            
                  << "CP rotation speed: "          << std::round(m_clippingPlane.get_rotation_speed())        << "    "     
                  << "CP constraint axis: "         << m_clippingPlane.get_constraint_axis_str()               << "\n\33[2K"     
                  << "Light color: ("               << m_ambientColor.x()                                      << ", " 
                                                    << m_ambientColor.y()                                      << ", " 
                                                    << m_ambientColor.z()                                      << ")\n\33[2K"     
                  << "Size of vertices: "           << m_sizeVertices << "    Size of edges: " <<  m_sizeEdges << "    "            
                  << "\033[F\033[F\033[F\033[F\r" << std::flush;
      }
    }
  }

  inline
  void Basic_viewer::start_action(int action, const float deltaTime)
  {
    switch (action)
    {
    case ROTATE_CLIPPING_PLANE:
    case TRANSLATE_CLIPPING_PLANE:
    case TRANSLATE_CP_ALONG_CAMERA_DIRECTION:
    case TRANSLATE_CAMERA:
    case ROTATE_CAMERA:
      // glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
      break;
    }
  }

  inline 
  void Basic_viewer::exit_app() 
  {
    clear_application();

    std::cout << "\n\n\n\n\n\nAPPLICATION IS CLOSING" << std::endl;
    exit(EXIT_SUCCESS);
  }

    inline
  void Basic_viewer::double_click_event(int btn)
  {
    if (m_camera.is_orbiter()) 
    {
      if (btn == GLFW_MOUSE_BUTTON_RIGHT)
      {
        if (is_key_pressed(m_window, GLFW_KEY_LEFT_CONTROL)) 
        {
          m_camera.reset_orientation();
        }
        m_camera.reset_position();
      }
      else if (btn == GLFW_MOUSE_BUTTON_LEFT)
      {
        if (is_key_pressed(m_window, GLFW_KEY_LEFT_CONTROL))
        {
          m_camera.align_to_plane(m_clippingPlane.get_normal());
        }
        else
        {
          m_camera.align_to_nearest_axis();
        }
      }
      else if (btn == GLFW_MOUSE_BUTTON_MIDDLE)
      {
        m_camera.reset_size();
      }
    }
  }

  inline
  void Basic_viewer::scroll_event(const float deltaTime)
  {
    float yoffset = get_scroll_yOffset() / m_aspectRatio;

    if (is_key_pressed(m_window, GLFW_KEY_LEFT_SHIFT))
    {
      m_camera.increase_zoom_smoothness(yoffset);
    }
    else if (is_key_pressed(m_window, GLFW_KEY_Z) && !m_camera.is_orthographic())
    {
      m_camera.increase_fov(yoffset);
      m_clippingPlane.set_size(m_camera.get_size());
    }
    else if (is_key_pressed(m_window, GLFW_KEY_LEFT_CONTROL) && m_displayMode != DisplayMode::CLIPPING_PLANE_OFF)
    {
        m_clippingPlane.translation(8.f * yoffset * deltaTime);
    }
    else 
    {
      m_camera.move(8.f * yoffset * deltaTime);
      m_clippingPlane.set_size(m_camera.get_size());
    }
  }

  inline 
  void Basic_viewer::change_pivot_point() 
  {
    auto [mouseX, mouseY] = get_mouse_position();

    vec2f nc = utils::normalized_coordinates({mouseX, mouseY}, m_windowSize.x(), m_windowSize.y());

    vec3f cameraPosition = m_camera.get_position(); 
    float cameraX = cameraPosition.x();
    float cameraY = cameraPosition.y();

    float xValue = nc.x() + cameraX; 
    float yValue = nc.y() + cameraY;

    // vec4f test = {cameraPosition.x(), cameraPosition.y(), cameraPosition.z(), 1.0};

    // vec4f vppos = transform::viewport(m_windowSize.x(), m_windowSize.y) * 

    if (utils::inside_bounding_box_2d({xValue, yValue}, {m_boundingBox.first.x(), m_boundingBox.first.y()}, {m_boundingBox.second.x(), m_boundingBox.second.y()})) 
    {
      std::cout << "INSIDE\n";
      m_camera.set_center({nc.x(), nc.y(), 0});
    }
    else 
    {
      m_camera.set_center(utils::center(m_boundingBox.first, m_boundingBox.second));
      std::cout << "OUTSIDE\n";
    }

    // std::cout << "Mouse position : " << nc.x() << " " << nc.y() << ", Camera position : " << cameraPosition.transpose() << "\n";
    // std::cout << "BBOX : " << m_boundingBox.first.transpose() << " " << m_boundingBox.second.transpose() << "\n";
  }

  inline 
  void Basic_viewer::action_event(int action, const float deltaTime)
  {
    switch (action)
    {
    /*APPLICATION*/
    case EXIT: 
      exit_app();
    case PRINT_APPLICATION_STATE:
      m_printApplicationState = !m_printApplicationState; 
      std::cout << "\33[2K" << "\n\33[2K" << "\n\33[2K" << "\n\33[2K" << "\n\33[2K" << "\033[F\033[F\033[F\033[F\r" << std::flush;
      break;
    /*WINDOW*/
    case FULLSCREEN:
      fullscreen();
      break;
    case SCREENSHOT:
      capture_screenshot("./screenshot.png");
      break;
    /*SCENE*/
    case NORMALS_DISPLAY:
      m_drawNormals = !m_drawNormals;
      break;
    case VERTICES_DISPLAY:
      m_drawVertices = !m_drawVertices;
      break;
    case SPHERE_VERTEX_DISPLAY:
      m_drawSphereVertex = !m_drawSphereVertex;
      break;
    case FACES_DISPLAY:
      m_drawFaces = !m_drawFaces;
      break;
    case EDGES_DISPLAY:
      m_drawEdges = !m_drawEdges;
      break;
    case CYLINDER_EDGE_DISPLAY:
      m_drawCylinderEdge = !m_drawCylinderEdge;
      break;
    case SHADING_MODE:
      m_flatShading = !m_flatShading;
      m_areBuffersInitialized = false;
      break;
    case INVERSE_NORMAL:
      m_inverseNormal = !m_inverseNormal;
      m_scene->reverse_all_normals();
      m_areBuffersInitialized = false;
      break;
    case MONO_COLOR:
      m_useMonoColor = !m_useMonoColor;
      break;
    case NORMALS_MONO_COLOR: 
      m_useNormalMonoColor = !m_useNormalMonoColor;
      break;
    case INC_EDGES_SIZE:
      m_sizeEdges = std::min(25.f, m_sizeEdges + 10.0f*deltaTime);
      break;
    case DEC_EDGES_SIZE:
      m_sizeEdges = std::max(0.1f, m_sizeEdges - 10.0f*deltaTime);
      break;
    case INC_POINTS_SIZE:
      m_sizeVertices = std::min(50.f, m_sizeVertices + 10.0f*deltaTime);
      break;
    case DEC_POINTS_SIZE:
      m_sizeVertices = std::max(0.1f, m_sizeVertices - 10.0f*deltaTime);
      break;
    case INC_LIGHT_ALL:
      increase_light_all(deltaTime);
      break;
    case DEC_LIGHT_ALL:
      increase_light_all(-deltaTime);
      break;
    case INC_LIGHT_R:
      increase_red_component(deltaTime);
      break;
    case INC_LIGHT_G:
      increase_green_component(deltaTime);
      break;
    case INC_LIGHT_B:
      increase_blue_component(deltaTime);
      break;
    case DEC_LIGHT_R:
      increase_red_component(-deltaTime);
      break;
    case DEC_LIGHT_G:
      increase_green_component(-deltaTime);
      break;
    case DEC_LIGHT_B:
      increase_blue_component(-deltaTime);
      break;
    case WORLD_AXIS_DISPLAY:
      m_drawWorldAxis = !m_drawWorldAxis;
      break;
    case XY_GRID_DISPLAY:
      m_drawXYGrid = !m_drawXYGrid;
      break;
    /*CAMERA*/
    case UP:
      m_camera.move_up(deltaTime);
      break;
    case DOWN:
      m_camera.move_down(deltaTime);
      break;
    case LEFT:
      m_camera.move_left(deltaTime);
      break;
    case RIGHT:
      m_camera.move_right(deltaTime);
      break;
    case FORWARD:
      m_camera.move(deltaTime);
      break;
    case BACKWARDS:
      m_camera.move(-deltaTime);
      break;
    case INC_CAMERA_ROTATION_SMOOTHNESS:
      m_camera.increase_rotation_smoothness(deltaTime);
      break;
    case DEC_CAMERA_ROTATION_SMOOTHNESS:
      m_camera.decrease_rotation_smoothness(deltaTime);
      break;
    case INC_CAMERA_TRANSLATION_SMOOTHNESS:
      m_camera.increase_translation_smoothness(deltaTime);
      break;
    case DEC_CAMERA_TRANSLATION_SMOOTHNESS:
      m_camera.decrease_translation_smoothness(deltaTime);
      break;
    case SWITCH_CAMERA_MODE:
      m_camera.toggle_mode();
      break;
    case SWITCH_CAMERA_TYPE:
      m_camera.toggle_type();
      break;
    case SWITCH_CAMERA_CONSTRAINT_AXIS: 
      m_camera.switch_constraint_axis();
      break;
    case INC_CAMERA_TRANSLATION_SPEED:
      m_camera.increase_translation_speed(deltaTime);
      break;
    case DEC_CAMERA_TRANSLATION_SPEED:
      m_camera.decrease_translation_speed(deltaTime);
      break;
    case INC_CAMERA_ROTATION_SPEED:
      m_camera.increase_rotation_speed(deltaTime);
      break;
    case DEC_CAMERA_ROTATION_SPEED:
      m_camera.decrease_rotation_speed(deltaTime);
      break;
    case ROTATE_CAMERA:
      rotate_camera();
      break;
    case TRANSLATE_CAMERA:
      translate_camera(deltaTime);
      break;
    case RESET_CAMERA_AND_CP:
      reset_camera_and_clipping_plane();
      break;
    case CHANGE_PIVOT_POINT:
      change_pivot_point();
      break;
    /*CLIPPING PLANE*/
    case INC_CP_TRANSLATION_SPEED:
      m_clippingPlane.increase_translation_speed(deltaTime);
      break;
    case DEC_CP_TRANSLATION_SPEED:
      m_clippingPlane.decrease_translation_speed(deltaTime);
      break;
    case INC_CP_ROTATION_SPEED:
      m_clippingPlane.increase_rotation_speed(deltaTime);
      break;
    case DEC_CP_ROTATION_SPEED:
      m_clippingPlane.decrease_rotation_speed(deltaTime);
      break;
    case SWITCH_CLIPPING_PLANE_DISPLAY:
      m_drawClippingPlane = !m_drawClippingPlane;
      break;
    case SWITCH_CLIPPING_PLANE_MODE:
      switch_display_mode();
      break;
    case ROTATE_CLIPPING_PLANE:
      rotate_clipping_plane();
      break;
    case TRANSLATE_CLIPPING_PLANE:
      translate_clipping_plane(deltaTime);
      break;
    case TRANSLATE_CP_ALONG_CAMERA_DIRECTION:
      translate_clipping_plane(deltaTime, true);
      break;
    case SWITCH_CP_CONSTRAINT_AXIS:
      m_clippingPlane.switch_constraint_axis();
      break;
    case RESET_CLIPPING_PLANE: 
      reset_clipping_plane();
      break;
    /*ANIMATION*/
    case SAVE_KEY_FRAME:
      save_key_frame();
      break;
    case RUN_OR_STOP_ANIMATION:
      run_or_stop_animation();
      break;
    case CLEAR_ANIMATION: 
      m_animationController.clear_buffer();
      break;
    }
  }

  inline 
  void Basic_viewer::increase_red_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_ambientColor.x() += speed;
    if (m_ambientColor.x() > 1.f)
      m_ambientColor.x() = 1.f;
    if (m_ambientColor.x() < 0.f)
      m_ambientColor.x() = 0.f;

    m_diffuseColor.x() += speed;
    if (m_diffuseColor.x() > 1.f)
      m_diffuseColor.x() = 1.f;
    if (m_diffuseColor.x() < 0.f)
      m_diffuseColor.x() = 0.f;

    m_specularColor.x() += speed;
    if (m_specularColor.x() > 1.f)
      m_specularColor.x() = 1.f;
    if (m_specularColor.x() < 0.f)
      m_specularColor.x() = 0.f;
  }

  inline 
  void Basic_viewer::increase_green_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_ambientColor.y() += speed;
    if (m_ambientColor.y() > 1.f)
      m_ambientColor.y() = 1.f;
    if (m_ambientColor.y() < 0.f)
      m_ambientColor.y() = 0.f;

    m_diffuseColor.y() += speed;
    if (m_diffuseColor.y() > 1.f)
      m_diffuseColor.y() = 1.f;
    if (m_diffuseColor.y() < 0.f)
      m_diffuseColor.y() = 0.f;

    m_specularColor.y() += speed;
    if (m_specularColor.y() > 1.f)
      m_specularColor.y() = 1.f;
    if (m_specularColor.y() < 0.f)
      m_specularColor.y() = 0.f;
  }

  inline 
  void Basic_viewer::increase_blue_component(const float deltaTime)
  {
    float speed = deltaTime;
    m_ambientColor.z() += speed;
    if (m_ambientColor.z() > 1.f)
      m_ambientColor.z() = 1.f;
    if (m_ambientColor.z() < 0.f)
      m_ambientColor.z() = 0.f;

    m_diffuseColor.z() += speed;
    if (m_diffuseColor.z() > 1.f)
      m_diffuseColor.z() = 1.f;
    if (m_diffuseColor.z() < 0.f)
      m_diffuseColor.z() = 0.f;

    m_specularColor.z() += speed;
    if (m_specularColor.z() > 1.f)
      m_specularColor.z() = 1.f;
    if (m_specularColor.z() < 0.f)
      m_specularColor.z() = 0.f;
  }


  inline 
  void Basic_viewer::increase_light_all(const float deltaTime)
  {
    increase_red_component(deltaTime);
    increase_green_component(deltaTime);
    increase_blue_component(deltaTime);
  }

  inline 
  void Basic_viewer::save_key_frame() 
  {
    if (m_camera.is_orbiter())
    {
      m_animationController.add_key_frame(
        m_camera.get_position(),
        m_camera.get_orientation()
      );
    }
  }

  inline 
  void Basic_viewer::run_or_stop_animation() 
  {
    if (m_camera.is_orbiter())
    {
      if (m_animationController.is_running()) 
      {
        m_animationController.stop(m_animationController.get_frame());
      }
      else 
      {
        m_animationController.start();
      }
    }
  }

  inline
  void Basic_viewer::reset_camera_and_clipping_plane()
  {
    m_camera.reset_all();
    m_clippingPlane.reset_all();
    m_clippingPlane.set_size(m_camera.get_size());
  }

  inline
  void Basic_viewer::reset_clipping_plane()
  {
    m_clippingPlane.reset_all();
    m_clippingPlane.set_size(m_camera.get_size());
  }

  inline
  void Basic_viewer::switch_display_mode()
  {
    switch(m_displayMode) 
    {
    case DisplayMode::CLIPPING_PLANE_OFF:
      m_displayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_TRANSPARENT_HALF:
      m_displayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_WIRE_HALF: 
      m_displayMode = DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY; 
      break;
    case DisplayMode::CLIPPING_PLANE_SOLID_HALF_ONLY: 
      m_displayMode = DisplayMode::CLIPPING_PLANE_OFF; 
      break;
    }
  }

  inline
  void Basic_viewer::rotate_clipping_plane()
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    m_clippingPlane.set_right_axis(m_camera.get_right());
    m_clippingPlane.set_up_axis(m_camera.get_up());

    m_clippingPlane.rotation(
      mouseDeltaX / m_aspectRatio, 
      mouseDeltaY / m_aspectRatio
    );
  }

  inline
  void Basic_viewer::translate_clipping_plane(const float deltaTime, bool useCameraForward)
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    float deltaX = deltaTime * mouseDeltaX / m_aspectRatio;
    float deltaY = deltaTime * mouseDeltaY / m_aspectRatio;

    if (useCameraForward)
    {
      vec3f forwardDirection = m_camera.get_forward();

      float s = abs(deltaY) > abs(deltaX) ? -deltaY : deltaX;

      m_clippingPlane.translation(forwardDirection, s);
    }
    else 
    {
      m_clippingPlane.translation(-deltaX, deltaY);
    }
  }

  inline
  void Basic_viewer::rotate_camera()
  {
    auto mouseDelta = get_mouse_delta();

    if (m_camera.get_constraint_axis_str() == "Forward")
    {
      mouseDelta = get_roll_mouse_delta();
    }

    float mouseDeltaX = mouseDelta.first;
    float mouseDeltaY = mouseDelta.second; 
    
    m_camera.rotation(
      mouseDeltaX / m_aspectRatio, 
      mouseDeltaY / m_aspectRatio
    );
  }

  inline
  void Basic_viewer::translate_camera(const float deltaTime)
  {
    auto [mouseDeltaX, mouseDeltaY] = get_mouse_delta();

    m_camera.translation(
      deltaTime * -mouseDeltaX / m_aspectRatio,
      deltaTime * mouseDeltaY / m_aspectRatio
    );
  }

  inline
  void Basic_viewer::fullscreen()
  {
    m_isFullscreen = !m_isFullscreen;

    int count;
    GLFWmonitor *monitor = glfwGetMonitors(&count)[0];
    const GLFWvidmode *mode = glfwGetVideoMode(monitor);

    if (m_isFullscreen) 
    {
      m_oldWindowSize = m_windowSize;

#if !defined(GLFW_USE_WAYLAND)
      glfwGetWindowPos(m_window, &m_oldWindowPosition.x(), &m_oldWindowPosition.y()); 
#endif 
      glfwSetWindowMonitor(m_window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
      glViewport(0, 0, mode->width, mode->height);

      m_windowSize.x() = mode->width;
      m_windowSize.y() = mode->height;
    }
    else 
    {
      m_windowSize = m_oldWindowSize;
      glfwSetWindowMonitor(m_window, nullptr, m_oldWindowPosition.x(), m_oldWindowPosition.y(), m_windowSize.x(), m_windowSize.y(), mode->refreshRate);
      glViewport(0, 0, m_windowSize.x(), m_windowSize.y());
    }
  }

  inline
  void Basic_viewer::capture_screenshot(const std::string &filepath)
  {
    // https://lencerf.github.io/post/2019-09-21-save-the-opengl-rendering-to-image-file/ (thanks)
    // https://github.com/nothings/stb/
    // The stb lib used here is from glfw/deps

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_STENCIL_TEST);

    const GLsizei NB_CHANNELS = 4;
    GLsizei stride = NB_CHANNELS * m_windowSize.x();
    stride += (stride % 4) ? (4 - stride % 4) : 0; // stride must be a multiple of 4
    GLsizei bufferSize = stride * m_windowSize.y();

    std::vector<char> buffer(bufferSize);
    m_shaderFace.use(); 
    glPixelStorei(GL_PACK_ALIGNMENT, 4);

    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, m_windowSize.x(), m_windowSize.y(), GL_RGBA, GL_UNSIGNED_BYTE, buffer.data());

    stbi_flip_vertically_on_write(true);
    stbi_write_png(filepath.data(), m_windowSize.x(), m_windowSize.y(), NB_CHANNELS, buffer.data(), stride);
  }
 
  inline
  void Basic_viewer::end_action(int action, const float deltaTime)
  {
    switch (action)
    {
    case TRANSLATE_CAMERA:
    case ROTATE_CAMERA:
      glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
      break;
    }
  }

  inline
  void Basic_viewer::initialize_keys_actions()
  {
    /*APPLICATION*/
    add_keyboard_action({GLFW_KEY_ESCAPE                  }, InputMode::RELEASE, EXIT);
    add_keyboard_action({GLFW_KEY_Q, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, EXIT);
    add_keyboard_action({GLFW_KEY_T                       }, InputMode::RELEASE, PRINT_APPLICATION_STATE);

    /*WINDOW*/
    add_keyboard_action({GLFW_KEY_ENTER, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, FULLSCREEN);
    add_keyboard_action({GLFW_KEY_F2                      }, InputMode::RELEASE, SCREENSHOT);

    /*SCENE*/
    add_keyboard_action({GLFW_KEY_N, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, NORMALS_DISPLAY);
    add_keyboard_action({GLFW_KEY_M, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, NORMALS_MONO_COLOR);

    add_keyboard_action({GLFW_KEY_A}, InputMode::RELEASE, WORLD_AXIS_DISPLAY);
    add_keyboard_action({GLFW_KEY_G}, InputMode::RELEASE, XY_GRID_DISPLAY);

    add_keyboard_action({GLFW_KEY_W}, InputMode::RELEASE, FACES_DISPLAY);
    add_keyboard_action({GLFW_KEY_V}, InputMode::RELEASE, VERTICES_DISPLAY);
    add_keyboard_action({GLFW_KEY_E}, InputMode::RELEASE, EDGES_DISPLAY);

    add_keyboard_action({GLFW_KEY_S}, InputMode::RELEASE, SHADING_MODE);
    add_keyboard_action({GLFW_KEY_N}, InputMode::RELEASE, INVERSE_NORMAL);
    add_keyboard_action({GLFW_KEY_M}, InputMode::RELEASE, MONO_COLOR);

    add_keyboard_action({GLFW_KEY_E, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, CYLINDER_EDGE_DISPLAY);
    add_keyboard_action({GLFW_KEY_V, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, SPHERE_VERTEX_DISPLAY);

    add_keyboard_action({GLFW_KEY_EQUAL, GLFW_KEY_LEFT_SHIFT, GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, INC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_KP_ADD, GLFW_KEY_LEFT_CONTROL                    }, InputMode::HOLD, INC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_6, GLFW_KEY_LEFT_CONTROL                         }, InputMode::HOLD, DEC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_KP_SUBTRACT, GLFW_KEY_LEFT_CONTROL               }, InputMode::HOLD, DEC_POINTS_SIZE);
    add_keyboard_action({GLFW_KEY_EQUAL, GLFW_KEY_LEFT_SHIFT                       }, InputMode::HOLD, INC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_KP_ADD                                           }, InputMode::HOLD, INC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_6                                                }, InputMode::HOLD, DEC_EDGES_SIZE);
    add_keyboard_action({GLFW_KEY_KP_SUBTRACT                                      }, InputMode::HOLD, DEC_EDGES_SIZE);

    add_keyboard_action({GLFW_KEY_PAGE_UP                          }, InputMode::HOLD, INC_LIGHT_ALL);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN                        }, InputMode::HOLD, DEC_LIGHT_ALL);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_SHIFT  }, InputMode::HOLD, INC_LIGHT_R);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_SHIFT  }, InputMode::HOLD, DEC_LIGHT_R);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_ALT    }, InputMode::HOLD, INC_LIGHT_G);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_ALT    }, InputMode::HOLD, DEC_LIGHT_G);
    add_keyboard_action({GLFW_KEY_PAGE_UP,    GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, INC_LIGHT_B);
    add_keyboard_action({GLFW_KEY_PAGE_DOWN,  GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, DEC_LIGHT_B);

    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, RESET_CAMERA_AND_CP);

    /*CAMERA*/
    add_keyboard_action({GLFW_KEY_UP,   GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, FORWARD);
    add_keyboard_action({GLFW_KEY_DOWN, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, BACKWARDS);

    add_keyboard_action({GLFW_KEY_UP   }, InputMode::HOLD, UP);
    add_keyboard_action({GLFW_KEY_DOWN }, InputMode::HOLD, DOWN);
    add_keyboard_action({GLFW_KEY_LEFT }, InputMode::HOLD, LEFT);
    add_keyboard_action({GLFW_KEY_RIGHT}, InputMode::HOLD, RIGHT);

    
    add_keyboard_action({GLFW_KEY_SPACE}, InputMode::RELEASE, SWITCH_CAMERA_TYPE);
    add_keyboard_action({GLFW_KEY_O    }, InputMode::RELEASE, SWITCH_CAMERA_MODE); 

    add_keyboard_action({GLFW_KEY_A, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, SWITCH_CAMERA_CONSTRAINT_AXIS);

    add_keyboard_action({GLFW_KEY_X                     }, InputMode::HOLD, INC_CAMERA_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R                     }, InputMode::HOLD, INC_CAMERA_ROTATION_SPEED);
    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CAMERA_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CAMERA_ROTATION_SPEED);

    add_mouse_action({GLFW_MOUSE_BUTTON_LEFT }, InputMode::HOLD, ROTATE_CAMERA);
    add_mouse_action({GLFW_MOUSE_BUTTON_RIGHT}, InputMode::HOLD, TRANSLATE_CAMERA);

    add_keyboard_action({GLFW_MOUSE_BUTTON_RIGHT, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, CHANGE_PIVOT_POINT);

    /*CLIPPING PLANE*/
    add_keyboard_action({GLFW_KEY_C, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, SWITCH_CLIPPING_PLANE_DISPLAY);
    add_keyboard_action({GLFW_KEY_C                   }, InputMode::RELEASE, SWITCH_CLIPPING_PLANE_MODE);

    add_keyboard_action({GLFW_KEY_A, GLFW_KEY_LEFT_CONTROL}, InputMode::RELEASE, SWITCH_CP_CONSTRAINT_AXIS);

    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_CONTROL                     }, InputMode::HOLD, INC_CP_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_X, GLFW_KEY_LEFT_CONTROL, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CP_TRANSLATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL                     }, InputMode::HOLD, INC_CP_ROTATION_SPEED);
    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_LEFT_CONTROL, GLFW_KEY_LEFT_SHIFT}, InputMode::HOLD, DEC_CP_ROTATION_SPEED);

    add_mouse_action({GLFW_MOUSE_BUTTON_LEFT,   GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, ROTATE_CLIPPING_PLANE);
    add_mouse_action({GLFW_MOUSE_BUTTON_RIGHT,  GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, TRANSLATE_CLIPPING_PLANE);
    add_mouse_action({GLFW_MOUSE_BUTTON_MIDDLE, GLFW_KEY_LEFT_CONTROL}, InputMode::HOLD, TRANSLATE_CP_ALONG_CAMERA_DIRECTION);

    add_keyboard_action({GLFW_KEY_R, GLFW_KEY_TAB}, InputMode::RELEASE, RESET_CLIPPING_PLANE);

    /*ANIMATION*/
    add_keyboard_action({GLFW_KEY_F1                               }, InputMode::RELEASE, RUN_OR_STOP_ANIMATION);
    add_keyboard_action({GLFW_KEY_F1, GLFW_KEY_LEFT_ALT            }, InputMode::RELEASE, SAVE_KEY_FRAME);
    add_keyboard_action({GLFW_KEY_F1, GLFW_KEY_D, GLFW_KEY_LEFT_ALT}, InputMode::RELEASE, CLEAR_ANIMATION);

    /*===================== BIND DESCRIPTIONS ============================*/

    set_action_description({
      {"Animation", {
        {get_binding_text_from_action(SAVE_KEY_FRAME),        "Add a key frame to animation (only available in orbiter)"},
        {get_binding_text_from_action(CLEAR_ANIMATION),       "Delete animation (only available in orbiter)"},
        {get_binding_text_from_action(RUN_OR_STOP_ANIMATION), "Start/Stop animation (only available in orbiter)"},
      }},
      {"Clipping plane", {
        {get_binding_text_from_action(SWITCH_CLIPPING_PLANE_MODE),    "Switch clipping plane display mode"},
        {get_binding_text_from_action(SWITCH_CLIPPING_PLANE_DISPLAY), "Toggle clipping plane rendering on/off"},
        {get_binding_text_from_action(SWITCH_CP_CONSTRAINT_AXIS),     "Switch constraint axis for clipping plane rotation"},

        {"", ""},

        {get_binding_text_from_action(INC_CP_ROTATION_SPEED),    "Increase rotation speed"},
        {get_binding_text_from_action(DEC_CP_ROTATION_SPEED),    "Decrease rotation speed"},
        {get_binding_text_from_action(INC_CP_TRANSLATION_SPEED), "Increase translation speed"},
        {get_binding_text_from_action(DEC_CP_TRANSLATION_SPEED), "Decrease translation speed"},

        {"", ""},

        {get_binding_text_from_action(ROTATE_CLIPPING_PLANE),               "Rotate the clipping plane when enabled"},
        {get_binding_text_from_action(TRANSLATE_CLIPPING_PLANE),            "Translate the clipping plane when enabled"},
        {get_binding_text_from_action(TRANSLATE_CP_ALONG_CAMERA_DIRECTION), "Translate the clipping plane along camera direction axis when enabled"},
        {"[LCTRL+WHEEL]",                                                   "Translate the clipping plane along its normal when enabled"},

        {"", ""},

        {get_binding_text_from_action(RESET_CLIPPING_PLANE), "Reset clipping plane"},
      }},
      {"Camera", {
        {get_binding_text_from_action(FORWARD),                      "Move camera forward"},
        {get_binding_text_from_action(BACKWARDS),                    "Move camera backwards"},
        {get_binding_text_from_action(UP),                           "Move camera up"},
        {get_binding_text_from_action(DOWN),                         "Move camera down"},
        {get_binding_text_from_action(RIGHT),                        "Move camera right"},
        {get_binding_text_from_action(LEFT),                         "Move camera left"},
        {"",                                                         ""},
        {get_binding_text_from_action(SWITCH_CAMERA_MODE),           "Switch to Perspective/Orthographic view"},
        {get_binding_text_from_action(SWITCH_CAMERA_TYPE),           "Switch to Orbiter/Free-fly camera type"},
        {get_binding_text_from_action(SWITCH_CAMERA_CONSTRAINT_AXIS),"Switch constraint axis for camera rotation"},
        {"",                                                         ""},
        {get_binding_text_from_action(ROTATE_CAMERA),                "Rotate the camera"},
        {get_binding_text_from_action(TRANSLATE_CAMERA),             "Translate the camera"},
        {get_binding_text_from_action(RESET_CAMERA_AND_CP),                    "Reset camera"},
        {"",                                                         ""},
        {get_binding_text_from_action(INC_CAMERA_ROTATION_SPEED),    "Increase rotation speed"},
        {get_binding_text_from_action(DEC_CAMERA_ROTATION_SPEED),    "Decrease rotation speed"},
        {get_binding_text_from_action(INC_CAMERA_TRANSLATION_SPEED), "Increase translation speed"},
        {get_binding_text_from_action(DEC_CAMERA_TRANSLATION_SPEED), "Decrease translation speed"},
        {"",                                                         ""},
        {"[WHEEL]",                                                  "Zoom in/out"},
        {"[LEFT_DOUBLE_CLICK]",                                      "Aligns camera to nearest axis"},
        {"[LCTRL+LEFT_DOUBLE_CLICK]",                                "Aligns camera to clipping plane"},
        {"[RIGHT_DOUBLE_CLICK]",                                     "Center the camera to the object"},
        {"[LCTRL+RIGHT_DOUBLE_CLICK]",                               "Reset camera position & orientation"},
        {get_binding_text_from_action(CHANGE_PIVOT_POINT),           "Change pivot point (center of the scene)"},
        {"[Z+WHEEL]",                                                "Increase/Decrease FOV"},
      }},
      {"Scene", {
        {get_binding_text_from_action(INC_LIGHT_ALL),         "Increase light (all colors, use shift/alt/ctrl for one rgb component)"},
        {get_binding_text_from_action(DEC_LIGHT_ALL),         "Decrease light (all colors, use shift/alt/ctrl for one rgb component)"},
        {"",                                                  ""},
        {get_binding_text_from_action(NORMALS_DISPLAY),       "Toggle normals display"},
        {get_binding_text_from_action(VERTICES_DISPLAY),      "Toggle vertices display"},
        {get_binding_text_from_action(SPHERE_VERTEX_DISPLAY), "Toggle vertices display as sphere"},
        {get_binding_text_from_action(EDGES_DISPLAY),         "Toggle edges display"},
        {get_binding_text_from_action(CYLINDER_EDGE_DISPLAY), "Toggle edges display as cylinder"},
        {get_binding_text_from_action(FACES_DISPLAY),         "Toggle faces display"},
        {"",                                                  ""},
        {get_binding_text_from_action(WORLD_AXIS_DISPLAY),    "Toggle world axis display"},
        {get_binding_text_from_action(XY_GRID_DISPLAY),       "Toggle XY grid display"},
        {"",                                                  ""},
        {get_binding_text_from_action(INC_POINTS_SIZE),       "Increase size of vertices"},
        {get_binding_text_from_action(DEC_POINTS_SIZE),       "Decrease size of vertices"},
        {get_binding_text_from_action(INC_EDGES_SIZE),        "Increase size of edges"},
        {get_binding_text_from_action(DEC_EDGES_SIZE),        "Decrease size of edges"},
        {"",                                                  ""},
        {get_binding_text_from_action(MONO_COLOR),            "Toggle mono color"},
        {get_binding_text_from_action(NORMALS_MONO_COLOR),    "Toggle normals mono color"},
        {get_binding_text_from_action(INVERSE_NORMAL),        "Invert direction of normals"},
        {get_binding_text_from_action(SHADING_MODE),          "Switch between flat/Gouraud shading display"},
        {"",                                                  ""},
        {"[MIDDLE_DOUBLE_CLICK]",                             "Show entire scene"},
      }},         
      {"Window", {
        {get_binding_text_from_action(FULLSCREEN), "Switch to windowed/fullscreen mode"},
        {get_binding_text_from_action(SCREENSHOT), "Take a screenshot of the current view"},
      }},
      {"Application", {
        {get_binding_text_from_action(EXIT),                    "Exit program"},
        {"",                                                    ""},
        {get_binding_text_from_action(PRINT_APPLICATION_STATE), "Activate/Deactivate application state refreshing (FPS...)"},
      }},
    });

    print_help();
  }
} // end namespace GLFW 
} // end namespace CGAL

#endif // CGAL_BASIC_VIEWER_GLFW_IMPL_H