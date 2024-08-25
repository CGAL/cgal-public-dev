#pragma once

/*************GLFW WINDOW PARAMS*************/

#ifndef CGAL_WINDOW_WIDTH_INIT
#define CGAL_WINDOW_WIDTH_INIT 500
#endif

#ifndef CGAL_WINDOW_HEIGHT_INIT
#define CGAL_WINDOW_HEIGHT_INIT 450
#endif

#ifndef CGAL_WINDOW_SAMPLES
#define CGAL_WINDOW_SAMPLES 1
#endif

/*************GLFW SCENE PARAMS*************/
#ifndef CGAL_size_vertices
#define CGAL_size_vertices 7.0f
#endif

#ifndef CGAL_SIZE_EDGES
#define CGAL_SIZE_EDGES 1.1f
#endif

#ifndef CGAL_SIZE_RAYS
#define CGAL_SIZE_RAYS 3.1f
#endif

#ifndef CGAL_SIZE_LINES
#define CGAL_SIZE_LINES 3.1f
#endif

#ifndef CGAL_SIZE_NORMALS
#define CGAL_SIZE_NORMALS 0.2f
#endif 

#ifndef CGAL_NORMAL_HEIGHT_FACTOR
#define CGAL_NORMAL_HEIGHT_FACTOR 0.02f
#endif

#ifndef CGAL_NORMALS_MONO_COLOR
#define CGAL_NORMALS_MONO_COLOR { 220, 20, 20 }
#endif

#ifndef CGAL_LIGHT_POSITION
#define CGAL_LIGHT_POSITION { 0.0f, 0.0f, 0.0f, 0.0f }
#endif

#ifndef CGAL_AMBIENT_COLOR
#define CGAL_AMBIENT_COLOR { 0.6f, 0.5f, 0.5f, 0.5f }
#endif

#ifndef CGAL_DIFFUSE_COLOR
#define CGAL_DIFFUSE_COLOR { 0.9f, 0.9f, 0.9f, 0.9f }
#endif

#ifndef CGAL_SPECULAR_COLOR
#define CGAL_SPECULAR_COLOR { 0.0f, 0.0f, 0.0f, 1.0f }
#endif

#ifndef CGAL_SHININESS
#define CGAL_SHININESS 1.0f
#endif

/************GLFW CAMERA PARAMS*************/
#ifndef CGAL_CAMERA_TRANSLATION_SPEED
#define CGAL_CAMERA_TRANSLATION_SPEED 1.f
#endif

#ifndef CGAL_CAMERA_ROTATION_SPEED
#define CGAL_CAMERA_ROTATION_SPEED 270.f
#endif

#ifndef CGAL_CAMERA_RADIUS
#define CGAL_CAMERA_RADIUS 5.f
#endif

#ifndef CGAL_CAMERA_FOV
#define CGAL_CAMERA_FOV 45.f
#endif

#ifndef CGAL_CAMERA_TRANSLATION_SMOOTHNESS
#define CGAL_CAMERA_TRANSLATION_SMOOTHNESS 0.25f
#endif

#ifndef CGAL_CAMERA_ROTATION_SMOOTHNESS
#define CGAL_CAMERA_ROTATION_SMOOTHNESS 0.25f
#endif

#ifndef CGAL_CAMERA_ZOOM_SMOOTHNESS
#define CGAL_CAMERA_ZOOM_SMOOTHNESS 0.1f
#endif

/*********GLFW CLIPPING PLANE PARAMS********/
#ifndef CGAL_CLIPPING_PLANE_RENDERING_TRANSPARENCY
#define CGAL_CLIPPING_PLANE_RENDERING_TRANSPARENCY 0.5f
#endif

#ifndef CGAL_CLIPPING_PLANE_TRANSLATION_SPEED
#define CGAL_CLIPPING_PLANE_TRANSLATION_SPEED 1.f
#endif

#ifndef CGAL_CLIPPING_PLANE_ROTATION_SPEED
#define CGAL_CLIPPING_PLANE_ROTATION_SPEED 270.f
#endif

#ifndef CGAL_CLIPPING_PLANE_TRANSLATION_SMOOTHNESS
#define CGAL_CLIPPING_PLANE_TRANSLATION_SMOOTHNESS 1.f
#endif

#ifndef CGAL_CLIPPING_PLANE_ROTATION_SMOOTHNESS
#define CGAL_CLIPPING_PLANE_ROTATION_SMOOTHNESS 0.35f
#endif

/***********GLFW ANIMATION PARAMS***********/
