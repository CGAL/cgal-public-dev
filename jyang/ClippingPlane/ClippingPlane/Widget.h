#pragma once
#include <qopenglwidget.h>
#include <qmatrix4x4.h>
#include <qopenglshaderprogram.h>
#include <qopengltexture.h>
#include <qopenglbuffer.h>
#include <qpainter.h>
#include <qdebug.h>
#include <iostream>

class Widget :
	public QOpenGLWidget
{
public:
	Widget(QWidget* parent = 0);
	~Widget();

protected:
	void init();

	void initGL();
	void resizeGL(int w, int h);
	void paintGL();

	void initShaders();
	void initGeometry(float w);

private:
	QMatrix4x4 pMatrix;
	QMatrix4x4 mMatrix;
	QOpenGLShaderProgram shaderProgram;
	QOpenGLTexture* texture;
	QOpenGLBuffer arrayBuffer;
	QOpenGLBuffer indexBuffer;
};

struct Vertex {
	Vertex() {};
	Vertex(QVector3D position, QVector2D texCoord, QVector3D normal) :
		position(position), texCoord(texCoord), normal(normal) {};
	QVector3D position;
	QVector2D texCoord;
	QVector3D normal;
};