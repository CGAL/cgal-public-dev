#pragma once
#include <qopenglwidget.h>
#include <qmatrix4x4.h>
#include <qopenglshaderprogram.h>
#include <qopenglbuffer.h>
#include <qopengltexture.h>

class Widget :
	public QOpenGLWidget
{
public:
	Widget(QWidget* parent);
	~Widget();

protected:
	void initializeGL();
	void resizeGL(int w, int h);
	void paintGL();

	void initShader();
	void initGeometry(float w);

private:
	QMatrix4x4 mMatrix;
	QMatrix4x4 pMatrix;
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