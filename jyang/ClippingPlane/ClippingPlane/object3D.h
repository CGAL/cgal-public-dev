#pragma once
#include <qmatrix4x4.h>
#include <qopenglbuffer.h>
#include <qopengltexture.h>
#include <qopenglfunctions.h>
#include <qopenglshaderprogram.h>

struct Vertex {
	Vertex() {};
	Vertex(QVector3D position, QVector2D texCoord, QVector3D normal) :
		position(position), texCoord(texCoord), normal(normal) {};
	QVector3D position;
	QVector2D texCoord;
	QVector3D normal;
};

class Object3D
{
public:
	Object3D(const QVector<Vertex>& vertices, const QVector<GLuint>& indices, const QImage& image);
	~Object3D();

	void init(const QVector<Vertex>& vertices, const QVector<GLuint>& indices, const QImage& image);
	void draw(QOpenGLShaderProgram* shaderProgram, QOpenGLFunctions* functions);
	void translate(const QVector3D& translate);

private:
	QMatrix4x4 mMatrix;
	QOpenGLTexture* texture;
	QOpenGLBuffer vertexBuffer;
	QOpenGLBuffer indexBuffer;
};