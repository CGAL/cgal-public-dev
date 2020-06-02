#include "object3D.h"

Object3D::Object3D(const QVector<Vertex>& vertices, const QVector<GLuint>& indices, const QImage& image) :
	indexBuffer(QOpenGLBuffer::IndexBuffer), texture(0)
{
	init(vertices, indices, image);
}

Object3D::~Object3D()
{
}

void Object3D::init(const QVector<Vertex>& vertices, const QVector<GLuint>& indices, const QImage& image)
{
	if (vertexBuffer.isCreated()) vertexBuffer.destroy();
	if (indexBuffer.isCreated()) indexBuffer.destroy();

	vertexBuffer.create();
	vertexBuffer.bind();
	vertexBuffer.allocate(vertices.constData(), vertices.size() * sizeof(Vertex));
	vertexBuffer.release();

	indexBuffer.create();
	indexBuffer.bind();
	indexBuffer.allocate(indices.constData(), indices.size() * sizeof(GLuint));
	indexBuffer.release();

	texture = new QOpenGLTexture(image.mirrored());
	texture->setMinificationFilter(QOpenGLTexture::Linear);
	texture->setMagnificationFilter(QOpenGLTexture::Linear);
	texture->setWrapMode(QOpenGLTexture::Repeat);

	// model matrix
	mMatrix.setToIdentity();
}

void Object3D::draw(QOpenGLShaderProgram* shaderProgram, QOpenGLFunctions* functions)
{
	if (!vertexBuffer.isCreated() || !indexBuffer.isCreated()) return;

	texture->bind(0);
	shaderProgram->setUniformValue("u_texture", 0);
	shaderProgram->setUniformValue("u_mMatrix", mMatrix);

	vertexBuffer.bind();

	int offset = 0;

	int verLoc = shaderProgram->attributeLocation("a_position");
	shaderProgram->enableAttributeArray(verLoc);
	shaderProgram->setAttributeBuffer(verLoc, GL_FLOAT, offset, 3, sizeof(Vertex));

	offset += sizeof(QVector3D);

	int texLoc = shaderProgram->attributeLocation("a_texcoord");
	shaderProgram->enableAttributeArray(texLoc);
	shaderProgram->setAttributeBuffer(texLoc, GL_FLOAT, offset, 2, sizeof(Vertex));

	offset += sizeof(QVector2D);

	int normLoc = shaderProgram->attributeLocation("a_normal");
	shaderProgram->enableAttributeArray(normLoc);
	shaderProgram->setAttributeBuffer(normLoc, GL_FLOAT, offset, 3, sizeof(Vertex));

	indexBuffer.bind();

	functions->glDrawElements(GL_TRIANGLES, indexBuffer.size(), GL_UNSIGNED_INT, 0);

	vertexBuffer.release();
	indexBuffer.release();
	texture->release();
}

void Object3D::translate(const QVector3D& translate)
{
	mMatrix.translate(translate);
}
