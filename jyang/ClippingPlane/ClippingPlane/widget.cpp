#include "widget.h"

Widget::Widget(QWidget* parent): 
	QOpenGLWidget(parent), indexBuffer(QOpenGLBuffer::IndexBuffer)
{
}

Widget::~Widget()
{
}

void Widget::initializeGL()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	initShader();
	initGeometry(1.0);
}

void Widget::resizeGL(int w, int h)
{
	float aspect = w / (float)h;

	pMatrix.setToIdentity();
	pMatrix.perspective(45, aspect, 0.1f, 10.0f);
}

void Widget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	QMatrix4x4 mMatrix;
	mMatrix.setToIdentity();
	mMatrix.translate(0.0, 0.0, -5.0);
	mMatrix.rotate(30, 1.0, 0.0, 0.0);
	mMatrix.rotate(30, 0.0, 1.0, 0.0);

	texture->bind();

	shaderProgram.bind();
	shaderProgram.setUniformValue("qt_modelViewProjectionMatrix", pMatrix * mMatrix);
	shaderProgram.setUniformValue("qt_texture0", 0);

	arrayBuffer.bind();

	int offset = 0;

	int verLoc = shaderProgram.attributeLocation("qt_vertex");
	shaderProgram.enableAttributeArray(verLoc);
	shaderProgram.setAttributeBuffer(verLoc, GL_FLOAT, offset, 3, sizeof(Vertex));

	offset += sizeof(QVector3D);

	int texLoc = shaderProgram.attributeLocation("qt_texCoord");
	shaderProgram.enableAttributeArray(texLoc);
	shaderProgram.setAttributeBuffer(texLoc, GL_FLOAT, offset, 2, sizeof(Vertex));

	indexBuffer.bind();

	glDrawElements(GL_TRIANGLES, indexBuffer.size(), GL_UNSIGNED_INT, 0);
}

void Widget::initShader()
{
	if (!shaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "./vertex.glsl"))
		close();
	if (!shaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "./fragment.glsl"))
		close();
	if (!shaderProgram.link())
		close();
}

void Widget::initGeometry(float w)
{
	QVector<Vertex> vertices;
	vertices <<
		Vertex(QVector3D(-w, w, w), QVector2D(0.0, 1.0), QVector3D(0.0, 0.0, 1.0)) <<
		Vertex(QVector3D(-w, -w, w), QVector2D(0.0, 0.0), QVector3D(0.0, 0.0, 1.0)) <<
		Vertex(QVector3D(w, w, w), QVector2D(1.0, 1.0), QVector3D(0.0, 0.0, 1.0)) <<
		Vertex(QVector3D(w, -w, w), QVector2D(1.0, 0.0), QVector3D(0.0, 0.0, 1.0)) <<

		Vertex(QVector3D(w, w, w), QVector2D(0.0, 1.0), QVector3D(1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(w, -w, w), QVector2D(0.0, 0.0), QVector3D(1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(w, w, -w), QVector2D(1.0, 1.0), QVector3D(1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(w, -w, -w), QVector2D(1.0, 0.0), QVector3D(1.0, 0.0, 0.0)) <<

		Vertex(QVector3D(w, w, w), QVector2D(0.0, 1.0), QVector3D(0.0, 1.0, 0.0)) <<
		Vertex(QVector3D(w, w, -w), QVector2D(0.0, 0.0), QVector3D(0.0, 1.0, 0.0)) <<
		Vertex(QVector3D(-w, w, w), QVector2D(1.0, 1.0), QVector3D(0.0, 1.0, 0.0)) <<
		Vertex(QVector3D(-w, w, -w), QVector2D(1.0, 0.0), QVector3D(0.0, 1.0, 0.0)) <<

		Vertex(QVector3D(w, w, -w), QVector2D(0.0, 1.0), QVector3D(0.0, 0.0, -1.0)) <<
		Vertex(QVector3D(w, -w, -w), QVector2D(0.0, 0.0), QVector3D(0.0, 0.0, -1.0)) <<
		Vertex(QVector3D(-w, w, -w), QVector2D(1.0, 1.0), QVector3D(0.0, 0.0, -1.0)) <<
		Vertex(QVector3D(-w, -w, -w), QVector2D(1.0, 0.0), QVector3D(0.0, 0.0, -1.0)) <<

		Vertex(QVector3D(-w, w, w), QVector2D(0.0, 1.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(-w, w, -w), QVector2D(0.0, 0.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(-w, -w, w), QVector2D(1.0, 1.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(-w, -w, -w), QVector2D(1.0, 0.0), QVector3D(-1.0, 0.0, 0.0)) <<

		Vertex(QVector3D(-w, -w, w), QVector2D(0.0, 1.0), QVector3D(0.0, -1.0, 0.0)) <<
		Vertex(QVector3D(-w, -w, -w), QVector2D(0.0, 0.0), QVector3D(0.0, -1.0, 0.0)) <<
		Vertex(QVector3D(w, -w, w), QVector2D(1.0, 1.0), QVector3D(0.0, -1.0, 0.0)) <<
		Vertex(QVector3D(w, -w, -w), QVector2D(1.0, 0.0), QVector3D(0.0, -1.0, 0.0));

	QVector<GLuint> indices;
	indices << 0 << 1 << 2 << 2 << 1 << 3;
	for (int i = 0; i < 24; i += 4)
		indices << i + 0 << i + 1 << i + 2 << i + 2 << i + 1 << i + 3;

	arrayBuffer.create();
	arrayBuffer.bind();
	arrayBuffer.allocate(vertices.constData(), vertices.size() * sizeof(Vertex));
	arrayBuffer.release();

	indexBuffer.create();
	indexBuffer.bind();
	indexBuffer.allocate(indices.constData(), indices.size() * sizeof(GLuint));
	indexBuffer.release();

	QImage image(100, 100, QImage::Format_ARGB32_Premultiplied);
	image.fill(QColor(255, 0, 0, 255));

	texture = new QOpenGLTexture(image);

	// Set nearest filtering mode for texture minification
	texture->setMinificationFilter(QOpenGLTexture::Linear);

	// Set bilinear filtering mode for texture magnification
	texture->setMagnificationFilter(QOpenGLTexture::Linear);

	// Wrap texture coordinates by repreating
	// f. ex. texture coordinate (1.1, 1.2) is same as 0.1, 0.2;
	texture->setWrapMode(QOpenGLTexture::Repeat);
}
