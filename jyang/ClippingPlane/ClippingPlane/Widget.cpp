#include "Widget.h"

Widget::Widget(QWidget* parent)
	: QOpenGLWidget(parent), texture(0), indexBuffer(QOpenGLBuffer::IndexBuffer)
{
	init();
}

Widget::~Widget()
{
}

void Widget::init()
{
	initGL();
	initShaders();
	//initGeometry(1.0);
}

void Widget::initGL()
{
	// clear the screen with black
	//glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	//glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);
}

void Widget::resizeGL(int w, int h)
{
	float aspect = w / (float)h;

	pMatrix.setToIdentity();
	pMatrix.perspective(45, aspect, 0.1f, 10.0f);
}

void Widget::paintGL()
{
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set transformation for model
	mMatrix.setToIdentity();
	mMatrix.translate(0.0, 0.0, -5.0);
	mMatrix.rotate(30, 1.0, 0.0, 0.0);
	mMatrix.rotate(30, 0.0, 1.0, 0.0);

	//texture->bind();

	//shaderProgram.bind();
	//shaderProgram.setUniformValue("qt_modelViewProjectionMatrix", pMatrix * mMatrix);
	//shaderProgram.setUniformValue("qt_texture", 0);

	//int offset = 0;

	//int verLoc = shaderProgram.attributeLocation("qt_vertex");
	//shaderProgram.enableAttributeArray(verLoc);
	//shaderProgram.setAttributeBuffer(verLoc, GL_FLOAT, offset, 3, sizeof(Vertex));

	//int texLoc = shaderProgram.attributeLocation("qt_texCoord");
	//shaderProgram.enableAttributeArray(texLoc);
	//shaderProgram.setAttributeBuffer(texLoc, GL_FLOAT, offset, 2, sizeof(Vertex));

	//indexBuffer.bind();

	//glDrawElements(GL_TRIANGLES, indexBuffer.size(), GL_UNSIGNED_INT, 0);
}

void Widget::initShaders()
{
	if (!shaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "./VertexShader.glsl")) {
		close();
	}
	if (!shaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "./FragmentShader.fsh")) {
		qDebug() << shaderProgram.log();
		close();
	}
	if (!shaderProgram.link()) {
		qDebug() << shaderProgram.log();
		close();
	}
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
	QPainter painter(&image);
	painter.fillRect(image.rect(), Qt::green);

	//texture = new QOpenGLTexture(image);

	// Set nearest filtering mode for texture minification
	//texture->setMinificationFilter(QOpenGLTexture::Linear);

	// Set bilinear filtering mode for texture magnification
	//texture->setMagnificationFilter(QOpenGLTexture::Linear);

	// Wrap texture coordinates by repreating
	// f. ex. texture coordinate (1.1, 1.2) is same as 0.1, 0.2;
	//texture->setWrapMode(QOpenGLTexture::Repeat);
}
