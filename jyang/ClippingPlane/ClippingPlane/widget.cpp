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
	//glEnable(GL_CULL_FACE);
	glEnable(GL_BLEND);

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

	// model matrix
	mMatrix.setToIdentity();
	mMatrix.translate(0.0, 0.0, -5.0);
	mMatrix.rotate(30, 1.0, 0.0, 0.0);
	mMatrix.rotate(30, 0.0, 1.0, 0.0);
	mMatrix.rotate(rotation);

	// view matrix
	vMatrix.setToIdentity();

	// bind texture
	texture->bind();

	shaderProgram.bind();
	shaderProgram.setUniformValue("u_mMatrix", mMatrix);
	shaderProgram.setUniformValue("u_pMatrix", pMatrix);
	shaderProgram.setUniformValue("u_vMatrix", vMatrix);
	shaderProgram.setUniformValue("u_texture", 0);

	// bind buffer
	arrayBuffer.bind();

	int offset = 0;
	int verLoc = shaderProgram.attributeLocation("a_position");
	shaderProgram.enableAttributeArray(verLoc);
	shaderProgram.setAttributeBuffer(verLoc, GL_FLOAT, offset, 3, sizeof(Vertex));

	offset += sizeof(QVector3D);
	int texLoc = shaderProgram.attributeLocation("a_texcoord");
	shaderProgram.enableAttributeArray(texLoc);
	shaderProgram.setAttributeBuffer(texLoc, GL_FLOAT, offset, 2, sizeof(Vertex));

	offset += sizeof(QVector2D);
	int normLoc = shaderProgram.attributeLocation("a_normal");
	shaderProgram.enableAttributeArray(normLoc);
	shaderProgram.setAttributeBuffer(normLoc, GL_FLOAT, offset, 3, sizeof(Vertex));

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
		// front
		Vertex(QVector3D(-w, w, w), QVector2D(0.0, 1.0), QVector3D(0.0, 0.0, 1.0)) <<
		Vertex(QVector3D(-w, -w, w), QVector2D(0.0, 0.0), QVector3D(0.0, 0.0, 1.0)) <<
		Vertex(QVector3D(w, w, w), QVector2D(1.0, 1.0), QVector3D(0.0, 0.0, 1.0)) <<
		Vertex(QVector3D(w, -w, w), QVector2D(1.0, 0.0), QVector3D(0.0, 0.0, 1.0)) <<

		// right
		Vertex(QVector3D(w, w, w), QVector2D(0.0, 1.0), QVector3D(1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(w, -w, w), QVector2D(0.0, 0.0), QVector3D(1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(w, w, -w), QVector2D(1.0, 1.0), QVector3D(1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(w, -w, -w), QVector2D(1.0, 0.0), QVector3D(1.0, 0.0, 0.0)) <<

		// top
		Vertex(QVector3D(w, w, w), QVector2D(0.0, 1.0), QVector3D(0.0, 1.0, 0.0)) <<
		Vertex(QVector3D(w, w, -w), QVector2D(0.0, 0.0), QVector3D(0.0, 1.0, 0.0)) <<
		Vertex(QVector3D(-w, w, w), QVector2D(1.0, 1.0), QVector3D(0.0, 1.0, 0.0)) <<
		Vertex(QVector3D(-w, w, -w), QVector2D(1.0, 0.0), QVector3D(0.0, 1.0, 0.0)) <<

		// back
		Vertex(QVector3D(w, w, -w), QVector2D(0.0, 1.0), QVector3D(0.0, 0.0, -1.0)) <<
		Vertex(QVector3D(w, -w, -w), QVector2D(0.0, 0.0), QVector3D(0.0, 0.0, -1.0)) <<
		Vertex(QVector3D(-w, w, -w), QVector2D(1.0, 1.0), QVector3D(0.0, 0.0, -1.0)) <<
		Vertex(QVector3D(-w, -w, -w), QVector2D(1.0, 0.0), QVector3D(0.0, 0.0, -1.0)) <<

		// left
		Vertex(QVector3D(-w, w, w), QVector2D(0.0, 1.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(-w, w, -w), QVector2D(0.0, 0.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(-w, -w, w), QVector2D(1.0, 1.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(-w, -w, -w), QVector2D(1.0, 0.0), QVector3D(-1.0, 0.0, 0.0)) <<

		// bottom
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

void Widget::mousePressEvent(QMouseEvent* event)
{
	if (event->buttons() == Qt::LeftButton)
		mousePosition = QVector2D(event->localPos());
	event->accept();
}

void Widget::mouseMoveEvent(QMouseEvent* event)
{
	if (event->buttons() != Qt::LeftButton) return;

	QVector2D diff = QVector2D(event->localPos()) - mousePosition;
	mousePosition = QVector2D(event->localPos());

	float angle = diff.length();

	QVector3D axis = QVector3D(diff.y(), diff.x(), 0.0);

	rotation = QQuaternion::fromAxisAndAngle(axis, angle) * rotation;

	update();
}
