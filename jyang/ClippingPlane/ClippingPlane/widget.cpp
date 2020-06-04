#include "widget.h"

Widget::Widget(QWidget* parent): 
	QOpenGLWidget(parent)
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

	// enable blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	initShader();
	initCube(1.0, QColor(255, 0, 0, 255));
	objects[0]->translate(QVector3D(0.0, 0.0, 0.0));
	initCube(1.0, QColor(0, 255, 0, 255));
	objects[1]->translate(QVector3D(0.0, 0.0, 10.0));
	initPlane(2.0);
	//objects[2]->translate(QVector3D(0.0, 0.0, 0.0));
	//objects[2]->rotate(45, QVector3D(0.0, 1.0, 0.0));
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

	// view matrix
	vMatrix.setToIdentity();
	vMatrix.translate(0.0, 0.0, -5.0);
	vMatrix.rotate(30, 1.0, 0.0, 0.0);
	vMatrix.rotate(-30, 0.0, 1.0, 0.0);
	vMatrix.rotate(rotation);

	shaderProgram.bind();
	shaderProgram.setUniformValue("u_pMatrix", pMatrix);
	shaderProgram.setUniformValue("u_vMatrix", vMatrix);
	shaderProgram.setUniformValue("u_texture", 0);

	// draw objects
	for (int i = 0; i < objects.size(); i++) {
		objects[i]->draw(&shaderProgram, context()->functions());
	}
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

void Widget::initCube(float w, QColor color)
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

	QImage image(100, 100, QImage::Format_ARGB32_Premultiplied);
	image.fill(color);

	objects << new Object3D(vertices, indices, image);
}

void Widget::initPlane(float w)
{
	QVector<Vertex> vertices;
	vertices <<
		// vertical
		Vertex(QVector3D(0.0, w, w), QVector2D(0.0, 1.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(0.0, w, -w), QVector2D(0.0, 0.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(0.0, -w, w), QVector2D(1.0, 1.0), QVector3D(-1.0, 0.0, 0.0)) <<
		Vertex(QVector3D(0.0, -w, -w), QVector2D(1.0, 0.0), QVector3D(-1.0, 0.0, 0.0));

	QVector<GLuint> indices;
	indices << 0 << 1 << 2 << 2 << 1 << 3;
	for (int i = 0; i < 24; i += 4)
		indices << i + 0 << i + 1 << i + 2 << i + 2 << i + 1 << i + 3;

	QImage image(100, 100, QImage::Format_ARGB32_Premultiplied);
	image.fill(QColor(0, 0, 255, 200));

	objects << new Object3D(vertices, indices, image);
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
