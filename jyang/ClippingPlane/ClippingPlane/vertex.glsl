attribute highp vec4 qt_vertex;
attribute highp vec2 qt_texCoord0;
uniform highp mat4 qt_modelViewProjectionMatrix;

varying highp vec2 qt_texCoord;

void main(void) {
	gl_Position = qt_modelViewProjectionMatrix * qt_vertex;
	qt_texCoord = qt_texCoord0;
}