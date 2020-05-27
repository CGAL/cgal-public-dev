uniform sampler2D qt_texture;
varying highp vec2 qt_texCoord;

void main(void) {
	gl_FragColor = texture2D(qt_texture, qt_texCoord);
}