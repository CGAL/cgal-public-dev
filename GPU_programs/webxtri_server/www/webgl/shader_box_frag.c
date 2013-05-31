
#ifdef GL_ES
precision highp float;
#endif

uniform sampler2D uSampler;

// varying vec4 vColor;
varying vec2 vTexCoord;

void main(void) {
    vec4 sample;
    sample = texture2D(uSampler, vec2(vTexCoord.s, vTexCoord.t));
    gl_FragColor = sample;
    //vec4(sample.x, sample.y, 0, 1.0);
}

