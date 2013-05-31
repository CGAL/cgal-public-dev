<html><head>
<meta http-equiv="content-type" content="text/html; charset=ISO-8859-1"></head><body>#ifdef GL_ES
    precision highp float;
    #endif
    
    uniform sampler2D uSampler;
    
    varying vec4 vColor;
    varying vec2 vTextureCoord;
    
    void main(void) {
        vec4 sample;
        sample = texture2D(uSampler, vec2(vTextureCoord.s,vTextureCoord.t));
        gl_FragColor = sample;
        
    }</body></html>