<html><head>
<meta http-equiv="content-type" content="text/html; charset=ISO-8859-1"></head><body>#ifdef GL_ES
    precision highp float;
    #endif
    
    uniform sampler2D uSampler;
    
    varying vec4 vColor;
    varying vec2 vTextureCoord;
    
    void main(void) {
        
        gl_FragColor = vec4(1.0,0.5,0.0,1.0);
        
    }</body></html>