<html><head>
<meta http-equiv="content-type" content="text/html; charset=ISO-8859-1"></head><body>#ifdef GL_ES
    precision highp float;
    #endif
    
    uniform float uTime;
    uniform float uStep;
    uniform sampler2D uSampler;
    
    varying vec2 vTextureCoord;
      
    void main(void) {
        float fish = uTime;
        vec4 sample0,sample1,sample2,sample3;
        sample0 = texture2D(uSampler, vec2(vTextureCoord.s+uStep ,vTextureCoord.t+uStep));
        sample1 = texture2D(uSampler, vec2(vTextureCoord.s-uStep,vTextureCoord.t+uStep));
        sample2 = texture2D(uSampler, vec2(vTextureCoord.s-uStep,vTextureCoord.t-uStep));
        sample3 = texture2D(uSampler, vec2(vTextureCoord.s+uStep,vTextureCoord.t-uStep));
        
        gl_FragColor = (sample0 + sample1 + sample2 + sample3) /4.0;
        
        
    }</body></html>