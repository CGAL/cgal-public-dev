<html><head>
<meta http-equiv="content-type" content="text/html; charset=ISO-8859-1"></head><body>attribute vec3 aVertexPosition;
    attribute vec4 aVertexColor;
    attribute vec2 aTextureCoord;
    
    uniform mat4 uMVMatrix;
    uniform mat4 uPMatrix;
   
    varying vec2 vTextureCoord;
    varying vec4 vColor;
    
   
    void main(void) {
        
        gl_Position = uPMatrix * uMVMatrix * vec4(aVertexPosition, 1.0);
        vColor = aVertexColor;
        vTextureCoord = aTextureCoord;
    }</body></html>