
// fBoxShader:

#ifdef GL_ES
    precision highp float;
    #endif

    uniform sampler2D uSampler;

    varying vec4 vColor;
    varying vec2 vTextureCoord;

    void main(void) {
        vec4 sample;
        sample = texture2D(uSampler, vec2(vTextureCoord.s,vTextureCoord.t));
        gl_FragColor = sample;

    }

// vBoxShader:

attribute vec3 aVertexPosition;
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
    }

// fBoxShaderSil:

#ifdef GL_ES
    precision highp float;
    #endif

    uniform sampler2D uSampler;

    varying vec4 vColor;
    varying vec2 vTextureCoord;

    void main(void) {

        gl_FragColor = vec4(1.0,0.5,0.0,1.0);

    }

vBoxShaderSil:

attribute vec3 aVertexPosition;
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
    }

// fSurfaceShader:

#ifdef GL_ES
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


    }

// vSurfaceShader:

attribute vec3 aVertexPosition;
  attribute vec4 aVertexColor;
  attribute vec2 aTextureCoord;

  uniform mat4 uMVMatrix;
  uniform mat4 uPMatrix;
  varying vec2 vTextureCoord;

  void main(void) {

    gl_Position = uPMatrix * uMVMatrix * vec4(aVertexPosition, 1.0);
    vTextureCoord = aTextureCoord;
  }


uniform sampler2D tex1;
float blur=0.0 ;
vec2 TexelKernelh[13];
float BlurWeights[13];
void main()
{
vec4 baseColor = vec4(0.0, 0.0, 0.0, 0.0);

       float d = blur * 100.0 -300.0;
       int i= 0;

       TexelKernelh[0] = vec2( 0, -6 );
       TexelKernelh[1] = vec2( 0, -5 );
       TexelKernelh[2] = vec2( 0, -4 );
       TexelKernelh[3] = vec2( 0, -3 );
       TexelKernelh[4] = vec2( 0, -2 );
       TexelKernelh[5] = vec2( 0, -1 );
       TexelKernelh[6] = vec2( 0, 0 );
       TexelKernelh[7] = vec2( 0, 1 );
       TexelKernelh[8] = vec2( 0, 2 );
       TexelKernelh[9] = vec2( 0, 3 );
       TexelKernelh[10] = vec2( 0, 4 );
       TexelKernelh[11] = vec2( 0, 5 );
       TexelKernelh[12] = vec2( 0, 6 );

       BlurWeights[0] = 0.002216;
       BlurWeights[1] = 0.008764;
       BlurWeights[2] = 0.026995;
       BlurWeights[3] = 0.064759;
       BlurWeights[4] = 0.120985;
       BlurWeights[5] = 0.176033;
       BlurWeights[6] = 0.199471;
       BlurWeights[7] = 0.176033;
       BlurWeights[8] = 0.120985;
       BlurWeights[9] = 0.064759;
       BlurWeights[10] = 0.026995;
       BlurWeights[11] = 0.008764;
       BlurWeights[12] = 0.002216;




for (int i = 0; i < 13; i++)

   {

      baseColor += texture2D( tex1, gl_TexCoord[0].xy + TexelKernelh[i].yx /d) * BlurWeights[i];

   }

gl_FragColor = baseColor;
}







  